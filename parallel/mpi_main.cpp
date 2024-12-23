#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <mpi.h>
#include <chrono>
//FEATURE FLAGS
const bool PERIODIC_BOUNDARIES = false;
const bool ARTIFICAL_VISCOSITY = true;

//SYSTEM CONSTANTS
const double dt = 0.005;
const double length = 1.2;
const double totalTime = 0.2;
const double numParticles = 400;
const double dx = 0.6/80.0;
const double x0 = 0.5*length;
const double h = 0.016;
const double gam = 1.4;
const double m = 0.001875;
const double alpha_pi = 1.0;
const double beta_pi = 1.0;
const double phi = 0.1*h;


// manually counted @sanket, could be wrong double check once
const int FLOPS_PER_KERNEL = 15;          // Approximate FLOPs for kernel calculation
const int FLOPS_PER_GRAD_KERNEL = 20;     // Approximate FLOPs for gradient kernel
const int FLOPS_PER_UPDATE = 30;          // Basic update operations per particle pair
const int FLOPS_PER_VISCOSITY = 25;   

// struct to profile code
struct PerformanceMetrics {
    long long totalFlops;
    std::chrono::high_resolution_clock::time_point startTime;
    std::chrono::high_resolution_clock::time_point endTime;
    
    PerformanceMetrics() : totalFlops(0) {}
    
    void start() {
        startTime = std::chrono::high_resolution_clock::now();
    }
    
    void stop() {
        endTime = std::chrono::high_resolution_clock::now();
    }
    
    double getElapsedSeconds() {
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
        return duration.count() / 1e6;
    }
    
    double getMflops() {
        return (totalFlops / 1e6) / getElapsedSeconds();
    }
};


// each particle in the system
struct Particle {
    double x, v, a, rho, P, e, pwr;
};

// MPI datatype for Particle
MPI_Datatype MPI_PARTICLE;

void createParticleType() {
    const int nitems = 7;
    int blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[7] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
                            MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[7];
    // tell mpi where every field is relative to start of the memory block
    // macro defined as  offsetof(type, member)
    offsets[0] = offsetof(Particle, x);
    offsets[1] = offsetof(Particle, v);
    offsets[2] = offsetof(Particle, a);
    offsets[3] = offsetof(Particle, rho);
    offsets[4] = offsetof(Particle, P);
    offsets[5] = offsetof(Particle, e);
    offsets[6] = offsetof(Particle, pwr);
    // create MPI-TYPE
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_PARTICLE);
    MPI_Type_commit(&MPI_PARTICLE);
}


// initialize the system 
void init(std::vector<Particle>& localSystem, int rank, int size) {
    int particlesPerProcess = numParticles / size;
    int startIdx = rank * particlesPerProcess;
    int endIdx = (rank == size - 1) ? numParticles : (rank + 1) * particlesPerProcess;
    
    for(int i = startIdx; i < endIdx; i++) {
        Particle p;
        if(i < 320) {
            p.x = i*(dx/4.0)-x0+dx/4.0;
            p.v = 0;
            p.a = 0;
            p.rho = 1.0;
            p.P = 1.0;
            p.e = 2.5;
        } else {
            p.x = (i-320)*dx+0.5*dx;
            p.v = 0;
            p.a = 0;
            p.rho = 0.25;
            p.P = 0.1795;
            p.e = 1.795;
        }
        localSystem.push_back(p);
    }
}


// smoothing kernel
double kernel(double r) {
    double q = r/h;
    double w = 0;
    double sig = 2.0/(3.0*h);
    if(q>0.0 && q<=1.0){
        w = sig*(1.0-1.5*(q*q)+0.75*(q*q*q));
    }else if(q>1.0 && q<=2.0){
        w =  0.25*sig*(2-q)*(2-q)*(2-q);
    }
    return w;
}

double gradKernel(double r) {
    double sign = -r/fabs(r);
    r = fabs(r);
    double q = r/h;
    double gw = 0;
    double sig = 2.0/(3.0*h);
    if(q>0.0 && q<=1.0) {
        gw = sign*(2*q-1.5*(q*q))/(h*h);
    }else if(q>1.0 && q<=2.0){
        gw = sign*0.5*(q-2)*(q-2)/(h*h);
    }
    return gw;
}


// update the state of the local system
void update(std::vector<Particle>& localSystem, int rank, int size,PerformanceMetrics &metrics) {
    std::vector<Particle> allParticles;
    
    // Gather all particles to all processes
    int localSize = localSystem.size();
    std::vector<int> sizes(size), displacements(size);
    
    MPI_Allgather(&localSize, 1, MPI_INT, sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    displacements[0] = 0;
    for(int i = 1; i < size; i++) {
        displacements[i] = displacements[i-1] + sizes[i-1];
    }
    
    allParticles.resize(numParticles);
    MPI_Allgatherv(localSystem.data(), localSize, MPI_PARTICLE,
                   allParticles.data(), sizes.data(), displacements.data(),
                   MPI_PARTICLE, MPI_COMM_WORLD);
    
    // Calculate density and pressure
    for(auto& particle : localSystem) {
        particle.rho = m*2.0/(3.0*h);
        for(const auto& other : allParticles) {
            if(&other != &particle) {
                double r = fabs(particle.x - other.x);
                particle.rho += m*kernel(r);
                // count kernel flops
                metrics.totalFlops += FLOPS_PER_KERNEL;
            }
        }
        particle.P = (gam-1.0)*particle.rho*particle.e;
        // account for the floating point ops above
        metrics.totalFlops += 3;
    }
    
    // Update acceleration and power
    for(auto& particle : localSystem) {
        double a = 0;
        double pwr = 0;
        
        for(const auto& other : allParticles) {
            if(&other != &particle) {
                double r = particle.x - other.x;
                
                double Pi = particle.P;
                double rhoi = particle.rho;
                double Pj = other.P;
                double rhoj = other.rho;
                
                double delV = particle.v - other.v;
                double delP = (Pi/(rhoi*rhoi))+(Pj/(rhoj*rhoj));
                metrics.totalFlops += FLOPS_PER_UPDATE;
                if(ARTIFICAL_VISCOSITY) {
                    double phi_ij = (h*delV*r)/((fabs(r)*fabs(r))+(phi*phi));
                    double rhobar = 0.5*(rhoi+rhoj);
                    double cbar = 0.5*(sqrt(gam*(Pi/rhoi))+sqrt(gam*(Pj/rhoj)));
                    double Pi_ij = ((-alpha_pi*cbar*phi_ij+beta_pi*phi_ij*phi_ij)/rhobar)*(delV*r>0);
                    delP += Pi_ij;
                    metrics.totalFlops += FLOPS_PER_VISCOSITY;
                }
                
                a -= m*delP*gradKernel(r);
                pwr += 0.5*m*delP*delV*gradKernel(r);
                metrics.totalFlops += FLOPS_PER_GRAD_KERNEL;
            }
        }
        particle.a = a;
        particle.pwr = pwr;
    }
}

void record(const std::vector<Particle>& localSystem, int step, int rank, int size) {
    std::vector<Particle> allParticles;
    
    // Gather all particles to rank 0
    if(rank == 0) {
        allParticles.resize(numParticles);
    }
    
    int localSize = localSystem.size();
    std::vector<int> sizes(size), displacements(size);
    
    MPI_Gather(&localSize, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if(rank == 0) {
        displacements[0] = 0;
        for(int i = 1; i < size; i++) {
            displacements[i] = displacements[i-1] + sizes[i-1];
        }
    }
    
    MPI_Gatherv(localSystem.data(), localSize, MPI_PARTICLE,
                allParticles.data(), sizes.data(), displacements.data(),
                MPI_PARTICLE, 0, MPI_COMM_WORLD);
    
    if(rank == 0) {
        std::string filename = "./output/data_" + std::to_string(step) + ".dat";
        std::ofstream dat(filename);
        for(const auto& particle : allParticles) {
            dat << particle.x << '\t' << particle.e << '\t' 
                << particle.rho << '\t' << particle.P << '\t' 
                << particle.v << std::endl;
        }
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // initialize mflops counter struct
    PerformanceMetrics metrics;
    createParticleType();
    
    std::vector<Particle> localSystem;
    init(localSystem, rank, size);
    
    record(localSystem, 0, rank, size);
    if(rank == 0) std::cout << "0" << std::endl;
    
    metrics.start();  // Start timing
    
    for(int i = 1; i*dt <= totalTime; i++) {
        if(rank == 0) std::cout << i*dt << std::endl;
        
        // Leapfrog Integration
        for(auto& particle : localSystem) {
            // half-kick
            particle.v += 0.5*particle.a*dt;
            particle.e += 0.5*particle.pwr*dt;
            // drift
            particle.x += particle.v*dt;
            if(PERIODIC_BOUNDARIES) {
                if(particle.x < 0) particle.x += length;
                if(particle.x >= length) particle.x -= length;
            }
        }
        
        update(localSystem, rank, size,metrics);
        
        for(auto& particle : localSystem) {
            // second half-kick
            particle.v += 0.5*particle.a*dt;
            particle.e += 0.5*particle.pwr*dt;
        }
        
        update(localSystem, rank, size,metrics);
        record(localSystem, i, rank, size);
    }
    metrics.stop();  // Stop timing
    
    // Gather FLOPS from all processes
    long long totalFlops;
    MPI_Reduce(&metrics.totalFlops, &totalFlops, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if(rank == 0) {
        double totalTime = metrics.getElapsedSeconds();
        double mflops = (totalFlops / 1e6) / totalTime;
        std::cout << "Simulation completed in " << totalTime << " seconds" << std::endl;
        std::cout << "Total Performance: " << mflops << " MFLOPS" << std::endl;
        std::cout << "Per-process Performance: " << mflops/size << " MFLOPS" << std::endl;
    }
    MPI_Type_free(&MPI_PARTICLE);
    MPI_Finalize();
    return 0;
}
