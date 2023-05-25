#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

//FEATURE FLAGS
const bool PERIODIC_BOUNDARIES = false;
const bool ARTIFICAL_VISCOSITY = true;

//SYSTEM CONSTANTS
const double dt = 0.005; //time step
const double length = 1.2; //length of tube
const double totalTime = 0.2;
const double numParticles = 400;
const double dx = 0.6/80.0; //spatial partition factor
const double x0 = 0.5*length;
const double h = 0.016; //smoothing length
const double gam = 1.4; //adiabatic index
const double m = 0.001875; //particle mass
const double alpha_pi = 1.0; //artificial viscosity constant
const double beta_pi = 1.0; //artificial viscosity constant
const double phi = 0.1*h; //prevents divergence in artificial viscosity

struct Particle{
	double x, v, a, rho, P, e, pwr;
	//Position, Velocity, Acceleration, Density, Pressure, Internal Energy, Internal Power
};

std::vector<Particle> System;

//Init Sod shock tube from G. R. Liu, M. B. Liu - Smoothed Particle Hydrodynamic
void init(){
	//Left side (packed side)
	for(int i=0;i<320;i++){
		System.emplace_back();
		System[i].x = i*(dx/4.0)-x0+dx/4.0;
		System[i].v = 0;
		System[i].a = 0;
		System[i].rho = 1.0;
		System[i].P = 1.0;
		System[i].e = 2.5;
	}
	//Right side (unpacked side)
	for(int i=320;i<numParticles;i++){
		System.emplace_back();
		System[i].x = (i-320)*dx+0.5*dx;
		System[i].v = 0;
		System[i].a = 0;
		System[i].rho = 0.25;
		System[i].P = 0.1795;
		System[i].e = 1.795;
	}
}

//Cubic-Spline Kernel normalized for one-dimension
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
 
//Gradient of the Cubic-Spline Kernel
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

void update(){
	//Calculate density and pressure for each particle
	for(int i=0;i<numParticles;i++){
		System[i].rho = m*2.0/(3.0*h);
		for(int j=0;j<numParticles;j++){
			if(i!=j){
				double r = System[i].x-System[j].x; //vector distance
				r=fabs(r); //scalar distance
				System[i].rho += m*kernel(r);
			}
		}
		System[i].P = (gam-1.0)*System[i].rho*System[i].e;
	}
	//Use pressure and density to calculate the time differential of velocity and energy with artificial viscosity
	for(int i=0;i<numParticles;i++){
		double a = 0;
		double pwr = 0;
		for(int j=0;j<numParticles;j++){
			if(j!=i){
				double r = System[i].x-System[j].x; //vector distance
				
				double Pi = System[i].P;     //pressure i
				double rhoi = System[i].rho; //density i
				double Pj = System[j].P;	 //pressure j
				double rhoj = System[j].rho; //density j

				double delV = System[i].v-System[j].v;           //velocity spatial derivative
				double delP = (Pi/(rhoi*rhoi))+(Pj/(rhoj*rhoj)); //pressure spatial derivative
				if(ARTIFICAL_VISCOSITY){
					double phi_ij = (h*delV*r)/((fabs(r)*fabs(r))+(phi*phi));
					double rhobar = 0.5*(rhoi+rhoj);
					double cbar = 0.5*(sqrt(gam*(Pi/rhoi))+sqrt(gam*(Pj/rhoj)));
					double Pi_ij = ((-alpha_pi*cbar*phi_ij+beta_pi*phi_ij*phi_ij)/rhobar)*(delV*r>0);
					delP += Pi_ij;
				}

				a -= m*delP*gradKernel(r);            //partial acceleration
				pwr += 0.5*m*delP*delV*gradKernel(r); //partial energy time derivative
			}
		}
		System[i].a = a;   //total acceleration
		System[i].pwr = pwr; //total energy time derivative
	}
}

//outputs variables: Position, Energy, Density, Pressure, Velocity
void record(int step){
	std::string filename = "data_" + std::to_string(step) + ".dat";
	std::ofstream dat;
	dat.open(filename);
	for(int j=0;j<numParticles;j++){
		dat << System[j].x << '\t' << System[j].e << '\t' << System[j].rho << '\t' << System[j].P << '\t' << System[j].v << std::endl;
	}
	dat.close();
}

int main(){
	init();
	record(0);
	std::cout << 0 << std::endl;
	for(int i=1;i*dt<=totalTime;i++){
		std::cout << i*dt << std::endl;
		//Leapfrog Integration
		for(int j=0;j<numParticles;j++){
			//half-kick
			System[j].v += 0.5*System[j].a*dt;
			System[j].e += 0.5*System[j].pwr*dt;
			//drift with PBCs
			System[j].x += System[j].v*dt;
			if(PERIODIC_BOUNDARIES){
				if(System[i].x < 0){System[i].x += length;}
				if(System[i].x >= length){System[i].x -= length;}
			}
		}
		update();
		for(int j=0;j<numParticles;j++){
			//second half-kick
			System[j].v += 0.5*System[j].a*dt;
			System[j].e += 0.5*System[j].pwr*dt;
		}
		update();
		record(i);
	}
	return 0;
}
