#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

//SYSTEM CONSTANTS
const double dt = 0.005;
const double length = 1.2;
const double totalTime = 0.2;
const double numParticles = 400;
const double dx = 0.6/80.0;
const double x0 = 0.6;
const double h = 0.032; //smoothing length
const double gam = 1.4;
const double m = 0.001875;

struct Particle{
	double m, x, v, vold, a, rho, P, e, eold, H;
};

std::vector<Particle> System;

//Init System sin wave
/*
void init(){
	//Change this function for desired Initial Conditions
	for(int i=0;i<numParticles;i++){
		System.emplace_back();
		double x = i*(length/numParticles);
		System[i].m = m;
		System[i].x = x-sin(M_PI*x)*(length/numParticles);
		System[i].v = 0;
		System[i].a = 0;
		System[i].rho = 0;
		System[i].P = 0;
	}
	//System[0].x += 0.5*(length/numParticles);
}
*/

//Init Sod shock tube from G. R. Liu, M. B. Liu - Smoothed Particle Hydrodynamic
void init(){
	//Left side (packed side)
	for(int i=0;i<320;i++){
		System.emplace_back();
		double x = i*(dx/4.0)-x0+dx/4.0;
		System[i].m = m;
		System[i].x = x;
		System[i].v = 0;
		System[i].a = 0;
		System[i].rho = 1.0;
		System[i].P = 1.0;
		System[i].e = 2.5;
	}
	for(int i=320;i<numParticles;i++){
		System.emplace_back();
		double x = (i-320)*dx+0.5*dx;
		System[i].m = m;
		System[i].x = x;
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
    //std::cout << q << "q" <<std::endl;
    double w = 0;
    double sig = 2.0/(3.0*h);
    if (q>0.0 && q<=1.0) {
        w = sig*(1.0-1.5*(q*q)+0.75*(q*q*q));
        //std::cout << w << std::endl;
    } else if (q>1.0 && q<=2.0) {
        w =  0.25*sig*(2-q)*(2-q)*(2-q);
    }
    //std::cout << w << std::endl;
    return w;
}
 
//Gradient of the Cubic-Spline Kernel
double gradKernel(double r) {
	double sign = r/fabs(r);
	r = fabs(r);
    double q = r/h;
    double gw = 0;
    double sig = 2.0/(3.0*h);
    if (q>0.0 && q<=1.0) {
        gw = sign*(2*q-1.5*(q*q))/(h*h);
    } else if (q>1.0 && q<=2.0) {
        gw = sign*0.5*(q-2)*(q-2)/(h*h);
    }
    return gw;
}


//Calculate the densities and pressure for each particle
void calculateRhoPressure(){
	for(int i=0;i<numParticles;i++){
		System[i].rho = m*2.0/(3.0*h);
		for(int j=0;j<numParticles;j++){
			if(i!=j){
				double r = System[i].x-System[j].x;
				r=fabs(r);
				System[i].rho += System[j].m*kernel(r);
			}
		}
		System[i].P = (gam-1.0)*System[i].rho*System[i].e;
		//std::cout << System[i].P << 'p' << std::endl;
		
	}
}

void calculateAcc(){
	for(int i=0;i<numParticles;i++){
		double a = 0;
		double h = 0;
		for(int j=0;j<numParticles;j++){
			if (j!=i){
				double r = System[i].x-System[j].x;
				double Pi = System[i].P;
				double rhoi = System[i].rho;
				double Pj = System[j].P;
				double rhoj = System[j].rho;
				double delp = m*((Pi/(rhoi*rhoi))+(Pj/(rhoj*rhoj)));
				a += delp*gradKernel(r);
				h += 0.5*delp*(System[i].v-System[j].v)*gradKernel(r);
			}
		}
		System[i].a = a;
		System[i].H = h;
		//std::cout << System[i].a<< 'a' << std::endl;
	}
}

int main(){
	std::ofstream xyz;
	xyz.open("SPH.xyz");
	std::ofstream dat;
	dat.open("density.dat");
	init();
	xyz << numParticles << std::endl;
	xyz << "test" << std::endl;
	/*
	for(int i=0;i<30;i++){
		double delx = 2*h/30.0;
	 	double r = i*delx;
	 	std::cout <<r << '\t'<< gradKernel(r) << std::endl;
	}
	return 0;
	*/
	for(int i=0;i<numParticles;i++){
		xyz << 1 << '\t' << System[i].x << '\t' << 0 << '\t' << 0 << std::endl;
	}		
	std::cout << 0 << std::endl;
	//return 0;
	int nskip = 20;
	for(int i=1;i*dt<=totalTime;i++){
		std::cout << i*dt << std::endl;
		if(i>1){
			for(int j=nskip;j<numParticles-nskip;j++){
				System[j].vold = System[j].v;
				System[j].eold = System[j].e;
				System[j].v += 0.5*System[j].a*dt;
				System[j].e += 0.5*System[j].H*dt;
			}
		}
		calculateRhoPressure();
		calculateAcc();
		for(int j=nskip;j<numParticles-nskip;j++){
			if(i==1){
				//half-kick
				System[j].v += 0.5*System[j].a*dt;
				System[j].e += 0.5*System[j].H*dt;
				//drift
				System[j].x += System[j].v*dt;
			}else{
				System[j].v = System[j].vold + System[j].a*dt;
				System[j].e = System[j].eold + System[j].e*dt;
				//drift
				System[j].x += System[j].v*dt;
			}
		}
		xyz << numParticles << std::endl;
		xyz << "test" << std::endl;
		/*
		for(int j=0;j<numParticles;j++){
			xyz << 1 << '\t' << System[j].x << '\t' << 0 << '\t' << 0 << std::endl;
		}
		*/		
	}
	for(int j=0;j<numParticles;j++){
			dat << System[j].x << '\t' << System[j].rho << std::endl;
	}
	return 0;
}