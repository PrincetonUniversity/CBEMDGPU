#include "system.h"
#include "potential.h"
#include "integrator.h"
#include "nvt.h"
#include <iostream>
#include "utils.h"
#include <omp.h>
#include <stdlib.h>

#ifndef NOGPU
#include "cudaHelper.h"
#endif

#include <math.h>

/*!
 * Invoke the program as 
 * $ ./md numThreads > log 2> err
 */ 
int main (int argc, char* argv[]) {
    if (argc !=2){
        // catch incorrect number of arguments
        printf("USAGE: %s <nthreads> \n",argv[0]);
        exit(1);
    }
    static int nthreads = atoi(argv[1]);
    omp_set_num_threads(nthreads);
    const int rngSeed = 3145;
    
    float Temp = 0.5;
    float timestep = 0.005;
	const int nAtoms = 50;
    
    systemDefinition a;
	a.setBox(12, 12, 12);
	a.setTemp(Temp);
    a.setMass(1.0);
	a.setRskin(1.0);
    a.setRcut(2.5);	// if slj needs to incorporate "delta" shift already so cell list is properly created
	a.initThermal(nAtoms, 1.01*Temp, rngSeed, 1.2);
	
    // select shifted lennard jones potential function
	#ifdef NVCC
	pointFunction_t pp = dev_slj;
	#else
	pointFunction_t pp = slj;
	#endif

	#ifdef NVCC
	systemProps cudaProps;
	cudaProps.displayAllProps();
	if (cudaProps.numDevices() > 0) {
		const int cudaDevID = 0;
		a.cudaThreads = cudaProps.maxThreadsPerBlock(cudaDevID);
		if (a.cudaThreads > 512) {
			a.cudaThreads = 512;
		}
		a.cudaBlocks = (int) ceil(a.numAtoms()/(1.0*a.cudaThreads));
		if (a.cudaBlocks < 1) a.cudaBlocks = 1;
	} else {
		return -1;
	}
	#endif

	a.setPotential(pp);
	std::vector <float> args(5);
	args[0] = 1.0; // epsilon
	args[1] = 1.0; // sigma
	args[2] = 0.0; // delta 
	args[3] = 0.0; // ushift
	a.setPotentialArgs(args);

    nvt_NH integrate (1.0);     // damping constant for thermostat = 1.0
    integrate.setTimestep(timestep);

    const int nSteps = 3000;
    const int report = nSteps/1000;

	for (unsigned int long step = 0; step < nSteps; ++step) {
		integrate.step(a);
		if (step%report == 0) {
			std::cout << step << "\t" << a.KinE() << "\t" << a.PotE() << "\t" << a.instantT() << "\t" << a.KinE() + a.PotE() << std::endl;
			a.writeSnapshot();
		}
	}

    return 0;
}
