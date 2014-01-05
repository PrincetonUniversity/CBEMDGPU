#include "system.h"
#include "potential.h"
#include "integrator.h"
#include "nvt.h"
#include <iostream>
#include "utils.h"
#include <omp.h>
#include <stdlib.h>
#include <math.h>

#ifndef NOGPU
#include "cudaHelper.h"
#endif

/*!
 * Invoke the program as 
 * $ ./timing numThreads nAtoms rs nsteps 
 */ 
int main (int argc, char* argv[]) {

	double t1, t2;
	t1 = omp_get_wtime();

    	if (argc !=5){
		// catch incorrect number of arguments
		printf("USAGE: %s <nthreads> <natoms> <rs> <nsteps> \n",argv[0]);
		exit(1);
    	}

    	static int nthreads = atoi(argv[1]);
    	const int nAtoms = atoi(argv[2]);
	const double rs = atof(argv[3]);
    	const int nSteps = atoi(argv[4]);

	const double L = pow(2.0*nAtoms, 1.0/3.0);
	const double rCut = 2.5;

	omp_set_num_threads(nthreads);
    	const int rngSeed = 3145;
    	float Temp = 0.5; 
    	float timestep = 0.005;
	
    	systemDefinition a;
    	a.setBox(L, L, L);
	a.setTemp(Temp);
    	a.setMass(1.0); 
	a.setRskin(rs);
    	a.setRcut(rCut);	// if slj needs to incorporate "delta" shift already so cell list is properly created
   	a.initThermal(nAtoms, Temp, rngSeed, 1.2);

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
		a.cudaBlocks = (int) ceil(a.numAtoms()/a.cudaThreads);
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

    	nvt_NH integrate (1.0);
    	integrate.setTimestep(timestep);


	for (unsigned int long step = 0; step < nSteps; ++step) {
		integrate.step(a);
	}
	
	t2 = omp_get_wtime();
	double time_diff = t2 - t1;

	std::cout << nthreads << " " << nAtoms << " " << rs << " " << nSteps << " " << time_diff << std::endl;

	return 0;
}
