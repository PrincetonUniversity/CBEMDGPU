#include "system.h"
#include "potential.h"
#include "integrator.h"
#include "nve.h"
#include <iostream>
#include "utils.h"
#include <omp.h>
#include <stdlib.h>
#include <math.h>

/*!
 * Invoke the program as 
 * $ ./test_nve numThreads 
 */

// this function runs 400 atoms at T=0.71. 
// we use these results verify that total energy is conserved with NVE
int main (int argc, char* argv[]) {
    	if (argc !=2){
		// catch incorrect number of arguments
		printf("USAGE: %s <nthreads> \n",argv[0]);
		exit(1);
    	}
    	static int nthreads = atoi(argv[1]);
    	omp_set_num_threads(nthreads);
    	const int rngSeed = 3145;
    	float Temp = 0.71; 
    	float timestep = 0.005; 
	const int nAtoms = 400;
	const double L = 16.796;
    	systemDefinition a;
    	a.setBox(L, L, L);
	a.setTemp(Temp);
    	a.setMass(1.0); 
	a.setRskin(0);
    	a.setRcut(2.5);	// if slj needs to incorporate "delta" shift already so cell list is properly created
   	a.initThermal(nAtoms, Temp, rngSeed, 2.0);

	pointFunction_t pp = slj;

	a.setPotential(pp);
	std::vector <float> args(5);
	args[0] = 1.0; // epsilon
	args[1] = 1.0; // sigma
	args[2] = 0.0; // delta 
	args[3] = 0.0; // ushift
	a.setPotentialArgs(args);

    	nve integrate;
	integrate.setTimestep(timestep);

    	const int nSteps = 10000;
    	const int report = 1; 

	for (unsigned int long step = 0; step < nSteps; ++step) {
		integrate.step(a);
		if (step%report == 0) {
			//std::cout << a.numAtoms() << "\t" << step << "\t" << a.KinE() << "\t" << a.PotE() << "\t" << a.instantT() << "\t" << a.KinE() + a.PotE() << std::endl;
			printf("%u \t %2.2f \t %2.2f \t %2.4f \t %2.2f \n", step, a.KinE(), a.PotE(), a.instantT(), a.KinE()+a.PotE());
			//a.writeSnapshot();
		}
	}

    return 0;
}
