#include "system.h"
#include "potential.h"
#include "integrator.h"
#include "nvt.h"
#include <iostream>
#include "utils.h"
#include <omp.h>
#include <stdlib.h>


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
    	float Temp = 0.5; // this should be input to main
    	float timestep = 0.005; // this should be input to main
    	systemDefinition a;
    	a.setBox(12, 12, 12);
	a.setTemp(Temp);
    	a.setMass(1.0); 
	a.setRskin(1.0);
    	a.setRcut(2.5);	// if slj needs to incorporate "delta" shift already so cell list is properly created
   	a.initThermal(50, 1.01*Temp, rngSeed, 1.2);

    	/*pointFunction_t pp = pairUF;
 	std::vector <float> args(1);
	args[0] = 1.0; // epsilon
	*/
	pointFunction_t pp = slj;
	a.setPotential(pp);
	std::vector <float> args(5);
	args[0] = 1.0; // epsilon
	args[1] = 1.0; // sigma
	args[2] = 0.0; // delta 
	args[3] = 0.0; // ushift
	a.setPotentialArgs(args);

    	nvt_NH integrate (1.0);
    	integrate.setTimestep(timestep);

    	const int nSteps = 3000;
    	const int report = 10; //nSteps/1000;

	for (unsigned int long step = 0; step < nSteps; ++step) {
	    integrate.step2(a);
	    if (step%report == 0) {
		std::cout << step << "\t" << a.KinE() << "\t" << a.PotE() << "\t" << a.instantT() << "\t" << a.KinE() + a.PotE() << std::endl;
		a.writeSnapshot();
	    }
	}

    return 0;
}
