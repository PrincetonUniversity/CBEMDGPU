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
    int tid; // thread number
    const int rngSeed = 3145;
    float Temp = 0.5; // this should be input to main
    float timestep = 0.005; // this should be input to main
    float tau2 = timestep*timestep;
    systemDefinition a;
    a.setBox(10, 10, 10);
    a.setTemp(Temp);
    a.setMass(1.0);
    a.setRskin(1.0);
    a.setRcut(1.0);
    a.initThermal(50, 1.01*Temp, rngSeed, 1.2);

    //pointFunction_t pp = pairUF;
    pointFunction_t pp = slj;
    a.setPotential(pp);
    nvt_NH integrate (1.0);
    integrate.setTimestep(timestep);

    const int nSteps = 3000;
    const int report = 10; //nSteps/1000;
#pragma omp parallel 
    {
	tid = omp_get_thread_num();
	if (tid == 0) {
	    std::cout <<  "# Number of threads = " << nthreads <<std::endl;
	}
	printf("Thread %d starting...\n",tid);
    }
	for (unsigned int long step = 0; step < nSteps; ++step) {
	    integrate.step2(a);
	    if (step%report == 0) {
		std::cout << step << "\t" << a.KinE() << "\t" << a.PotE() << "\t" << a.instantT() << "\t" << a.KinE() + a.PotE() << std::endl;
		a.writeSnapshot();
	    }
	}

    return 0;
}
