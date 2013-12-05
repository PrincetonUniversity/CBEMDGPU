#include "system.h"
#include "potential.h"
#include "integrator.h"
#include "nvt.h"
#include "potential.h"
#include <iostream>
#include "utils.h"
#include <omp.h>
#include <stdlib.h>

/*!
 * Invoke the program as 
 * $ ./md numThreads > log 2> err
 */ 
int main (int argc, char* argv[]) {
	omp_set_num_threads(atoi(argv[1]));
		
	const int rngSeed = 1024;
	systemDefinition a;
	a.setBox(10, 10, 10);
	a.setTemp(1.0);
	a.setMass(1.0);
	a.setRskin(1.0);
	a.setRcut(1.0);
	a.initRandom(5000, rngSeed);

	pointFunction_t pp = pairUF;
	a.setPotential(pp);

	nvt_NH integrate (1000.0);
	integrate.setTimestep(0.005);

	const int nSteps = 1000000;
	const int report = 10; //nSteps/1000;
	for (unsigned int long step = 0; step < nSteps; ++step) {
		integrate.step(a);
		if (step%report == 0) {
			std::cout << step << "\t" << a.KinE() << "\t" << a.PotE() << "\t" << a.instantT() << std::endl;
			a.writeSnapshot();
		}
	}

	return 0;
}
