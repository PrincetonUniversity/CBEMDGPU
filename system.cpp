/*!
 * System methods
 * \author Nathan A. Mahynski
 * \date 11/17/13
 */

#include "system.h"
#include "common.h"
#include <stdlib.h>
#include "potential.h"
#include <math.h>
#include <boost/tr1/random.hpp>

#define RNG (1.0*rand())/RAND_MAX

/*!
 * Takes a __device__ potential function and sets the host equivalent so that the host kernal can
 * pass the fiunction to the device if using CUDA, else just sets the "host" function.
 */ 
void systemDefinition::setPotential (pointFunction_t pp) {
#ifdef NVCC
	dev_potential = pp;
	cudaMemcpyFromSymbol(&potential, dev_potential, sizeof(pointFunction_t));
#else
	potential = pp;
#endif
}

/*!
 * Initialize a system of N atoms with random velocities and positions.
 * Net momeentum is automatically initialized to zero.
 *
 * \param [in] N Number of atoms to create
 * \param [in] rngSeed Random number generator seed
 */
void systemDefinition::initRandom (const int N, const int rngSeed) {
	float3 totMomentum;
	totMomentum.x = 0; totMomentum.y = 0; totMomentum.z = 0;

	if (N < 1) {
		throw customException ("N must be > 0");
		return;
	}

	srand(rngSeed);
	int chunk = OMP_CHUNK;
	unsigned int i;
	atoms.resize(N);
#pragma omp parallel private(i)
	{
#pragma	omp for schedule(dynamic, chunk)
	for (i = 0; i < N; ++i) {
		if (i < N-1) {
			atoms[i].vel.x = (RNG-0.5);
			atoms[i].vel.y = (RNG-0.5);
			atoms[i].vel.z = (RNG-0.5);
			totMomentum.x += atoms[i].vel.x;
			totMomentum.y += atoms[i].vel.y;
			totMomentum.z += atoms[i].vel.z;
		} else {
			atoms[i].vel.x = -totMomentum.x;
			atoms[i].vel.y = -totMomentum.y;
			atoms[i].vel.z = -totMomentum.z;
		}
		atoms[i].pos.x = (RNG)*box_.x;
		atoms[i].pos.y = (RNG)*box_.y;
		atoms[i].pos.z = (RNG)*box_.z;
		atoms[i].acc.x = 0;
		atoms[i].acc.y = 0;
		atoms[i].acc.z = 0;
	}
	}
}

/*!
* Initialize a system of N atoms with random velocities and positions.
* Net momeentum is automatically initialized to zero.
* 
* \param [in] N Number of atoms to create
* \param [in] rngSeed Random number generator seed
*/
void systemDefinition::initThermal (const int N, const float Tset, const int rngSeed) {
    float3 totMomentum;
    totMomentum.x = 0; totMomentum.y = 0; totMomentum.z = 0;

    if (N < 1) {
	throw customException ("N must be > 0");
	return;
    }

    srand(rngSeed);
    atoms.resize(N);

    // maxwell boltzmann distribution has mean 0 stdev kT/m in each dimension
    typedef std::tr1::linear_congruential<int, 16807, 0, (int)((1U << 31) -1 ) > Myceng;
    Myceng eng;
    float sig = sqrt(Tset/mass_);
    std::tr1::normal_distribution<float> distribution(0.0,sig);
    float rannum = 0.0;
    float tmpT = 0.0;
    for (unsigned int i = 0; i < N; ++i) {
	atoms[i].pos.x = (RNG)*box_.x;
	atoms[i].pos.y = (RNG)*box_.y;
	atoms[i].pos.z = (RNG)*box_.z;
	atoms[i].acc.x = 0;
	atoms[i].acc.y = 0;
	atoms[i].acc.z = 0;
	if (i < N-1) {
	    atoms[i].vel.x = distribution(eng);
	    atoms[i].vel.y = distribution(eng);
	    atoms[i].vel.z = distribution(eng);
	    totMomentum.x += atoms[i].vel.x;
	    totMomentum.y += atoms[i].vel.y;
	    totMomentum.z += atoms[i].vel.z;
	} else {
	    atoms[i].vel.x = -totMomentum.x;
	    atoms[i].vel.y = -totMomentum.y;
	    atoms[i].vel.z = -totMomentum.z;
	}
	tmpT += (atoms[i].vel.x*atoms[i].vel.x + atoms[i].vel.y*atoms[i].vel.y + atoms[i].vel.z*atoms[i].vel.z);
    }

    // do velocity rescaling to get exactly the right T
    tmpT *= mass_/(3.0*(N-1));
    tmpT /= Tset;
    tmpT = 1.0/tmpT;
    for (unsigned int i = 0; i < N; i ++ ) {
	atoms[i].vel.x = atoms[i].vel.x*sqrt(tmpT);
	atoms[i].vel.y = atoms[i].vel.y*sqrt(tmpT);
	atoms[i].vel.z = atoms[i].vel.z*sqrt(tmpT);
    }
}

/*!
 * Write instantaneous snapshot of the system to a file called "trajectory.xyz"
 * This file is appended not overwritten consecutively.
 */
void systemDefinition::writeSnapshot () {
	static int snapNum = 0;
	if (snapFile_ == NULL) {
		snapFile_ = fopen("trajectory.xyz", "w");
	}
	
	fprintf(snapFile_, "%d\nSnapshot #%d\n", atoms.size(), snapNum);

	for (unsigned int i = 0; i < atoms.size(); ++i) {
		fprintf(snapFile_, "%s\t%g\t%g\t%g\n", "A", atoms[i].pos.x, atoms[i].pos.y, atoms[i].pos.z);
	}

	snapNum++;
}
