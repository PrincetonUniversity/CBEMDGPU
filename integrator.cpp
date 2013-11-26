/*!
 * Integration (velocity) verlet step
 * \author Nathan A. Mahynski
 * \date 11/19/13
 */

#include "system.h"
#include "dataTypes.h"
#include <exception>
#include "common.h"
#include "integrator.h"
#include <omp.h>

/*!
 * Update the positions and velocities using velocity verlet integration.
 * This uses the accelerations stored on each atom.
 *
 * \param [in] sys System definition
 */ 
void integrator::verletStep_ (systemDefinition &sys) {
	static int start = 1;
	if (start) {
		try {
			lastAccelerations_.resize(sys.numAtoms());
		} catch (std::exception &e) {
			std::cerr << e.what() << std::endl;
			throw customException("Failed to initialize integrator due to memory constraints");
			return;
		}
		// calculate the forces initially (sets Up)
		calcForce_ (sys);
		
		// also get initial T
		float Uk = 0.0, tmp = 0.0;
		#pragma omp parallel
		{
			#pragma omp for shared(sys.atoms) schedule(dynamic, OMP_CHUNK) nowait
			for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
				sys.atoms[i].vel.x += dt_*0.5*(lastAccelerations_[i].x+sys.atoms[i].acc.x);
				sys.atoms[i].vel.y += dt_*0.5*(lastAccelerations_[i].y+sys.atoms[i].acc.y);
				sys.atoms[i].vel.z += dt_*0.5*(lastAccelerations_[i].z+sys.atoms[i].acc.z);
			}
			#pragma omp reduction(+:Uk) schedule(dynamic, OMP_CHUNK) nowait
			for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
				Uk += (sys.atoms[i].vel.x*sys.atoms[i].vel.x)+(sys.atoms[i].vel.y*sys.atoms[i].vel.y)+(sys.atoms[i].vel.z*sys.atoms[i].vel.z);
			}
		}
		Uk *= sys.mass();
		tmp = Uk;
		Uk *= 0.5;
		tmp /= (3.0*(sys.numAtoms()-1.0));
		sys.updateInstantTemp(tmp);
		sys.setKinE(Uk);
		start = 0;
	}
	
	// update positions based on current positions
	#pragma omp parallel
	{
		#pragma omp for shared(sys.atoms) schedule(dynamic, OMP_CHUNK) nowait
		for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
			sys.atoms[i].pos.x += dt_*(sys.atoms[i].vel.x+0.5*dt_*sys.atoms[i].acc.x);
			sys.atoms[i].pos.y += dt_*(sys.atoms[i].vel.y+0.5*dt_*sys.atoms[i].acc.y);
			sys.atoms[i].pos.z += dt_*(sys.atoms[i].vel.z+0.5*dt_*sys.atoms[i].acc.z);
			lastAccelerations_[i] = sys.atoms[i].acc;
		}	
	}

	// calculate new forces at new positions
	calcForce_ (sys);

	// update velocities and get Uk and kinetic temperature
	float tmp = 0.0, Uk = 0.0;
	#pragma omp parallel
	{
		#pragma omp for shared(sys.atoms) schedule(dynamic, OMP_CHUNK) nowait
		for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
			sys.atoms[i].vel.x += dt_*0.5*(lastAccelerations_[i].x+sys.atoms[i].acc.x);
			sys.atoms[i].vel.y += dt_*0.5*(lastAccelerations_[i].y+sys.atoms[i].acc.y);
			sys.atoms[i].vel.z += dt_*0.5*(lastAccelerations_[i].z+sys.atoms[i].acc.z);
		}
		#pragma omp reduction(+:Uk) schedule(dynamic, OMP_CHUNK) nowait
		for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
			Uk += (sys.atoms[i].vel.x*sys.atoms[i].vel.x)+(sys.atoms[i].vel.y*sys.atoms[i].vel.y)+(sys.atoms[i].vel.z*sys.atoms[i].vel.z);
		}
	}
	Uk *= sys.mass();
	tmp = Uk;
	Uk *= 0.5;
	tmp /= (3.0*(sys.numAtoms()-1.0));
	sys.updateInstantTemp(tmp);
	sys.setKinE(Uk);
}
