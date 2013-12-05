/*!
 * Do NVE integration
 * \author Nathan A. Mahynski
 * \date 11/18/13
 */

#include "system.h"
#include "nvt.h"
#include "cellList.h"
#include <exception>
#include "common.h"
#include <math.h>
#include <vector>
#include <omp.h>

/*!
 * Integrate a single timestep forward using Velocity-Verlet integration scheme.
 * Update the positions and velocities using velocity verlet integration.
 * This uses the accelerations stored on each atom.
 * Creates a cell list the first time it is called.
 * \param [in, out] sys System definition
 */  
void nve::step (systemDefinition &sys) {
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
		calcForce (sys);
		
		// also get initial T
        #pragma omp parallel
		{
            #pragma omp for shared(sys.atoms) schedule(dynamic, OMP_CHUNK) 
			for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
				sys.atoms[i].vel.x += dt_*0.5*(lastAccelerations_[i].x+sys.atoms[i].acc.x);
				sys.atoms[i].vel.y += dt_*0.5*(lastAccelerations_[i].y+sys.atoms[i].acc.y);
				sys.atoms[i].vel.z += dt_*0.5*(lastAccelerations_[i].z+sys.atoms[i].acc.z);
			}
        }
        
        float tmp = 0.0, Uk = 0.0;
        #pragma omp parallel
		{
            #pragma omp reduction(+:Uk) schedule(dynamic, OMP_CHUNK) 
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
        #pragma omp for shared(sys.atoms) schedule(dynamic, OMP_CHUNK) 
		for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
			sys.atoms[i].pos.x += dt_*(sys.atoms[i].vel.x+0.5*dt_*sys.atoms[i].acc.x);
			sys.atoms[i].pos.y += dt_*(sys.atoms[i].vel.y+0.5*dt_*sys.atoms[i].acc.y);
			sys.atoms[i].pos.z += dt_*(sys.atoms[i].vel.z+0.5*dt_*sys.atoms[i].acc.z);
			lastAccelerations_[i] = sys.atoms[i].acc;
		}
	}
    
	// calculate new forces at new positions
	calcForce (sys);
    
	// update velocities and get Uk and kinetic temperature
    #pragma omp parallel
	{
        #pragma omp for shared(sys.atoms) schedule(dynamic, OMP_CHUNK) 
		for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
			sys.atoms[i].vel.x += dt_*0.5*(lastAccelerations_[i].x+sys.atoms[i].acc.x);
			sys.atoms[i].vel.y += dt_*0.5*(lastAccelerations_[i].y+sys.atoms[i].acc.y);
			sys.atoms[i].vel.z += dt_*0.5*(lastAccelerations_[i].z+sys.atoms[i].acc.z);
		}
    }
    
    // get temperature and kinetic energy
    float tmp = 0.0, Uk = 0.0;
    #pragma omp parallel
	{
        #pragma omp reduction(+:Uk) schedule(dynamic, OMP_CHUNK) 
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


