/*!
 * Do NVE integration
 * \author Nathan A. Mahynski
 * \date 11/18/13
 */

#include "system.h"
#include "nve.h"
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
    int chunk = OMP_CHUNK;
    if (start_) {
        try {
            cellList_cpu tmpCL (sys.box(), sys.rcut(), sys.rskin());
            cl_ = tmpCL;
        } catch (std::exception &e) {
            std::cerr << e.what() << std:: endl;
            throw customException("Failed to integrate on first step");
        }

        // get initial temperature
        calcForce(sys);
        float tmp=0.0, Uk = 0.0;
        #pragma omp parallel for reduction(+:Uk)
        for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
            Uk += (sys.atoms[i].vel.x*sys.atoms[i].vel.x)+(sys.atoms[i].vel.y*sys.atoms[i].vel.y)+(sys.atoms[i].vel.z*sys.atoms[i].vel.z);
        }
        Uk *= sys.mass();
        tmp = Uk;
        Uk *= 0.5;
        tmp /= (3.0*(sys.numAtoms()-1.0));
        sys.updateInstantTemp(tmp);
        sys.setKinE(Uk);
        start_ = 0;
    }
    
    // (1) evolve particle velocities
    #pragma omp parallel shared(sys)
    {
        #pragma omp for schedule(dynamic,OMP_CHUNK)
        for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
            sys.atoms[i].vel.x += 0.5*dt_*(sys.atoms[i].acc.x);
            sys.atoms[i].vel.y += 0.5*dt_*(sys.atoms[i].acc.y);
            sys.atoms[i].vel.z += 0.5*dt_*(sys.atoms[i].acc.z);
        }

        // (2) evolve particle positions
        #pragma omp for schedule(dynamic,OMP_CHUNK)
        for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
            sys.atoms[i].pos.x += sys.atoms[i].vel.x*dt_;
            sys.atoms[i].pos.y += sys.atoms[i].vel.y*dt_;
            sys.atoms[i].pos.z += sys.atoms[i].vel.z*dt_;
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


