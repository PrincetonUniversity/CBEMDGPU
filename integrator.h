/*!
 * Integrator
 * \author Nathan A. Mahynski
 * \date 11/17/13
 */

#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

#include "system.h"
#include "cellList.h"
#include <vector>

//! Base class for integrators such as NVT (Nose-Hoover) or NVE ensembles
class integrator {
	public:
		void setTimestep (const float dt) {dt_ = dt;}   //!< Set the integrator timestep
        void calcForce (systemDefinition &sys); //!< Calculate the forces on each atom
		virtual void step (systemDefinition &sys) = 0; //!< Move the system forward a step in time
    
    protected:
		cellList_cpu cl_;   //!< Cell or neighbor list
		std::vector <float3> lastAccelerations_;    //!< Acceleration of particles on previous timestep (useful for NVE integrator)
		float dt_;      //!< Timestep size
		int start_;     //!< Flag for whether this object has been initialized or not
};


#endif
