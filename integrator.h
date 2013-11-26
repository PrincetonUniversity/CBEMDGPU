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

class integrator {
	public:
		void setTimestep (const float dt) {dt_ = dt;}
		virtual void step (systemDefinition &sys) = 0;
	protected:
		cellList_cpu cl_;
		virtual void calcForce_ (systemDefinition &sys) = 0;
		void verletStep_ (systemDefinition &sys);
		std::vector <float3> lastAccelerations_;
		float dt_;
};


#endif