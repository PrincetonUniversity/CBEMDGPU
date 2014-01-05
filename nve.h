/*!
 * NVE integration.
 * \date 11/18/13
 */

#ifndef __NVE_H__
#define __NVE_H__

#include "system.h"
#include "integrator.h"

//! Integration scheme that preserves total energy of the system
class nve : public integrator {
public:
    nve () {start_ = 1;}
    ~nve () {}
    void step (systemDefinition &sys);
};

#endif
