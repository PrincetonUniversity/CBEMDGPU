/*!
 * NVT integration with Nose-Hoover thermostat.
 * \author Nathan A. Mahynski
 * \date 11/18/13
 */

#ifndef __NVT_NH_H__
#define __NVT_NH_H__

#include "system.h"
#include "integrator.h"

class nvt_NH : public integrator {
public:
    nvt_NH (const float Q);
    ~nvt_NH () {}
    void step (systemDefinition &sys);
private:
    float Q_, gamma_, tau2_, gammadot_, gammadd_;
};

#endif

