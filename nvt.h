/*!
 * NVT integration with Nose-Hoover thermostat.
 * \date 11/18/13
 */

#ifndef __NVT_NH_H__
#define __NVT_NH_H__

#include "system.h"
#include "integrator.h"

//! Uses Nose-Hover integration method to thermostat a system (constant T rather than E)
class nvt_NH : public integrator {
public:
    nvt_NH (const float Q);
    ~nvt_NH () {}
    void step (systemDefinition &sys);
private:
    float Q_;           //!< Thermostat's 'mass'
    float gamma_;       //!< Thermostat 'position' (it is essentially a spring)
    float tau2_;        //!< Square of damping constant, tau
    float gammadot_;    //!< First derivative of gamma (thermostat 'velocity')
    float gammadd_; //!< Second derivative of gamma (thermostat 'acceleration')
};

#endif

