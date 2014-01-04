/*!
 * Do NVT integration with Nose-Hoover thermostat.
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
 * Initialize integrator
 *
 * \param [in] Q Thermal mass for thermostat
 */ 
nvt_NH::nvt_NH (const float Q) {
    Q_ = Q; 
    gamma_ = 0.0;
	start_ = 1; 
}

/*!
 * Integrate a single timestep forward using Velocity-Verlet integration scheme.
 * This is NOT identical since some intermediate bookkeeping needs to be handled for thermostat.
 * Creates a cell list the first time it is called.
 * \param [in, out] sys System definition
 */
void nvt_NH::step (systemDefinition &sys) {
    int chunk = OMP_CHUNK;
    if (start_) {
	try {
		cellList_cpu tmpCL (sys.box(), sys.rcut(), sys.rskin());
		cl_ = tmpCL;
	} catch (std::exception &e) {
	    std::cerr << e.what() << std:: endl;
	    throw customException("Failed to integrate on first step");
	}

	gamma_ = 0.0;
	for (unsigned int i = 0; i < sys.numAtoms(); ++i)  {
	    gamma_ += (sys.atoms[i].vel.x*sys.atoms[i].vel.x)+(sys.atoms[i].vel.y*sys.atoms[i].vel.y)+(sys.atoms[i].vel.z*sys.atoms[i].vel.z);
	}
	gamma_ -= (3.0*(sys.numAtoms()-1.0))*sys.instantT();
	gamma_ /= Q_;

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
	gammadot_ = 0.0;
	gammadd_ = 0.0;

    }
    tau2_ = Q_/ ((3.0*(sys.numAtoms()-1.0))*sys.targetT());
    gammadd_ = 1/tau2_*(sys.instantT()/sys.targetT()-1);
    // (1) update thermostat velocity and thermostat position
    // velocity half-step
    gammadot_ += dt_*0.5*gammadd_;
    // position step
    gamma_ += gammadot_*dt_;
    // (2) evolve particle velocities
#pragma omp parallel shared(sys)
    {
#pragma omp for schedule(dynamic,OMP_CHUNK)
    for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
	sys.atoms[i].vel.x = sys.atoms[i].vel.x*exp(-gammadot_*dt_*0.5) + 0.5*dt_*(sys.atoms[i].acc.x);
	sys.atoms[i].vel.y = sys.atoms[i].vel.y*exp(-gammadot_*dt_*0.5) + 0.5*dt_*(sys.atoms[i].acc.y);
	sys.atoms[i].vel.z = sys.atoms[i].vel.z*exp(-gammadot_*dt_*0.5) + 0.5*dt_*(sys.atoms[i].acc.z);
    }
    // (3) evolve particle positions
#pragma omp for schedule(dynamic,OMP_CHUNK)
    for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
	sys.atoms[i].pos.x += sys.atoms[i].vel.x*dt_;
	sys.atoms[i].pos.y += sys.atoms[i].vel.y*dt_;
	sys.atoms[i].pos.z += sys.atoms[i].vel.z*dt_;
    }
    // (4) calc force
    }
    calcForce(sys);
    // (5) evolve particle velocities
#pragma omp parallel shared(sys)
    {
#pragma omp for schedule(dynamic,OMP_CHUNK)
    for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
	sys.atoms[i].vel.x = (sys.atoms[i].vel.x+sys.atoms[i].acc.x*dt_*0.5)*exp(-gammadot_*dt_*0.5);
	sys.atoms[i].vel.y = (sys.atoms[i].vel.y+sys.atoms[i].acc.y*dt_*0.5)*exp(-gammadot_*dt_*0.5);
	sys.atoms[i].vel.z = (sys.atoms[i].vel.z+sys.atoms[i].acc.z*dt_*0.5)*exp(-gammadot_*dt_*0.5);
    }
    }
    float Uk = 0.0;
    float tmp = 0.0;
#pragma omp parallel for reduction(+:Uk)
    // get temperature and kinetic energy
    for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
	Uk += (sys.atoms[i].vel.x*sys.atoms[i].vel.x)+(sys.atoms[i].vel.y*sys.atoms[i].vel.y)+(sys.atoms[i].vel.z*sys.atoms[i].vel.z);
    }
    Uk *= sys.mass();
    tmp = Uk;
    Uk *= 0.5;
    tmp /= (3.0*(sys.numAtoms()-1.0));
    sys.updateInstantTemp(tmp);
    sys.setKinE(Uk);
    // std::cout << "KINETIC ENERGY IS " << Uk << std::endl;
    // (6) update thermostat velocity
    gammadd_ = 1/tau2_*(sys.instantT()/sys.targetT()-1);
    gammadot_ += dt_*0.5*gammadd_;

}

