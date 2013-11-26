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
}

/*!
 * Integrate a single timestep forward using Velocity-Verlet integration
 * Creates a cell list the first time it is called.
 * \param [in, out] sys System definition
 */  
void nvt_NH::step (systemDefinition &sys) {
	static int start = 1;
	if (start) {
		try {
			cellList_cpu tmpCL (sys.box(), sys.rcut(), sys.rskin());
			cl_ = tmpCL;
		} catch (std::exception &e) {
			std::cerr << e.what() << std::endl;
			throw customException("Failed to integrate on first step");
		}
		start = 0;
	}
	verletStep_ (sys);
}

/*!
 * Calculate the pairwise forces in a system.  This also calculates the potential energy of a system.
 * The kinetic energy is calculated during the verlet integration.
 *
 * \param [in, out] sys System definition
 */ 
void nvt_NH::calcForce_ (systemDefinition &sys) {
	// For cache coeherency allocate new space for calculations
	float3 empty;
	empty.x = 0; empty.y = 0; empty.z = 0;
	std::vector <float3> acc (sys.numAtoms(), empty);	
	const float rc = sys.rcut();

	// every time, check if the cell list needs to be updated first
	cl_.checkUpdate(sys);

	// Find pairwise forces to get accelerations (cell lists)
	// Get Up and Uk at this time
	// Keep loop "linear" so OMP can handle this loop best
	float Up = 0.0;
	const float3 box = sys.box();
	const float invMass = 1.0/sys.mass();
	
	// brute force for comparison with cell lists (make sure to comment out checkUpdate above too when doing speed calc)
	/*for (int i = 0; i < sys.numAtoms(); ++i) {
		for (int j = i+1; j < sys.numAtoms(); ++j) {
			float3 pf;
			Up += sys.potential (&sys.atoms[i].pos, &sys.atoms[j].pos, &pf, &box, &rc);
			acc[i].x -= pf.x*invMass;
			acc[i].y -= pf.y*invMass;
			acc[i].z -= pf.z*invMass;
			acc[j].x += pf.x*invMass;
			acc[j].y += pf.y*invMass;
			acc[j].z += pf.z*invMass;
		}
	}*/
	
	#pragma omp parallel
	{
		#pragma omp for shared(acc, Up) schedule(dynamic, OMP_CHUNK) nowait
		for (unsigned int cellID = 0; cellID < cl_.nCells.x*cl_.nCells.y*cl_.nCells.z; ++cellID) {
			int atom1 = cl_.head(cellID);
			while (atom1 >= 0) {
				std::vector < int > neighbors = cl_.neighbors(cellID);
				for (unsigned int index = 0; index < neighbors.size(); ++index) {
					const int cellID2 = neighbors[index];
					int atom2 = cl_.head(cellID2);
					// only compute for i > j (prevents i == i and uses 3rd Law)
					while (atom2 >= 0) {
						if (atom1 > atom2) {
							float3 pf;
							Up += sys.potential (&sys.atoms[atom1].pos, &sys.atoms[atom2].pos, &pf, &box, &rc);
							acc[atom1].x -= pf.x*invMass;
							acc[atom1].y -= pf.y*invMass;
							acc[atom1].z -= pf.z*invMass;
							acc[atom2].x += pf.x*invMass;
							acc[atom2].y += pf.y*invMass;
							acc[atom2].z += pf.z*invMass;
						}
						atom2 = cl_.list(atom2);
					}
				}
				atom1 = cl_.list(atom1);
			}
		}
	}

	// apply thermal friction
	#pragma omp parallel
	{
		#pragma omp for shared(acc, sys.atoms) schedule(dynamic, OMP_CHUNK) nowait
		for (unsigned int atom1 = 0; atom1 < sys.numAtoms(); ++atom1) {
			acc[atom1].x -= gamma_*sys.atoms[atom1].vel.x;
			acc[atom1].y -= gamma_*sys.atoms[atom1].vel.y;
			acc[atom1].z -= gamma_*sys.atoms[atom1].vel.z;
			sys.atoms[atom1].acc = acc[atom1];
		}
	}

    	// update gamma for the NEXT time (finite differences)
    	const float gamma_dot = (sys.KinE()*2.0 - (3.0*sys.numAtoms()+1.0)*sys.targetT())/Q_;
    	gamma_ += gamma_dot*dt_;

    	// set Up
    	sys.setPotE(Up);
}