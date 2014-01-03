/*!
 * Integration
 * \author Nathan A. Mahynski
 * \date 11/19/13
 */

#include "system.h"
#include "dataTypes.h"
#include <exception>
#include <math.h>
#include "common.h"
#include "integrator.h"
#include <omp.h>
#include <vector>

#ifndef NVCC
/*!
 * Calculate the pairwise forces in a system.  This also calculates the potential energy of a system.
 * The kinetic energy is calculated during the verlet integration.
 *
 * \param [in, out] sys System definition
 */
void integrator::calcForce (systemDefinition &sys) {
	// For cache coeherency allocate new space for calculations
	float3 empty;
	empty.x = 0; empty.y = 0; empty.z = 0;
	std::vector <float3> acc (sys.numAtoms(), empty);
	const float rc = sys.rcut();

	// every time, check if the cell list needs to be updated first
	cl_.checkUpdate(sys);

	float Up = 0.0;
	const float3 box = sys.box();
	const float invMass = 1.0/sys.mass();
	
	// traverse cell list and calculate total system potential energy 
	std::vector <float> args = sys.potentialArgs();
	#pragma omp parallel for reduction(+:Up)
	for (unsigned int cellID = 0; cellID < cl_.nCells.x*cl_.nCells.y*cl_.nCells.z; ++cellID) {
		int atom1 = cl_.head(cellID);
		while (atom1 >= 0) {
			std::vector < int > neighbors = cl_.neighbors(cellID);
			for (int index = 0; index < neighbors.size(); ++index) {
				const int cellID2 = neighbors[index];
				int atom2 = cl_.head(cellID2);
				while (atom2 >= 0) {
					if (atom1 > atom2) {
						float3 pf;
						Up += sys.potential (&sys.atoms[atom1].pos, &sys.atoms[atom2].pos, &pf, &box, &args[0], &rc);
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
	
	// save acceleration in array of atoms in system
	#pragma omp parallel for
	for (int i = 0; i < sys.numAtoms(); ++i) {
		sys.atoms[i].acc.x = -acc[i].x;
		sys.atoms[i].acc.y = -acc[i].y;
		sys.atoms[i].acc.z = -acc[i].z;
	}
	
	// set Up
	sys.setPotE(Up);
}

#endif
