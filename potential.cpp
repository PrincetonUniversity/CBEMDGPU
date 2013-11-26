/*!
 * Pair Potentials
 * \author Nathan A. Mahynski
 * \date 11/18/13
 */

#include "potential.h"
#include "utils.h"
#include <iostream>
#include <math.h>

#ifdef NVCC
__device__ float pairUF (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *rcut) {
	return 0.0;
}
#else

/*!
 * Pairwise interaction between 2 atoms
 *
 * \param [in] p1 Pointer to atom 1's position
 * \param [in] p2 Pointer to atom 2's position
 * \param [in, out] pairForce Force atom 2 experiences due to atom 1
 * \param [in] box Pointer to box dimensions
 * \param [in] rcut Cutoff distance
 */
float pairUF (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *rcut) {
	// for now use the same potential as in the HW
	const float eps = 1.0, rc = *rcut;
	float3 dr;
	float r2 = pbcDist2 (*p1, *p2, dr, *box);
	if (r2 < rc*rc) {
		const float r = sqrt(r2);
		const float factor = -eps*rc*(1.0-r/rc);
		pairForce->x = -factor*dr.x/r;
		pairForce->y = -factor*dr.y/r;
		pairForce->z = -factor*dr.z/r;
		return 0.5*eps*rc*(1-r/rc)*(1-r/rc);
	} else {
		pairForce->x = 0.0;
		pairForce->y = 0.0;
		pairForce->z = 0.0;
		return 0.0;
	}
}
#endif

