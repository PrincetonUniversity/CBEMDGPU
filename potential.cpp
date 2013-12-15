/*!
 * Pair Potentials
 * \author Nathan A. Mahynski
 * \date 11/18/13
 */

#include "potential.h"
#include "utils.h"
#include <iostream>
#include "common.h"
#include <math.h>

#ifdef NVCC
__device__ float pairUF (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *args, const float *rcut) {
	return 0.0;
}
__device__ float slj (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *args, const float *rcut) {
	return 0.0
}
#else

/*!
 * Pairwise interaction between 2 atoms (DPD)
 *
 * \param [in] p1 Pointer to atom 1's position
 * \param [in] p2 Pointer to atom 2's position
 * \param [in, out] pairForce Force atom 2 experiences due to atom 1
 * \param [in] box Pointer to box dimensions
 * \param [in] args Additional arguments, in this case epsilon
 * \param [in] rcut Cutoff distance
 */
float pairUF (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *args, const float *rcut) {
	// for now use the same potential as in the HW
	const float eps = args[0], rc = *rcut;
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


/*!
 * Pairwise interaction between 2 atoms (shifted lennard-jones)
 *
 * \param [in] p1 Pointer to atom 1's position
 * \param [in] p2 Pointer to atom 2's position
 * \param [in, out] pairForce Force atom 2 experiences due to atom 1
 * \param [in] box Pointer to box dimensions
 * \param [in] args Additional arguments, in this case {epsilon, sigma, delta, ushift}
 * \param [in] rcut Cutoff distance; NOTE: this MUST already incorporate the delta shift, ie. if r < ((rc' + delta) = rc), then force is computed
 */
float slj (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *args, const float *rcut) {
	float3 dr;
	float r2 = pbcDist2(*p1, *p2, dr, *box);
	// check that r > delta and throw/catch
	if (r2 < args[2]*args[2]) {
		throw customException("dr < delta");
		pairForce->x = 0.0;
		pairForce->y = 0.0;
		pairForce->z = 0.0;
		return 0.0;
	}
	// If (r-delta)^2 < rcut^2 compute
	if (r2 < (*rcut)*(*rcut)) {
		float r = sqrt(r2);
        	float x = r - args[2];
		float b = 1.0/x, a = args[1]*b, a2 = a*a, a6 = a2*a2*a2, factor;
		factor = 24.0*args[0]*a6*(2.0*a6-1.0)*b/r;
		pairForce->x = -factor*dr.x;
		pairForce->y = -factor*dr.y;
		pairForce->z = -factor*dr.z;
		return 4.0*args[0]*(a6*a6-a6)+args[3];
	} else {
		pairForce->x = 0.0;
		pairForce->y = 0.0;
		pairForce->z = 0.0;
		return 0.0;
	}
}
	
	
#endif

