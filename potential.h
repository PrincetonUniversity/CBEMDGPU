/*!
 * Pair Potentials
 * \author Nathan A. Mahynski
 * \date 11/18/13
 */

#ifndef __POTENTIAL_H__
#define __POTENTIAL_H__

#include "dataTypes.h"

typedef float(*pointFunction_t)(const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *rcut);

#ifdef NVCC
__device__ float pairUF (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *rcut);
#else
float pairUF (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *rcut);
#endif

#endif
