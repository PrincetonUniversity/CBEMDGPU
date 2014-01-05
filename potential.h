/*!
 * Pair Potentials
 * \date 11/18/13
 */

#ifndef __POTENTIAL_H__
#define __POTENTIAL_H__

#ifdef NVCC
// These functions actually exist in integrator.cu because of how they must be compiled
__device__ float dev_pairUF (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *args, const float *rcut);
__device__ float dev_slj (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *args, const float *rcut);
#else
#include "dataTypes.h"
float pairUF (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *args, const float *rcut);
float slj (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *args, const float *rcut);
#endif

//!< Function pointer that all force calculation (pair potentials) must follow
typedef float(*pointFunction_t)(const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *args, const float *rcut);

#endif
