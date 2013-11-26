/*!
 * Utility functions
 * \author Nathan A. Mahynski
 * \date 11/19/13
 */ 

#ifndef __UTILS_H__
#define __UTILS_H__

#include "dataTypes.h"

#ifdef NVCC
__device__ float dev_pbcDist2 (const float3 *p1, const float3 *p2, float3 *dr, const float3 *box);
#endif
float pbcDist2 (const float3 &p1, const float3 &p2, float3 &dr, const float3 &box);
float3 pbc (const float3 &p1, const float3 &box);

#endif
