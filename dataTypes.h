/*!
 * Requisite 'optimal' data types.
 * \date 11/21/13
 */

#ifndef __DATA_TYPES_H__
#define __DATA_TYPES_H__

#ifdef NVCC
// this should be included automatically, but just in case...
#include "/usr/local/cuda/include/builtin_types.h"
#else

//! 3 floating point numbers, same as defined for GPUs
struct float3 {
		float x, y, z;
};
typedef float3 float3;

//! 3 integers, same as defined for GPUs
struct int3 {
	int x, y, z;
};
typedef int3 int3;
#endif

//! Structure of an atom.
struct atom {
    float3 pos; //!< Position
    float3 vel; //!< Velocity
    float3 acc; //!< Acceleration
};
typedef atom atom;

#endif
