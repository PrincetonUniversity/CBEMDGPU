

#ifndef __DATA_TYPES_H__
#define __DATA_TYPES_H__

#ifdef NVCC
#include "/usr/local/cuda/include/builtin_types.h"
#else
struct float3 {
		float x, y, z;
};
typedef float3 float3;
struct int3 {
	int x, y, z;
};
typedef int3 int3;
#endif

struct atom {
		float3 pos, vel, acc;
};
typedef atom atom;

#endif
