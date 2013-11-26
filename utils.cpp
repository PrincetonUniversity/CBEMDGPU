/*!
 * Utility functions
 * \author Nathan A. Mahynski
 * \date 11/18/13
 */

#include "dataTypes.h"
#include <math.h>

#ifdef NVCC
__device__ float dev_pbcDist2 (const float3 *p1, const float3 *p2, float3 *dr, const float3 *box) {
	float d = 0.0;
	
	dr->x = p2->x - p1->x;
	while (dr->x > box->x/2.0) {
		dr->x -= box->x;
	}
	while (dr->x <= box->x/2.0) {
		dr->x += box->x;
	}
	d += dr->x*dr->x;
	
	dr->y = p2->y - p1->y;
	while (dr->y > box->y/2.0) {
		dr->y -= box->y;
	} 
	while (dr->y <= box->y/2.0) {
		dr->y += box->y;
	}
	d += dr->y*dr->y;
	
	dr->z = p2.z - p1.z;
	while (dr->z > box.z/2.0) {
		dr->z -= box.z;
	}
	while (dr->z <= -box.z/2.0) {
		dr->z += box.z;
	}
	d += dr->z*dr->z;
	
	return d;
}
#endif

float3 pbc (const float3 &p1, const float3 &box) {
	float3 ans = p1;
	while (ans.x < 0.0) {
		ans.x += box.x;
	}
	while (ans.x >= box.x) {
		ans.x -= box.x;
	}
	while (ans.y < 0.0) {
		ans.y += box.y;
	}
	while (ans.y >= box.y) {
		ans.y -= box.y;
	}
	while (ans.z < 0.0) {
		ans.z += box.z;
	}
	while (ans.z >= box.z) {
		ans.z -= box.z;
	}
	return ans;
}
	
float pbcDist2 (const float3 &p1, const float3 &p2, float3 &dr, const float3 &box) {
	double d = 0.0;

	dr.x = p2.x - p1.x;
	while (dr.x > box.x/2.0) {
		dr.x -= box.x;
	}
	while (dr.x <= -box.x/2.0) {
		dr.x += box.x;
	}
	d += dr.x*dr.x;

	dr.y = p2.y - p1.y;
	while (dr.y > box.y/2.0) {
		dr.y -= box.y;
    	} 
	while (dr.y <= -box.y/2.0) {
		dr.y += box.y;
	}
	d += dr.y*dr.y;

	dr.z = p2.z - p1.z;
	while (dr.z > box.z/2.0) {
		dr.z -= box.z;
	}
	while (dr.z <= -box.z/2.0) {
		dr.z += box.z;
	}
	d += dr.z*dr.z;

	return d;
}
