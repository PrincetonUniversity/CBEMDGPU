/*!
 * Integration
 * \author Nathan A. Mahynski
 * \date 11/19/13
 */

#include "system.h"
#include <math.h>
#include "common.h"
#include "integrator.h"
#include <vector>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/system_error.h>
#include "potential.h"

// Cuda Potential must be compiled together with the integrator so they are here, not in their own file
__device__ float dev_pbcDist2 (const float3 *p1, const float3 *p2, float3 *dr, const float3 *box) {
	double d = 0.0;

	dr->x = p2->x - p1->x;
	while (dr->x > box->x/2.0) {
		dr->x -= box->x;
	}
	while (dr->x <= -box->x/2.0) {
		dr->x += box->x;
	}
	d += dr->x*dr->x;

	dr->y = p2->y - p1->y;
	while (dr->y > box->y/2.0) {
		dr->y -= box->y;
    	} 
	while (dr->y <= -box->y/2.0) {
		dr->y += box->y;
	}
	d += dr->y*dr->y;

	dr->z = p2->z - p1->z;
	while (dr->z > box->z/2.0) {
		dr->z -= box->z;
	}
	while (dr->z <= -box->z/2.0) {
		dr->z += box->z;
	}
	d += dr->z*dr->z;

	return d;
}

/*!
 * An arbitrary potential that could be changed to meet the user's needs in the future.
 */
__device__ float dev_pairUF (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *args, const float *rcut) {
        return 0.0;
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
__device__ float dev_slj (const float3 *p1, const float3 *p2, float3 *pairForce, const float3 *box, const float *args, const float *rcut) {
        float3 dr;
	float r2 = dev_pbcDist2(p1, p2, &dr, box);

	// check that r > delta and throw/catch
	if (r2 <= args[2]*args[2]) {
		// shift to just past the singularity and allow to run
		r2 = 1.0001*args[2]*args[2];
	}
	
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

/*!
 * From the host, call this kernel to loop over each atom's neighbor list
 *
 * \param [in] dev_atoms Atoms (sys.atoms)
 * \param [in] nlist Neighbor list of all atoms indexed in a linear array
 * \param [in] nlist_index Index to start at in nlist to find an atom's neighbors
 * \param [out] force Net force each atom experiences
 * \param [in] Box dimensions
 * \param [out] Up_each Potential energy of each atom in the field, the total is this summed divided by 2
 * \param [in] natoms Number of atoms
 * \param [in] args Pair potential arguments
 * \param [in] rcut Cutoff distance for potential
 * \param [in] pFlag Flags which potential function to use
 */
__global__ void loopOverNeighbors (atom* dev_atoms, int* nlist, int* nlist_index, float3* force, float3* box, float* Up_each, int* natoms, float* args, float* rcut, int* pFlag) {
	const int tid = threadIdx.x + blockIdx.x*blockDim.x;
	
	if (tid < *natoms) {
		float Up = 0.0;
		const int nn = nlist_index[tid]; // number of neighbors
		
		// choose the potential function
		pointFunction_t pFunc;
		if (*pFlag == 0) {
			pFunc = dev_slj;
		} else {
			pFunc = dev_pairUF;
		}

		float3 myforce;
		myforce.x = 0;
		myforce.y = 0;
		myforce.z = 0;

		// loop over this atom's neighbors
		for (unsigned int i = nlist_index[tid]+1; i < nlist_index[tid]+1+nn; ++i) {
			// compute potential between dev_atoms[nlist[i]] and dev_atoms[tid]
			float3 dummyForce;
			Up += pFunc (&dev_atoms[nlist[i]].pos, &dev_atoms[tid].pos, &dummyForce, box, args, rcut);

			// sign convention is to add the returned force on particle 2
			myforce.x += dummyForce.x;
			myforce.y += dummyForce.y;
			myforce.z += dummyForce.z;	
		}

		// carmeline suggested a sign change
		force[tid].x = -myforce.x;
		force[tid].y = -myforce.y;
		force[tid].z = -myforce.z;
		Up_each[tid] = Up;
	}	
}

/*!
 * Calculate the pairwise forces in a system.  This also calculates the potential energy of a system.
 * The kinetic energy is calculated during the verlet integration.
 *
 * \param [in, out] sys System definition
 */
void integrator::calcForce (systemDefinition &sys) {
	float Up = 0.0;
	const float invMass = 1.0/sys.mass();

	// set which potential to use
	std::vector < int > pFlag (1, -1);
	if (sys.potential == dev_slj) {
		pFlag[0] = 0;
	} else if (sys.potential == dev_pairUF) {
		pFlag[0] = 1;
	} else {
		throw customException ("Cannot understand which potential function to use");
	}
	thrust::device_vector < int > dev_pFlag (pFlag.begin(), pFlag.end());
	int* dev_pFlag_ptr = thrust::raw_pointer_cast(&dev_pFlag[0]);

	// box dimensions
	std::vector < float3 > box (1, sys.box());
	thrust::device_vector < float3 > sysbox (box.begin(), box.end());
	float3* dev_sysbox_ptr = thrust::raw_pointer_cast(&sysbox[0]);

	// check update for neighborlist
	cl_.checkUpdate(sys);
	
	std::cout << "cp00" << std::endl;

	// number of atoms
	thrust::device_vector < int > dev_natoms(1, sys.numAtoms());
	int* dev_natoms_ptr = thrust::raw_pointer_cast(&dev_natoms[0]);

	std::cout << "cp01" << std::endl;

	// copy sys.atoms to GPU
	thrust::device_vector < atom > dev_atoms (sys.atoms.begin(), sys.atoms.end());
	atom* dev_atoms_ptr = thrust::raw_pointer_cast(&dev_atoms[0]);

	std::cout << "cp02" << std::endl;
	
	// copy neighborlists to GPU
	thrust::device_vector < int > dev_neighbor_list (cl_.nlist.begin(), cl_.nlist.end());
	int* dev_neighbor_list_ptr = thrust::raw_pointer_cast(&dev_neighbor_list[0]);

	std::cout << "cp03" << std::endl;

	// copy location of where the neighbors for each atoms starts
	thrust::device_vector < int > dev_neighbor_index (cl_.nlist_index.begin(), cl_.nlist_index.end());
	int* dev_neighbor_index_ptr = thrust::raw_pointer_cast(&dev_neighbor_index[0]);
	
	std::cout << "cp04" << std::endl;
	
	// create acc and Up to store results in
	thrust::device_vector < float3 > dev_force (sys.numAtoms());
	float3* dev_force_ptr = thrust::raw_pointer_cast(&dev_force[0]);
	thrust::device_vector < float > dev_Up_each_atom (sys.numAtoms());
	float* dev_Up_each_atom_ptr = thrust::raw_pointer_cast(&dev_Up_each_atom[0]);
	
	std::cout << "cp05" << std::endl;

	// potential arguments
	std::vector < float > potArgs = sys.potentialArgs();
	std::vector < float > rcut (1, sys.rcut());
	thrust::device_vector < float > dev_args (potArgs.begin(), potArgs.end());
	float* dev_args_ptr = thrust::raw_pointer_cast(&dev_args[0]);
	thrust::device_vector < float > dev_rcut (rcut.begin(), rcut.end());
	float* dev_rcut_ptr = thrust::raw_pointer_cast(&dev_rcut[0]);

	std::cout << "cp06" << std::endl;

	// invoke kernel to compute
	loopOverNeighbors <<< sys.cudaBlocks, sys.cudaThreads >>> (dev_atoms_ptr, dev_neighbor_list_ptr, dev_neighbor_index_ptr, dev_force_ptr, dev_sysbox_ptr, dev_Up_each_atom_ptr, dev_natoms_ptr, dev_args_ptr, dev_rcut_ptr, dev_pFlag_ptr);

	std::cout << "cp07" << std::endl;

	// call a reduction to collect Up then divide by 2 since double counted
	/*Up = thrust::reduce(dev_Up_each_atom.begin(), dev_Up_each_atom.end(), (float) 0.0, thrust::plus<float>());
	Up /= 2.0;	// pairs are double counted
	*/

	//tmp
	std::vector < float > tmpUP (sys.numAtoms());
	thrust::copy(dev_Up_each_atom.begin(), dev_Up_each_atom.end(), tmpUP.begin());
	float sum = 0.0;
	for (int i = 0; i < sys.numAtoms(); ++i) {
		sum += tmpUP[i];
	}
	Up = sum/2.0;
	
	std::cout << "cp08" << std::endl;

	// store accelerations on atoms	
	std::vector < float3 > netForces (sys.numAtoms());
	thrust::copy(dev_force.begin(), dev_force.end(), netForces.begin());

	std::cout << "cp09" << std::endl;

	for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
		sys.atoms[i].acc.x = netForces[i].x*invMass;
		sys.atoms[i].acc.y = netForces[i].y*invMass;
		sys.atoms[i].acc.z = netForces[i].z*invMass;
	}	

	std::cout << "cp10" << std::endl;

	// set Up
	sys.setPotE(Up);
}
