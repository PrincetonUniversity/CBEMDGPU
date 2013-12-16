/*!
 * Integration
 * \author Nathan A. Mahynski
 * \date 11/19/13
 */

#include "system.h"
#include <exception>
#include <math.h>
#include "common.h"
#include "integrator.h"
#include <vector>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

/*!
 * From the host, call this kernel to loop over each atom's neighbor list
 *
 * \param [in] dev_atoms Atoms (sys.atoms)
 * \param [in] nlist Neighbor list of all atoms indexed in a linear array
 * \param [in] nlist_index Index to start at in nlist to find an atom's neighbors
 * \param [out] force Net force each atom experiences
 * \param [out] Up Potential energy of each atom in the field, the total is this summed divided by 2
 * \param [in] natoms Number of atoms
 */
__global__ void loopOverNeighbors (atom* dev_atoms, int* nlist, int* nlist_index, float3* force, float* Up_each, int* natoms) {
	const int tid = threadIdx.x + blockIdx.x*blockDim.x;
	if (tid < *natoms) {
		const int nn = nlist_index[tid]; // number of neighbors
		for (unsigned int i = nlist_index[tid]+1; i < nlist_index[tid]+1+nn; ++i) {
			// compute potential between dev_atoms[nlist[i]] and dev_atoms[tid]
			
		}
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
	const float3 box = sys.box();
	const float invMass = 1.0/sys.mass();
	
	// check update for neighborlist
	cl_.checkUpdate(sys);

	// number of atoms
	thrust::device_vector < int > dev_natoms(1, sys.numAtoms());
	int* dev_natoms_ptr = thrust::raw_pointer_cast(&dev_natoms[0]);

	// copy sys.atoms to GPU
	thrust::device_vector < atom > dev_atoms (sys.atoms.begin(), sys.atoms.end());
	atom* dev_atoms_ptr = thrust::raw_pointer_cast(&dev_atoms[0]);

	// copy neighborlists to GPU
	thrust::device_vector < int > dev_neighbor_list (cl_.nlist.begin(), cl_.nlist.end());
	int* dev_neighbor_list_ptr = thrust::raw_pointer_cast(&dev_neighbor_list[0]);

	// copy location of where the neighbors for each atoms starts
	thrust::device_vector < int > dev_neighbor_index (cl_.nlist_index.begin(), cl_.nlist_index.end());
	int* dev_neighbor_index_ptr = thrust::raw_pointer_cast(&dev_neighbor_index[0]);
	
	// create acc and Up to store results in
	thrust::device_vector < float3 > dev_force (sys.numAtoms());
	float3* dev_force_ptr = thrust::raw_pointer_cast(&dev_force[0]);
	thrust::device_vector < float > dev_Up_each_atom (sys.numAtoms());
	float* dev_Up_each_atom_ptr = thrust::raw_pointer_cast(&dev_Up_each_atom[0]);

	// invoke kernel to compute
	loopOverNeighbors <<< sys.cudaBlocks, sys.cudaThreads >>> (dev_atoms_ptr, dev_neighbor_list_ptr, dev_neighbor_index_ptr, dev_force_ptr, dev_Up_each_atom_ptr, dev_natoms_ptr);

	// call a reduction to collect Up then divide by 2 since double counted
	Up = thrust::reduce(dev_Up_each_atom.begin(), dev_Up_each_atom.end(), (float) 0.0, thrust::plus<float>());
	Up /= 2.0;	// pairs are double counted

	// store accelerations on atoms	
	std::vector < float3 > netForces (sys.numAtoms());
	thrust::copy(dev_force.begin(), dev_force.end(), netForces.begin());
	for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
		sys.atoms[i].acc.x = netForces[i].x*invMass;
		sys.atoms[i].acc.y = netForces[i].y*invMass;
		sys.atoms[i].acc.z = netForces[i].z*invMass;
	}	

	// set Up
	sys.setPotE(Up);
}
