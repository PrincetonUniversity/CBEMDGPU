/*!
 * Cell Lists
 * \author Nathan A. Mahynski
 * \date 11/19/13
 */

#ifndef __CELL_LIST_H__
#define __CELL_LIST_H__

#include <vector>
#include "dataTypes.h"
#include "system.h"

#ifdef NVCC
/*
 * Maintains a neighbor list for each atom on the CPU if using GPU implementation.
 */ 
class cellList_cpu {
    public:
        cellList_cpu () {}
        cellList_cpu (const float3 &box, const float rc, const float rs);
        ~cellList_cpu () {}
        void checkUpdate (const systemDefinition &sys); //!< Check if the neighbor list requires updating
        std::vector < int > nlist_index;    //!< Position in the neighbor list indicating where each particle's neighbors start from
        std::vector < int > nlist;          //!< Neighbor list containing the indices of each particles neighbors
	private:
        int start_;         //!< Flag indicating whether this list has been build before or not
        float rc_;          //!< Cutoff radius for pair potential
        float rs_;          //!< Skin radius for neighbor lists
        float3 box_;        //!< Simulation box size (x,y,z)
        float drMax1_;      //!< Largest displacement of a particle since the last build
        float drMax2_;      //!< Second largest displacement of a particle since the last build
        std::vector <float3> posAtLastBuild_;   //!< Coordinates of each particle since the last time the list was built
};
#else
/*!
 * Maintains a linked list to track cells on the CPU
 */ 
class cellList_cpu {
	public:
		cellList_cpu () {}
		cellList_cpu (const float3 &box, const float rc, const float rs);
		~cellList_cpu () {}
		void checkUpdate (const systemDefinition &sys); //!< Check if the neighbor list requires updating
		int cell (const float3 &pos);   //!< Calculate the cell in which a given coordinate is located
		int head (const int cell) const {return head_[cell];}   //!< Return the first atom (aka 'head') of each cell
		int list (const int index) const {return list_[index];} //!< Iterates through a linked list of particles, returns the index of the next atom in line
		int3 nCells; //!< Number of cells in each direction
		std::vector < int > neighbors (const int cellID) const {return neighbor_[cellID];}  //!< Returns the indices of a cell's neighboring cells
	private:
		int start_; //!< Flag indicating whether this list has been build before or not
		std::vector < std::vector < int > > neighbor_;  //!< Stores the indices of a cell's neighboring cells
        float rc_;      //!< Cutoff radius for pair potential
        float rs_;      //!< Skin radius for cell lists
        float3 lcell_;  //!< Length of a cell in each cartesian direction
        float3 box_;    //!< Simulation box size (x,y,z)
        float drMax1_;  //!< Largest displacement of a particle since the last build
        float drMax2_;  //!< Second largest displacement of a particle since the last build
		std::vector <float3> posAtLastBuild_;   //!< Coordinates of each particle since the last time the list was built
		std::vector <int> head_;    //!< Stores the first atom (aka 'head') of each cell
		std::vector <int> list_;    //!< Stores the linked list of particles in each cell
};

#endif

#endif
