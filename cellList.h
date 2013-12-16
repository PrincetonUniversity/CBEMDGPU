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
 * Maintains a neighbor list for each atom on the CPU
 */ 
class cellList_cpu {
        public:
                cellList_cpu () {}
                cellList_cpu (const float3 &box, const float rc, const float rs);
                ~cellList_cpu () {}
                void checkUpdate (const systemDefinition &sys);
		std::vector < int > neighborList (const int atomIndex) {return neighbors_[atomIndex];}
        private:
		int start_;
                float rc_, rs_;
                float3 box_;
                float drMax1_, drMax2_;
                std::vector <float3> posAtLastBuild_;
		std::vector < std::vector < int > > neighbors_;
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
		void checkUpdate (const systemDefinition &sys);
		int cell (const float3 &pos);
		int head (const int cell) const {return head_[cell];}
		int list (const int index) const {return list_[index];}
		int3 nCells; // this is commonly used and is easier to have as public
		std::vector < int > neighbors (const int cellID) const {return neighbor_[cellID];}
	private:
		int start_;
		std::vector < std::vector < int > > neighbor_;
		float rc_, rs_;
	        float3 lcell_, box_;
		float drMax1_, drMax2_;
		std::vector <float3> posAtLastBuild_;
		std::vector <int> head_;
		std::vector <int> list_;
};

#endif

#endif
