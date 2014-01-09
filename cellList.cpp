/*!
 * Cell Lists
 * \author Nathan A. Mahynski
 * \date 11/19/13
 */ 

#include "cellList.h"
#include "dataTypes.h"	// even in GPU implementation, not compiling with NVCC so these data types are required
#include "common.h"
#include <math.h>
#include "system.h"
#include "utils.h"
#include <stdlib.h>

#ifdef NVCC
/*!
 * Initialize a neighbor list. If using cuda, "cell lists" are actually neighbor lists instead but are still maintained on the cpu
 *
 * \param [in] box Box size
 * \param [in] rc Cutoff radius
 * \param [in] rs Skin Radius
 */ 
cellList_cpu::cellList_cpu (const float3 &box, const float rc, const float rs) {
	if (rc < 0.0) {
        	throw customException("Cutoff radius must be > 0");
        	return;
   	}
    	rc_ = rc;
   	if (rs < 0.0) {
        	throw customException("Skin radius must be > 0");
        	return;
    	}
    	rs_ = rs;
	start_ = 1;
    	box_ = box;
}

/*!
 * Checks to see if the cell list needs to be updated and if so, rebuilds
 * Uses linked list, so algorithm in O(N)
 *
 * \param [in] sys System definition
 */
void cellList_cpu::checkUpdate (const systemDefinition &sys) {
	int build = 0;
	if (start_) {
		// resize the neighborlists, and initially set to empty
		try {
			nlist_index.resize(sys.numAtoms());
		} catch (std::exception& e) {
			std::cerr << e.what() << std::endl;
			throw customException ("Unable to allocate memory for neighbor list");
			return;
		}
		try {
                        posAtLastBuild_.resize(sys.numAtoms());
                } catch (std::exception& e) {
                        std::cerr << e.what() << std::endl;
                        throw customException ("Unable to initialize position list for neighbor list");
                        return;
                }

		for (unsigned int i = 0; i < nlist_index.size(); ++i) {
			nlist_index[i] = -1;
		}

		// must build when initialized
		build = 1;
	} else {
		float drMax1_ = 0.0;
                float drMax2_ = 0.0;
                float3 dummy, box = sys.box();
                for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
			const float dr2 = pbcDist2 (sys.atoms[i].pos, posAtLastBuild_[i], dummy, box);
			if (dr2 > drMax1_*drMax1_) {
				drMax2_ = drMax1_;
				drMax1_ = sqrt(dr2);
			} else if (dr2 > drMax2_*drMax2_) {
				drMax2_ = sqrt(dr2);
			}
		}
                if (drMax1_+drMax2_ > rs_) {
                        build = 1;
                } else {
                        build = 0;
                }
	}

	// check to rebuild the neighbor list
	if (build) {
		std::vector < std::vector < int > > dummyNeighbors (sys.numAtoms());
		std::vector < int > empty (1, -1);
                for (unsigned int i = 0; i < dummyNeighbors.size(); ++i) {
                        dummyNeighbors[i] = empty;
                }
		const int N = sys.numAtoms();
		std::vector < int > nn(N, 0);
		const float cut2 = (rs_+rc_)*(rs_+rc_);
		float3 dummy, box = sys.box();
		int totalNeighbors = 0;
		for (unsigned int i = 0; i < N; ++i) {
			for (unsigned int j = i+i; j < N; ++j) {
				float dist2 = pbcDist2(sys.atoms[i].pos, sys.atoms[j].pos, dummy, box);
				if (dist2 < cut2) {
					nn[i]++;
					nn[j]++;
					totalNeighbors += 2;

					// double the size to reduce the number of times this has to happen
					if (dummyNeighbors[i].size() < nn[i]) {
						dummyNeighbors[i].resize(2*dummyNeighbors[i].size());
					}
					if (dummyNeighbors[j].size() < nn[j]) {
                                                dummyNeighbors[j].resize(2*dummyNeighbors[j].size());
                                        }
					dummyNeighbors[i][nn[i]-1] = j;
					dummyNeighbors[j][nn[j]-1] = i;
				}
			}
		}
		
		// make the list "linear"
		try {
			nlist.resize(totalNeighbors + N);
		} catch (std::exception &e) {
			std::cerr << e.what() << std::endl;
			return;
		}
		
		int counter = 0;
		for (unsigned int i = 0; i < N; ++i) {
                        //dummyNeighbors[i].resize(nn[i]); // free as much memory as possible
                        posAtLastBuild_[i] = sys.atoms[i].pos;
			nlist[counter] = nn[i];
			nlist_index[i] = counter;
			counter++;
			for (unsigned int j = 0; j < nn[i]; ++j) {
				nlist[counter] = dummyNeighbors[i][j];	
				counter++;	
			}
			dummyNeighbors[i].resize(0); // free as much memory as possible
                }
	}
}	

// if not using GPUs (CUDA) maintain cell lists
#else
/*!
 * Initialize a cell list
 *
 * \param [in] box Box size
 * \param [in] rc Cutoff radius
 * \param [in] rs Skin Radius
 */
cellList_cpu::cellList_cpu (const float3 &box, const float rc, const float rs) {
    if (rc < 0.0) {
	throw customException("Cutoff radius must be > 0");
	return;
    }
    rc_ = rc;
    if (rs < 0.0) {
	throw customException("Skin radius must be > 0");
	return;
    }
    rs_ = rs;

    box_ = box;

	start_ = 1;
    lcell_.x = 1.01*(rc+rs);
    nCells.x = (int) floor (box.x/lcell_.x);
    lcell_.x = (box.x/nCells.x);

    lcell_.y = 1.01*(rc+rs);
    nCells.y = (int) floor (box.y/lcell_.y);
    lcell_.y = (box.y/nCells.y);

    lcell_.z = 1.01*(rc+rs);
    nCells.z = (int) floor (box.z/lcell_.z);
    lcell_.z = (box.z/nCells.z);

    if (lcell_.x <= (rc+rs) || lcell_.y <= (rc+rs) || lcell_.z < (rc+rs)) {
	throw customException("Cell width must exceed sum of cutoff and skin radius");
	return;
    }

    if (lcell_.x < 1) {
	throw customException ("Box dimension x too small relative to skin and cutoff radius to use cell lists");
	return;
    }
    if (lcell_.y < 1) {
	throw customException ("Box dimension y too small relative to skin and cutoff radius to use cell lists");
	return;
    }
    if (lcell_.z < 1) {
	throw customException ("Box dimension z too small relative to skin and cutoff radius to use cell lists");
	return;
    }

    if (nCells.x < 3 || nCells.y < 3 || nCells.z < 3) {
	throw customException("Must be able to have at least 3 cells in each direction, change box size, skin, or cutoff radius");
	return;
    }

    try {
	head_.resize(nCells.x*nCells.y*nCells.z, -1);
    } catch (std::exception &e) {
	std::cerr << e.what() << std::endl;
	throw customException ("Unable to initialize head for cell list");
	return;
    }
    // build neighbors for each cell
	neighbor_.resize(nCells.x*nCells.y*nCells.z);
	for (unsigned int cellID = 0; cellID < nCells.x*nCells.y*nCells.z; ++cellID) {
	    const int zref = floor(cellID/(nCells.x*nCells.y));
	    const int yref = floor((cellID - zref*(nCells.x*nCells.y))/nCells.y);
	    const int xref = floor(cellID - zref*(nCells.x*nCells.y) - yref*nCells.x);
	    for (int dx = -1; dx <= 1; ++dx) {
		int xcell = xref + dx;
		if (xcell >= nCells.x) xcell = 0;
		if (xcell < 0) xcell = nCells.x-1;
		for (int dy = -1; dy <= 1; ++dy) {
		    int ycell = yref + dy;
		    if (ycell >= nCells.y) ycell = 0;
		    if (ycell < 0) ycell = nCells.y-1;
		    for (int dz = -1; dz <= 1; ++dz) {
			int zcell = zref + dz;
			if (zcell >= nCells.z) zcell = 0;
			if (zcell < 0) zcell = nCells.z-1;
			const int cellID2 = xcell + ycell*nCells.x + zcell*(nCells.x*nCells.y);
			neighbor_[cellID].push_back(cellID2);
		    }
		}
	    }
	}
    for (unsigned int cellID = 0; cellID < nCells.x*nCells.y*nCells.z; ++cellID) {
	if (neighbor_[cellID].size() != 27) {
	    throw customException ("Cell list initial build failed to find 27 neighbors (including self)");
	    return;
	}
    }
}

/*!
 * Calculates the cell (linear index of it) a position belongs to in a periodic box.
 * The position does not need to be "inside" the box, boundary conditions are applied.
 *
 * \param [in] pos Position of atom
 * \return cell
 */ 
int cellList_cpu::cell (const float3 &pos) {
    float3 inBox = pbc(pos, box_);
    const int x = (int) floor(inBox.x/lcell_.x);
    const int y = (int) floor(inBox.y/lcell_.y);
    const int z = (int) floor(inBox.z/lcell_.z);
    return x + y*nCells.x + z*nCells.x*nCells.y;
}

/*!
 * Checks to see if the cell list needs to be updated and if so, rebuilds
 * Uses linked list, so algorithm in O(N)
 *
 * \param [in] sys System definition
 */
void cellList_cpu::checkUpdate (const systemDefinition &sys) {
    int build = 0;
    if (start_) {
	const int natoms = sys.numAtoms();
	try {
			list_.resize(natoms, -1);
		} catch (std::exception& e) {
			std::cerr << e.what() << std::endl;
			throw customException ("Unable to initialize list for cell list");
			return;
		}
		try {
			posAtLastBuild_.resize(natoms);
		} catch (std::exception& e) {
			std::cerr << e.what() << std::endl;
			throw customException ("Unable to initialize position list for cell list");
			return;
		}
		start_ = 0;
		build = 1;
	} else {
		// check max displacements
		float drMax1_ = 0.0;
		float drMax2_ = 0.0;
		float3 dummy;
			for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
				const float dr2 = pbcDist2 (sys.atoms[i].pos, posAtLastBuild_[i], dummy, sys.box());
				if (dr2 > drMax1_*drMax1_) {
					drMax2_ = drMax1_;
					drMax1_ = sqrt(dr2);
				} else if (dr2 > drMax2_*drMax2_) {
					drMax2_ = sqrt(dr2);
				}
			}
		if (drMax1_+drMax2_ > rs_) {
			build = 1;
		} else {
			build = 0;
		}
	}

	if (build) {
		if (sys.numAtoms() != list_.size()) {
			throw customException ("Number of atoms in simulation has changed");
			return;
		}
		for (unsigned int i = 0; i < head_.size(); ++i) {
			head_[i] = -1;
		}
			for (unsigned int i = 0; i < sys.numAtoms(); ++i) {
				const int icell = cell(sys.atoms[i].pos);
				list_[i] = head_[icell];
				head_[icell] = i;
				posAtLastBuild_[i] = sys.atoms[i].pos;
			}
	} 

	return;
}
#endif
