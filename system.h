/*!
 * System definition
 * \date 11/17/13
 */

#ifndef __SYSTEM_DEFINITION_H__
#define __SYSTEM_DEFINITION_H__

#include <iostream>
#include <stdio.h>
#include <vector>
#include "dataTypes.h"
#include "potential.h"

//! Contains all information pertaining to a system being simulated.
class systemDefinition {
	public:
		systemDefinition () {mass_ = -1; instantT_ = 0; targetT_ = 0; snapFile_ = NULL; Uk_ = 0.0; Up_ = 0.0; rc_ = 0; rs_ = 0;}
		~systemDefinition () {if (snapFile_ != NULL) fclose(snapFile_);}
		void initRandom (const int N, const int rngSeed);
		void initThermal (const int N, const float Tset, const int rngSeed, const float dx);
		void updateInstantTemp (const float T) {instantT_ = T;} //!< Manually set the instantaneous temperature
		void setTemp (const float T) {targetT_ = T;}    //!< Assign the target temperature for NVT simulations
		void setMass (const float m) {mass_ = m;}       //!< Assign the mass of each particle
		void setRcut (const float rc) {rc_ = rc;}       //!< Assign the cutoff radius of the pair potential
		void setRskin (const float rs) {rs_ = rs;}      //!< Assign the skin radius as a buffer for the neighbor/cell lists
		void setBox(const float x, const float y, const float z) {box_.x = x; box_.y = y; box_.z = z;}  //!< Assign the simulation box size
		void printBox() {std::cout << box_.x << "\t" << box_.y << "\t" << box_.z << std::endl;} //!< Print the box dimesions to stdout
		float3 box() const {return box_;}   //!< Report the box dimensions
		float instantT() const {return instantT_;}  //!< Report the instantaneous temperature
		float targetT() const {return targetT_;}    //!< Report the target temperature for NVT simulations
		float mass() const {return mass_;}          //!< Report the particle's mass
		float PotE() const {return Up_;}            //!< Report the instantaneous potential energy of the system
		void setPotE(const float Up) {Up_ = Up;}    //!< Assign the potential energy
		void setKinE(const float Uk) {Uk_ = Uk;}    //!< Assign the kinetic energy
		float KinE() const {return Uk_;}            //!< Report the instantaneous kinetic energy of the system
		float rskin() const {return rs_;}           //!< Report the skin radius for the neighbor/cell lists
		float rcut() const {return rc_;}            //!< Report the pair potential function cutoff radius
		void writeSnapshot ();                      //!< Print a snapshot of the system
		int numAtoms() const {return atoms.size();} //!< Report the pair potential function cutoff radius
		void setPotentialArgs (const std::vector <float> args) {potentialArgs_ = args;} //!< Assign additional arguments to the pair potential function
		std::vector <float> potentialArgs () {return potentialArgs_;}   //!< Report additional arguments to the pair potential function
		
		#ifdef NVCC
		int cudaBlocks, cudaThreads;            //!< Block and thread size if using GPUs
		#endif

		pointFunction_t potential;              //!< Pointer to pair potential function
		void setPotential (pointFunction_t pp); //!< Assign pair potential function
		std::vector <atom> atoms;               //!< Vector of atoms in the system
    
	private:
        float rc_;              //!< Cutoff radius of the pair potential function
        float rs_;              //!< Skin radius for neighbor/cell lists
		FILE *snapFile_;        //!< File to record system's trajectory to
		float3 box_;            //!< Box dimensions
        float targetT_;         //!< Target temperature for NVT simulations
        float instantT_;        //!< Instantaneous (kinetic) temperature of the system
		float mass_;            //!< Particle mass
        float Uk_;              //!< Potential energy
        float Up_;              //!< Kinetic energy
		std::vector <float> potentialArgs_; //!< Additional arguments to the pair potential function
};

#endif
