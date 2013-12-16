/*!
 * System definition
 * \author Nathan A. Mahynski
 * \date 11/17/13
 */

#ifndef __SYSTEM_DEFINITION_H__
#define __SYSTEM_DEFINITION_H__

#include <iostream>
#include <stdio.h>
#include <vector>
#include "dataTypes.h"
#include "potential.h"

class systemDefinition {
	public:
		systemDefinition () {mass_ = -1; instantT_ = 0; targetT_ = 0; snapFile_ = NULL; Uk_ = 0.0; Up_ = 0.0; rc_ = 0; rs_ = 0;}
		~systemDefinition () {if (snapFile_ != NULL) fclose(snapFile_);}
		void initRandom (const int N, const int rngSeed);
		void initThermal (const int N, const float Tset, const int rngSeed, const float dx);
		void updateInstantTemp (const float T) {instantT_ = T;}
		void setTemp (const float T) {targetT_ = T;}
		void setMass (const float m) {mass_ = m;}
		void setRcut (const float rc) {rc_ = rc;}
		void setRskin (const float rs) {rs_ = rs;}
		void setBox(const float x, const float y, const float z) {box_.x = x; box_.y = y; box_.z = z;}
		void printBox() {std::cout << box_.x << "\t" << box_.y << "\t" << box_.z << std::endl;}
		float3 box() const {return box_;}
		float instantT() const {return instantT_;}
		float targetT() const {return targetT_;}
		float mass() const {return mass_;}
		float PotE() const {return Up_;}
		void setPotE(const float Up) {Up_ = Up;}
		void setKinE(const float Uk) {Uk_ = Uk;}
		float KinE() const {return Uk_;}
		float rskin() const {return rs_;}
		float rcut() const {return rc_;}
		void writeSnapshot ();
		int numAtoms() const {return atoms.size();}
		void setPotentialArgs (const std::vector <float> args) {potentialArgs_ = args;}
		std::vector <float> potentialArgs () const {return potentialArgs_;}
//#ifdef NVCC
//		__device__ pointFunction_t dev_potential;
//#endif
		#ifdef NVCC
		int cudaBlocks, cudaThreads;
		#endif

		pointFunction_t potential;
		void setPotential (pointFunction_t pp);

		std::vector <atom> atoms;
	private:
		float rc_, rs_;
		FILE *snapFile_;
		float3 box_;
		float targetT_, instantT_;
		float mass_;
		float Uk_, Up_;
		std::vector <float> potentialArgs_;
};

#endif
