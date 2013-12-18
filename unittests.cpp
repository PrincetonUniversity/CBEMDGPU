#include "system.h"
#include "potential.h"
#include "integrator.h"
#include "nvt.h"
#include <iostream>
#include "utils.h"
#include <omp.h>
#include <stdlib.h>
#include "gtest/gtest.h"

class SystemTest : public ::testing::Test {
protected:
	
	SystemTest() : integrate(1.0) {}

	virtual void SetUp() {
		
		natoms = 2;
		mass = 1.0;
		L = 12.0;	
		eps = 1.0;
		sigma = 1.0;

		const int rngSeed = 3145;
		float Temp = 0.5; 
		float rcut = 2.5;
		
		a.setBox(L, L, L);
		a.setTemp(Temp);
		a.setMass(mass); 
		a.setRskin(1.0);
		a.setRcut(rcut);	
		a.initThermal(natoms, Temp, rngSeed, 1.2);
		
		pointFunction_t pp = slj;
		a.setPotential(pp);
		std::vector <float> args(5);
		args[0] = eps; // epsilon
		args[1] = sigma; // sigma
		args[2] = 0.0; // delta 
		args[3] = 0.0; // ushift
		a.setPotentialArgs(args);

		float timestep = 0.002;
		integrate.setTimestep(timestep);
		integrate.step(a);
	}

    	systemDefinition a;
	int natoms;
	float mass;
	float L;
	float eps;
	float sigma;
	nvt_NH integrate;
};	
	
TEST_F(SystemTest, NumAtoms) {
	ASSERT_EQ(natoms, a.numAtoms());
}

TEST_F(SystemTest, KineticEnergy) {
	atom atom1 = a.atoms[0];
	atom atom2 = a.atoms[1];
	
	float ke = 0.5 * mass * (atom1.vel.x*atom1.vel.x + atom1.vel.y*atom1.vel.y + atom1.vel.z*atom1.vel.z + atom2.vel.x*atom2.vel.x + atom2.vel.y*atom2.vel.y + atom2.vel.z*atom2.vel.z);
	
	ASSERT_FLOAT_EQ(ke, a.KinE());
}

	
TEST_F(SystemTest, PotentialEnergy) {
	float r = 0.5;
 
	float sig_r2 = (sigma / r) * (sigma / r);
	float sig_r6 = sig_r2 * sig_r2 * sig_r2;
	float sig_r12 = sig_r6 * sig_r6;
 
	a.atoms[0].pos.x = 0.0;
	a.atoms[0].pos.y = 0.0;
	a.atoms[0].pos.z = 0.0;

	a.atoms[1].pos.x = 0.0;
	a.atoms[1].pos.y = 0.0;
	a.atoms[1].pos.z = r;

	integrate.calcForce(a);

	float pe = 4.0 * eps * (sig_r12 - sig_r6);
	
	ASSERT_FLOAT_EQ(pe, a.PotE());
}

TEST_F(SystemTest, ChangeBox) {
	
	float3 box1 = a.box();
	
	ASSERT_FLOAT_EQ(L, box1.x);
	ASSERT_FLOAT_EQ(L, box1.y);
	ASSERT_FLOAT_EQ(L, box1.z);

	float L2 = 6.0;
	a.setBox(L2, L2, L2);
	box1 = a.box();
	
	ASSERT_FLOAT_EQ(L2, box1.x);
	ASSERT_FLOAT_EQ(L2, box1.y);
	ASSERT_FLOAT_EQ(L2, box1.z);

}

TEST_F(SystemTest, PBC) {
	float r = 0.5;
 
	a.atoms[0].pos.x = 0.0;
	a.atoms[0].pos.y = 0.0;
	a.atoms[0].pos.z = 0.0;

	a.atoms[1].pos.x = 0.0;
	a.atoms[1].pos.y = 0.0;
	a.atoms[1].pos.z = r;

	float3 dr;
	ASSERT_FLOAT_EQ(r * r, pbcDist2(a.atoms[0].pos, a.atoms[1].pos, dr, a.box()));

	a.atoms[1].pos.z = L - r;
	ASSERT_FLOAT_EQ(r * r, pbcDist2(a.atoms[0].pos, a.atoms[1].pos, dr, a.box()));
	
}

int main (int argc, char** argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();

}
