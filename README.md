CBEMDGPU
========

> CBE MD on GPUs

>> by Nathan A. Mahynski, George A. Khoury, and Carmeline J. Dsilva

See main.cpp to set parameters which are documented by example in this file.

To compile the CPU version, type 
$ make MD

To compile the GPU version, type
$ make -f Makefile_cuda

To compile tests, type
$ make TESTS

To compile the timing executable (used for scaling studies), type
$ make TIMING

To compile the program that runs a simulation to compare with LAMMPS output, type
$ make LMP_COMPARE

To compile the program that tests the NVE integrator, type
$ make TEST_NVE

In the Makefile, the PATHTOBOOST variable should point to the diretory where the C++ boost libraries are saved.

In the GTEST_DIR variable should point to the directory where the Google Tests libraries are saved.

Note that the Intel C++ compilers must be used; the GNU C++ compilers do not work with the OpenMP portion of our code (this is a known bug in the compiler).


Execution
====
main.cpp expects 1 input, the number of threads to use with OMP.
$ ./md numThreads > log 2> err

This also produces a trajectory.xyz file which can be visualized with VMD (if you have it installed)
$ vmd -xyz trajectory.xyz

Explanation of main.cpp
====
> How to use, change, and make your own in 10 steps.

1. main.cpp expects 1 input, the number of threads to use with OMP.
$ ./md numThreads > log 2> err

2. The random number generator seed is then set manually to ensure that results are reproducible.

3. The user must then specify basic properties about the system.  A systemDefinition object should be instantiated and then temperature, box size, particle mass, skin radius, and pair potential cutoff should be specified 
    
    $ systemDefinition a;
    
    $ a.setTemp(1.0); 
    
    $ a.setRskin(1.0); ...

4. The system can then be initialized with the command initThermal or initRandom (see doxygen documentation for more details).  For example, in the former:

    $ a.initThermal(nAtoms, Temp, rngSeed, separation);

5. Where "nAtoms" is the number of atoms, "separation" is the initial separation of the particles on the simple cubic lattice on which they are initialized, and the rest of the variables are self-explanatory.  This command completely initializes the system for simulation.

6. Next, the pair potential should be specified.  As an example in main.cpp, the preprocessor flag NVCC is used to select either the CPU version of the shifted lennard-jones potential (slj) or the GPU version (dev_slj).  The Makefile (compare to Makefile_cuda) will define this if necessary and allows the code to flow naturally and work in both cases.
    
    $ pointFunction_t pp = slj;
    
   $ a.setPotential(pp);

7. If the GPUs are used (NVCC compiler flag is defined) the number of threads and blocks the GPU kernel will be invoked with must be set.  This is done next.

8. Then additional arguments for the pair potential must be specified.  In the case of the slj function, variables epsilon, sigma, delta, and ushift must be given.  This is then set in the systemDefinition.

    $ std::vector <float> args(5);
    
    $ args[0] = 1.0; // epsilon
    
    $ args[1] = 1.0; // sigma
    
    $ args[2] = 0.0; // delta 
    
    $ args[3] = 0.0; // ushift
    
    $ a.setPotentialArgs(args);

9. Afterwards, the integrator (ensemble) must be specified. In the case of the NVT ensemble where we are using the Nose-Hoover thermostat, the integrator is given by the object nvt_NH.  This requires a damping constant, which for the slj potential should be about unity.  Then the numerical timestep should be given, for the slj this is usually efficient around 0.005.

    $ nvt_NH integrate (1.0);
    
    $ integrate.setTimestep(timestep);

10. Finally the simulation is ready to iterate.  A simple loop can be set to do this. An example for the case of the NVE ensemble is also provided in test_nve.cpp which can be compiled with make TEST_NVE (see test_nve.cpp)


FYI
====

> Known bugs, etc.


There are a few known instances of bugs related to compiler options, etc.

1. Using icpc instead of g++
	Our code uses OMP to parallelize many calculations.  Specifically, the atoms member (vector) in the systemDefinition object must be shared often.  
	However, because it is a member of a class g++ struggles to properly share this in memory.  
	This is a known bug in the g++ compiler which the intel (icpc) version handles rigorously. 
	As a result, our code will produce errors if compiled with the g++ compiler when more than 1 core is used.  
	To use this on tiger the module openmpi/intel-12.1/1.4.5/64 should be loaded to make the icpc compiler accessible.

2. Cuda toolkit 5.5
	CUDA is a finicky tool.  Different GPUs require different toolkits and versions to work properly.  
	In fact, compilation may succeed with a bad version but the run time behavior produces unexpected (incorrect) results.
	For the K20 cards on tiger, the latest toolkit (v5.5) must be loaded.
	To do so, load the module cudatoolkit/5.5.22 before attempting to compile the program.
	The Makefile must include flags consistent with the GPUs version of CUDA (which on tiger in 3.5) so the Makefile_cuda contains a flag "NVFLAGS = -gencode arch=compute_35,code=sm_35"  for the .cu files. 
	Furthermore, you will find that preprocessor flags NVCC and NOGPU are found throughout the code which act as switches to activate/deactivate GPU functionality throughout the compilation process.

3. Tiger GPUs
	Unfortunately it appears that the queueing system on tiger is having problems reserving all or some of the GPUs exclusively for single jobs by users.
	As a result, thrust (the GPU equivalent of the STL for C++) will have memory issues if one or more of the GPUs on a node are already in use.
	Because of the unusually high load on tiger over the past month, we have been unable to obtain good results on this cluster since usually this situation is encountered.  
	Our private GPU cluster was used to obtian our results instead, though our Makefile_cuda is set so this should compile properly on tiger if you want to check.
