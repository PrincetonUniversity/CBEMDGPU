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

