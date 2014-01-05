CBEMDGPU
========

CBE MD on GPUs

by Nathan A. Mahynski, George A. Khoury, and Carmeline J. Dsilva

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