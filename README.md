CBEMDGPU
========

CBE MD on GPUs

by Nathan A. Mahynski, George A. Khoury, and Carmeline J. Dsilva

See main.cpp to set parameters which are documented by example in this file.

To compile the CPU version, type 
$ make MD

To compile the GPU version, type
$ make -f Makefile_cuda

To compile tests, etc., type
$ make ${var}
where var is the suffix desired (see Makefile for details)