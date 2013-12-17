/*
 * \author Nathan A. Mahynski
 */ 

#ifndef __cudaHelper__
#define __cudaHelper__

#include <cuda.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <string.h>
#include <netdb.h>
#include <sys/utsname.h>

#if defined(__STRICT_ANSI__) && !defined(__alpha__)
extern int gethostname __P ((char *__name, size_t __len));
#else
#include <unistd.h>
#endif


class systemProps {
public:
	systemProps () {output_ = &std::cout;}
	systemProps (std::string fname);
	~systemProps ();
	void displayAllProps () {displayHost(); displayCudaProps();}
	__host__ void displayCudaProps ();
	void displayHost () {checkHost_(); *output_ << "# " << hostname_ << std::endl;}
	int maxThreadsPerBlock (const int devID) {return device_[devID].maxThreadsPerBlock;}
	int maxGridDimX (const int devID) {return device_[devID].maxGridSize[0];}
	int numDevices () const {return device_.size();}
private:
	std::vector < cudaDeviceProp > device_;
	std::string hostname_;
	std::ostream* output_;
	__host__ void checkCudaDevices_ ();
	void checkHost_ ();
};

#endif
