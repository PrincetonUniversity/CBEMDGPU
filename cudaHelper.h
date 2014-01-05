/*!
 * Functions to assist in using CUDA and/or GPUs
 * \date 11/21/13
 */

#ifndef __CUDAHELPER_H__
#define __CUDAHELPER_H__

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

//! Holds all information about the CPU and GPU the simulation is being performed on.
class systemProps {
public:
	systemProps () {output_ = &std::cout;}
	systemProps (std::string fname);
	~systemProps ();
	void displayAllProps () {displayHost(); displayCudaProps();}    //!< Display cpu host and device (if using GPUs) properties
	__host__ void displayCudaProps ();                              //!< Display the properties of the GPU being used
	void displayHost () {checkHost_(); *output_ << "# " << hostname_ << std::endl;}     //!< Display host name and other properties
	int maxThreadsPerBlock (const int devID) {return device_[devID].maxThreadsPerBlock;}    //!< Returns the maximum number of threads per block on the GPU
	int maxGridDimX (const int devID) {return device_[devID].maxGridSize[0];}               //!< Returns the maximum number of blocks per grid on the GPU
	int numDevices () const {return device_.size();}                                        //!< Returns the number of GPUs found attached to the CPU
private:
	std::vector < cudaDeviceProp > device_;                         //!< CUDA properties for each GPU attached to the CPU
	std::string hostname_;                                          //!< Hostname
	std::ostream* output_;                                          //!< Output file stream (stdout or can be a user specified file)
	__host__ void checkCudaDevices_ ();                             //!< Check on each GPU
	void checkHost_ ();                                             //!< Check on the CPU
};

#endif
