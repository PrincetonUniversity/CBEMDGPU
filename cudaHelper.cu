/*!
 * Functions to assist to identifying and characterizing CUDA devices available
 * \author Nathan A. Mahynski
 */

#include "cudaHelper.h"
#include <iostream>
using namespace std;

systemProps::systemProps (string fname) {
	output_ = new std::ofstream(fname.c_str());
}

systemProps::~systemProps () {
	if (output_ != &std::cout) {
		delete output_;
	}
}

/*!
 * Display properties of all CUDA capable devices
 */
__host__ void systemProps::displayCudaProps () {
	checkCudaDevices_();
	
	for (unsigned int i = 0; i < device_.size(); ++i) {
		*output_ << "# > Device Name: " <<  device_[i].name << endl;
		*output_ << "# -----------------------------------" << endl;
		*output_ << "# Total Global Memory: " << device_[i].totalGlobalMem/1024 << " KB" << endl;
		*output_ << "# Shared Memory available per Block: " << device_[i].sharedMemPerBlock/1024 << " KB" << endl;
		*output_ << "# Registers per Thread Block: " << device_[i].regsPerBlock << endl;
		*output_ << "# Warp Size: " << device_[i].warpSize << endl;
		*output_ << "# Memory Pitch: " << device_[i].memPitch << endl;
		*output_ << "# Maximum Threads per Block: " << device_[i].maxThreadsPerBlock << endl;
		*output_ << "# Maximum Thread Dimensions (Block): " << device_[i].maxThreadsDim[0] << "x" << device_[i].maxThreadsDim[1] << "x" << device_[i].maxThreadsDim[2] << endl;
		*output_ << "# Maximum Thread Dimensions (Grid): " << device_[i].maxGridSize[0] << "x" << device_[i].maxGridSize[1] << "x" << device_[i].maxGridSize[2] << endl;
		*output_ << "# Total Constant Memory: " << device_[i].totalConstMem << " bytes" << endl;
		*output_ << "# CUDA version: " << device_[i].major << "." << device_[i].minor << endl;
		*output_ << "# Clock Rate: " << device_[i].clockRate << " kHz" << endl;
		*output_ << "# Texture Alignment: " << device_[i].textureAlignment << " bytes"<< endl;
		if (device_[i].deviceOverlap == 0) {
			*output_ << "# Device Overlap: Not Allowed" << endl;
		} else {
			*output_ << "# Device Overlap: Allowed" << endl;
		}
		*output_ << "# Number of Multiprocessors: " << device_[i].multiProcessorCount << endl;
		*output_ << "# -----------------------------------"<< endl;
	}
}

/*!
 * Check CUDA devices and store all their properties
 */
__host__ void systemProps::checkCudaDevices_ () {
	int device_Count = 0;
	cudaError_t error = cudaGetDeviceCount(&device_Count);
	
	if (error != cudaSuccess) {
		*output_ << "No CUDA capable devices found." << endl;
		device_Count = 0;
	}
	
	device_.resize(device_Count);
	
	for (unsigned int i = 0; i < device_Count; ++i) {
		cudaGetDeviceProperties(&device_[i], i);
	}
}

/*!
 * Identify CPU host information
 */
void systemProps::checkHost_ () {
	struct utsname Uname;
	struct hostent *host;
	char ihname[64], ident[1024];
	
	if ((uname(&Uname) < 0) || (gethostname(ihname,64) != 0)) {
		hostname_ = "Unknown Hostname";
		return;
	}
	
	host=gethostbyname(ihname);
	strcpy(ident,Uname.sysname);
	strcat(ident," ");
	strcat(ident,Uname.release);
	strcat(ident," ");
	strcat(ident,Uname.machine);
	strcat(ident," [");
	strcat(ident,(*host).h_name);
	strcat(ident,"]");
	
	hostname_ = ident;
}
