/*!
 * Common values, structures, etc.
 * \author Nathan A. Mahynski
 * \date 11/17/13
 */

#ifndef __COMMON_H__
#define __COMMON_H__

#include <exception>
#include <string>

//! Error reporting
class customException : public std::exception {
	public:
		const char* what() const throw() {return msg_.c_str();}
		customException (std::string m="custom exception occurred"):msg_(m){}
		~customException () throw() {};
	private:
		std::string msg_;
};

#define OMP_CHUNK 100

#endif
