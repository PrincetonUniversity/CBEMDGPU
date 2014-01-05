/*!
 * Common values, structures, etc.
 * \date 11/17/13
 */

#ifndef __COMMON_H__
#define __COMMON_H__

#include <exception>
#include <string>

//! Error reporting
class customException : public std::exception {
	public:
		const char* what() const throw() {return msg_.c_str();} //!< Return a message pertaining to what went wrong
		customException (std::string m="custom exception occurred"):msg_(m){}   //!< Instantiate the object with a user specified message
		~customException () throw() {};
	private:
		std::string msg_;       //!< Message to report
};

//! Default OMP chunk size
#define OMP_CHUNK 100

#endif
