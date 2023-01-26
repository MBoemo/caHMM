//----------------------------------------------------------
// Copyright 2023 University of Cambridge
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under the MIT license.
// You should have received a copy of the license with
// this software.  If not, please Email the author.
//----------------------------------------------------------

#ifndef PENTHUS_ERROR_HANDLING_H
#define PENTHUS_ERROR_HANDLING_H

#include <exception>
#include <string.h>
#include <string>

class NegativeLog : public std::exception {
	public:
		virtual const char * what () const throw () {
			return "Negative value passed to natural log function.";
		}
};


class DivideByZero : public std::exception {
	public:
		virtual const char * what () const throw () {
			return "lnQuot: Cannot divide by zero.";
		}
};


class NotFinalised : public std::exception {
	public:
		virtual const char * what () const throw () {
			return "Finalise HMM before running a learning algorithm.";
		}
};


class NumericalInstability : public std::exception {
	public:
		virtual const char * what () const throw () {
			return "Aborted HMM training due to numerical instability.";
		}
};

class BadIndex : public std::exception {
	public:
		virtual const char * what () const throw () {
			return "Matrix subscript out of bounds.";
		}
};

#endif
