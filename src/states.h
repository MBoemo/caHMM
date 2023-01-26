//----------------------------------------------------------
// Copyright 2023 University of Cambridge
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under the MIT license.
// You should have received a copy of the license with
// this software.  If not, please Email the author.
//----------------------------------------------------------

#ifndef STATES_H
#define STATES_H

#include<string>
#include<vector>
#include "probability.h"
#include "error_handling.h"

class State {
	/*simple class for a state object that houses a state's distribution, name, and weight */

	public:
        Distribution *dist;
		std::string name;
       	std::string tether;
		double weight;
		std::vector< unsigned int > inTransitions;
		std::vector< unsigned int > outTransitions;
		std::string meta; //string for some meta information about what this state does (if you don't want to house it in the name)
		State( Distribution *, std::string, std::string, std::string, double );	
};

#endif
