//----------------------------------------------------------
// Copyright 2023 University of Cambridge
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under the MIT license.
// You should have received a copy of the license with
// this software.  If not, please Email the author.
//----------------------------------------------------------

#include "states.h"

State::State( Distribution *d, std::string n, std::string m = "", std::string t = "", double w = 1.0) {
	meta = m;
	dist = d;
	name = n;
	weight = w;
	tether = t;
}
