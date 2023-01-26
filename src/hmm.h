//----------------------------------------------------------
// Copyright 2023 University of Cambridge
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under the MIT license.
// You should have received a copy of the license with
// this software.  If not, please Email the author.
//----------------------------------------------------------

#ifndef HMM_H
#define HMM_H

#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <utility>
#include <cassert>
#include <assert.h>
#include <limits>

#include "states.h"
#include "error_handling.h"

#define STARTNAME "START"
#define ENDNAME "END"

#pragma omp declare reduction(DoubleLnSum : double : \
omp_out = lnSum( omp_out, omp_in ) ) \
initializer(omp_priv = omp_orig)

#pragma omp declare reduction(vecDoubleLnSum : std::vector< double > : \
std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), lnSumInc< double >() )) \
initializer(omp_priv = omp_orig)


class HiddenMarkovModel {
	/*hidden Markov model object that houses states, transitions, and any machine learning algorithms we may want to run */
	private:
		bool finalised = false;
		unsigned int startIdx, endIdx;
        std::map< Distribution *, std::vector< int > > tethers;


	public:
        unsigned int silentStart;
		std::map< std::pair< std::string, std::string >, double > transitions;
		std::vector< State > states; /*TOCHANGE: For testing and development.  This should really be private. */
		std::vector< State > silentStates; /*TOCHANGE: For testing and development.  This should really be private. */

		/*start and end */
        State start = State( NULL, STARTNAME, "", "", 1.0 );
        State end = State( NULL, ENDNAME, "", "", 1.0 );

		/*aux functions */
		void add_state( State & );
		void add_transition( State &, State &, double );
		void finalise( );
		void sort_silentStates ( );
		void summarise( );

		/*machine learning functions */
		double sequenceProbability( std::vector <double> & );
		std::pair< double, std::vector< std::vector< double > > > forward( std::vector< double > & );
		std::pair< double, std::vector< std::vector< double > > > backward( std::vector< double > & );
		std::pair< double, std::vector< std::string > > viterbi( std::vector< double > & );
		void BaumWelch( std::vector< std::vector< double > > &, double = 1.0, int = 250, double = 0.0, bool = true, int = 1, bool = true ); 
};

#endif
