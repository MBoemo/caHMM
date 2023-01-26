//----------------------------------------------------------
// Copyright 2023 University of Cambridge
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under the MIT license.
// You should have received a copy of the license with
// this software.  If not, please Email the author.
//----------------------------------------------------------

#ifndef PROBABILITY_H
#define PROBABILITY_H

#include <iostream>
#include <math.h>
#include <algorithm>
#include "error_handling.h"

#define _USE_MATH_DEFINES

/*
Functions for computing probabilities in log space to improve numerical stability.  Most algorithms 
inspired by Numerically Stable Hidden Markov Model Implementation by Tobias P. Mann
(see http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf).
*/

inline double eexp( double x ){
/*Map from log space back to linear space. */

	if ( std::isnan( x ) ) {
		return 0.0;
	}	
	else {
		return exp( x );
	}
}


inline double eln( double x ){
/*Map from linear space to log space. */    

	if (x == 0.0){
		return NAN;
	}
	else if (x > 0.0){
		return log( x );
	}
	else{
		throw NegativeLog();
	}
}


inline double lnSum( double ln_x, double ln_y ){
/*evalutes the quotient ln_x + ln_y */    
    
	if ( std::isnan( ln_x ) || std::isnan( ln_y ) ){

        	if ( std::isnan( ln_x ) && std::isnan( ln_y ) ){
			return NAN;
		}
		else if ( std::isnan( ln_x ) ){
			return ln_y;
		}
		else{
			return ln_x;
		}
	}
	else{
	
		if ( ln_x > ln_y ){
			return ln_x + eln( 1.0 + eexp( ln_y - ln_x ) );
		}
		else{
			return ln_y + eln( 1.0 + eexp( ln_x - ln_y ) );
		}
	}
}


inline double lnProd( double ln_x, double ln_y ){
/*evalutes the quotient ln_x*ln_y */

	if ( std::isnan( ln_x ) || std::isnan( ln_y ) ){
		return NAN;
	}
	else{
		return ln_x + ln_y;
	}
}


inline double lnQuot( double ln_x, double ln_y ){
/*evalutes the quotient ln_x/ln_y */

	if ( std::isnan( ln_y ) ){
        	throw DivideByZero();
	}
	else if ( std::isnan( ln_x ) ){
		return NAN;
	}
	else{
		return ln_x - ln_y;
	}
}


inline bool lnGreaterThan( double ln_x, double ln_y ){
/*evalutes whether ln_x is greater than ln_y, and returns a boolean */
    	
	if ( std::isnan( ln_x ) || std::isnan( ln_y ) ){
		
		if ( std::isnan( ln_x ) || std::isnan( ln_y ) == false ){
			return false;
		}
        	else if ( std::isnan( ln_x ) == false || std::isnan( ln_y ) ){
			return true;
		}
		else{
			return false;
		}
	}
	else{
		if ( ln_x > ln_y ){
			return true;
		}
		else{
			return false;
		}
	}
}


/*lnSum functor for openMP */
template <class T> struct lnSumInc {

	T operator() ( const T& x, const T& y ) const {return lnSum( x, y );}
	typedef T first_argument_type;
	typedef T second_argument_type;
	typedef T result_type;
};


class Distribution{

	protected:
		Distribution( double x, double y ){
			param1 = x;
			param2 = y;
			trained_param1 = x;
			trained_param2 = y;
		}
	public:
		double param1, param2, trained_param1, trained_param2;
		virtual double pdf( double ) = 0;
		virtual std::string identify( void ) = 0;
};


class NormalDistribution: public Distribution {
	
	public:
		NormalDistribution( double x, double y ): Distribution( x, y ) {}

		inline double pdf( double x ){
			return ( 1.0/sqrt( 2.0*pow( trained_param2, 2.0 )*M_PI ) )*exp( -pow( x - trained_param1 , 2.0 )/( 2.0*pow( trained_param2, 2.0 ) ) );
		}
		std::string identify( void ){
			return "NormalDistribution";
		}
};


class UniformDistribution: public Distribution{

	public:
		UniformDistribution( double x, double y ): Distribution( x, y ) {}

		inline double pdf( double x){
			double probability;	
	
			if ( x >= param1 && x <= param2 ){
				probability = 1.0/( trained_param2 - trained_param1 );
			}
			else {
				probability = 0.0;
			}
			return probability;
		}
		std::string identify( void ){
			return "UniformDistribution";
		}
};


class SilentDistribution: public Distribution{

	public:
		SilentDistribution( double x, double y ): Distribution( x, y ) {}

		inline double pdf( double x){
			std::cout << "Exited with error.  Cannot calculate the probability density function of a silent state - it has no distribution." << std::endl;
			exit( EXIT_FAILURE );
			return 0;
		}
		std::string identify( void ){
			return "SilentDistribution";
		}
};

#endif
