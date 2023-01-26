//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------

//Compile: g++ -o simpleHMM simpleHMM.cpp -std=c++11 -lPenthus -L../../Penthus -fopenmp


#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <unistd.h>
#include <libgen.h>
#include "../src/hmm.h"
#include "../src/states.h"
#include "../src/probability.h"

struct PoreParameters {

	double shift;
	double drift = 0.0;
	double scale;
	double var = 1.0;
};

//Initial transitions within modules (internal transitions)
static double internalM12I = 0.3475;
static double internalI2I = 0.5;
static double internalM12M1 = 0.4;

//Initial transitions between modules (external transitions)
static double externalD2D = 0.3;
static double externalD2M1 = 0.7;
static double externalI2M1 = 0.5;
static double externalM12D = 0.0025;
static double externalM12M1 = 0.25;

std::string getExePath(void){

	int PATH_MAX=1000;
	char result[ PATH_MAX ];
	ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
	const char *path;

	if (count != -1) path = dirname(result);

	std::string s(path);
	return s;
}

std::map< std::string, std::pair< double, double > > import_poreModel( std::string poreModelFilename ){

	std::string pathExe = getExePath();
	std::string modelPath = pathExe + "/"+ poreModelFilename;
	//std::cout << modelPath << std::endl;

	/*map that sends a 5mer or 6mer to the characteristic mean and standard deviation (a pair) */
	std::map< std::string, std::pair< double, double > > kmer2MeanStd;

	/*file handle, and delimiter between columns (a \t character in the case of ONT model files) */
	std::ifstream file( modelPath );

	if ( not file.is_open() ) std::cout << "could not open file" << std::endl;

	std::string line, key, mean, std;
	std::string delim = "\t";

	/*while there's a line to read */
	while ( std::getline( file, line ) ){

		/*and the line isn't part of the header */
		if ( line.substr(0,4) != "kmer" && line[0] != '#' ){ 

			/*the kmer, mean, and standard deviation are the first, second, and third columns, respectively. */
			/*take the line up to the delimiter (\t), erase that bit, and then move to the next one */
			key = line.substr( 0, line.find( delim ) );
			line.erase( 0, line.find( delim ) + delim.length() );

			mean = line.substr( 0, line.find( delim ) );
			line.erase( 0, line.find( delim ) + delim.length() );

			std = line.substr( 0, line.find( delim ) );

			/*key the map by the kmer, and convert the mean and std strings to doubles */
			kmer2MeanStd[ key ] = std::make_pair( atof( mean.c_str() ), atof( std.c_str() ) );
		}
	}

	/*if you need to print out the map for debugging purposes 
	for(auto it = kmer2MeanStd.cbegin(); it != kmer2MeanStd.cend(); ++it){
	    std::cout << it->first << " " << it->second.first << " " << it->second.second << "\n";
	}*/
	
	return kmer2MeanStd;
}

double sequenceProbability_Penthus( std::vector <double> &observations,
				std::string &sequence, 
				size_t windowSize, 
				bool useBrdU, 
				PoreParameters scalings,
				size_t BrdUStart,
				size_t BrdUEnd,
				std::map< std::string, std::pair< double, double > > &thymidineModel ){

	HiddenMarkovModel hmm = HiddenMarkovModel();

	/*STATES - vector (of vectors) to hold the states at each position on the reference - fill with dummy values */
	std::vector< std::vector< State > > states( 3, std::vector< State >( sequence.length() - 5, State( NULL, "", "", "", 1.0 ) ) );

	/*DISTRIBUTIONS - vector to hold normal distributions, a single uniform and silent distribution to use for everything else */
	std::vector< NormalDistribution > nd;
	nd.reserve( sequence.length() - 5 );

	SilentDistribution sd( 0, 0 );
	UniformDistribution ud( 0, 250.0 );

	std::string loc, sixMer;
		
	/*create make normal distributions for each reference position using the ONT 6mer model */
	for ( unsigned int i = 0; i < sequence.length() - 5; i++ ){

		sixMer = sequence.substr( i, 6 );
		//std::cout << sixMer << std::endl;
		//std::cout << thymidineModel.at(sixMer).first << " " << thymidineModel.at(sixMer).second << std::endl;
		//std::cout << scalings.shift + scalings.scale * thymidineModel.at(sixMer).first << " " << scalings.var * thymidineModel.at(sixMer).second << std::endl;

		nd.push_back( NormalDistribution( scalings.shift + scalings.scale * thymidineModel.at(sixMer).first, scalings.var * thymidineModel.at(sixMer).second ) );

	}

	/*the first insertion state after start */
	State firstI = State( &ud, "-1_I", "", "", 1.0 );
	hmm.add_state( firstI );

	/*add states to the model, handle internal module transitions */
	for ( unsigned int i = 0; i < sequence.length() - 5; i++ ){

		loc = std::to_string( i );
		sixMer = sequence.substr( i, 6 );

		states[ 0 ][ i ] = State( &sd,		loc + "_D", 	sixMer,	"", 		1.0 );		
		states[ 1 ][ i ] = State( &ud,		loc + "_I", 	sixMer,	"", 		1.0 );
		states[ 2 ][ i ] = State( &nd[i], 	loc + "_M1", 	sixMer,	loc + "_match", 1.0 );

		/*add state to the model */
		for ( unsigned int j = 0; j < 3; j++ ){

			states[ j ][ i ].meta = sixMer;
			hmm.add_state( states[ j ][ i ] );
		}

		/*transitions between states, internal to a single base */
		/*from I */
		hmm.add_transition( states[1][i], states[1][i], internalI2I );

		/*from M1 */
		hmm.add_transition( states[2][i], states[2][i], internalM12M1 );
		hmm.add_transition( states[2][i], states[1][i], internalM12I );
	}

	/*add transitions between modules (external transitions) */
	for ( unsigned int i = 0; i < sequence.length() - 6; i++ ){

		/*from D */
		hmm.add_transition( states[0][i], states[0][i + 1], externalD2D );
		hmm.add_transition( states[0][i], states[2][i + 1], externalD2M1 );

		/*from I */
		hmm.add_transition( states[1][i], states[2][i + 1], externalI2M1 );

		/*from M */
		hmm.add_transition( states[2][i], states[0][i + 1], externalM12D );
		hmm.add_transition( states[2][i], states[2][i + 1], externalM12M1 );
	}

	/*handle start states */
	hmm.add_transition( hmm.start, firstI, 0.25 );
	hmm.add_transition( hmm.start, states[0][0], 0.25 );
	hmm.add_transition( hmm.start, states[2][0], 0.5 );

	/*transitions from first insertion */
	hmm.add_transition( firstI, firstI, 0.25 );
	hmm.add_transition( firstI, states[0][0], 0.25 );
	hmm.add_transition( firstI, states[2][0], 0.5 );

	/*handle end states */
	hmm.add_transition( states[0][sequence.length() - 6], hmm.end, 1.0 );
	hmm.add_transition( states[1][sequence.length() - 6], hmm.end, externalI2M1 );
	hmm.add_transition( states[2][sequence.length() - 6], hmm.end, externalM12M1 + externalM12D );

	hmm.finalise();

	return hmm.sequenceProbability( observations );
}


int main(){

	//make some simple test data - uncomment one or the other and check against pomegranateCrosscheck.py
	
	std::string sequence = "CTCCACTCGTTACCCTGTCCCATTCA";
	std::vector<double> observations = {98.5804, 94.1981, 85.6842, 87.9788, 91.598, 109.494, 86.048, 72.0098, 84.7888, 97.4931, 95.5623, 100.11, 91.7845, 86.4678, 92.8339, 88.0668, 96.8873, 95.0026, 94.0932, 97.2692, 92.5841, 101.299, 83.0631, 98.8736, 109.134, 105.496, 95.7022, 97.3252, 94.303, 73.0359, 70.9371, 86.8595, 89.7138, 82.4568, 87.727, 93.1837, 87.3073};
	
	/*
	std::string sequence = "CACACACACATCCTAACACTACCCTA";
	std::vector<double> observations = {96.6349, 95.5902, 94.0932, 102.325, 99.7597, 102.418, 95.3664, 92.0644, 94.1911, 95.2824, 97.1946, 101.299, 104.563, 103.607, 100.627, 101.858, 97.8476, 103.178, 95.8621, 102.884, 96.6816, 99.0602, 95.8701, 100.646, 97.5677, 82.3502, 79.0056, 81.4868, 116.06, 107.175, 96.7282, 96.8215, 112.492, 110.487, 114.031, 105.776, 91.7845, 88.893, 93.9113, 92.0644, 100.179, 92.3442, 92.624, 95.2824, 96.3551, 101.802, 112.492, 108.994, 93.8633, 93.6034, 103.864, 99.7877, 97.5677, 100.515, 102.511};
	*/	

	//both the above samples use the same scalings
	PoreParameters scalings;
	scalings.shift = 4.55653;
	scalings.scale = 0.956433;
	scalings.var = 1.75943;

	//
	size_t windowSize = 10;
	size_t BrdUStart = windowSize;
	size_t BrdUEnd = sequence.substr(windowSize,6).rfind('T') + windowSize;

	std::map< std::string, std::pair< double, double > > poreModel = import_poreModel( "template_median68pA.model" );

	double prob = sequenceProbability_Penthus( observations, sequence,  windowSize,  false,  scalings, BrdUStart, BrdUEnd, poreModel );
	std::cout << "Log probability: " << prob << std::endl;
	return 0;
}

