//----------------------------------------------------------
// Copyright 2023 University of Cambridge
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under the MIT license.
// You should have received a copy of the license with
// this software.  If not, please Email the author.
//----------------------------------------------------------

//#define TEST_SORTING 1
//#define TEST_FORWARD 1
//#define TEST_BACKWARD 1
//#define TEST_BAUMWELCH 1


#include "hmm.h"

//-------------------------------------------------------------------------------------------------------------------------------
void HiddenMarkovModel::add_state( State &x ) {
	/* Add a state to the HMM.  States must have unique names.
	ARGUMENTS
	---------
	- x: state object
	  type: State  */

	if ( x.name == STARTNAME || x.name == ENDNAME ){

		std::cout << "Exiting with error.  Added state has the same name as the protected start or end state." << std::endl;
		exit( EXIT_FAILURE );
	}

	for ( unsigned int i = 0; i < states.size(); i++ ){
		
		if ( x.name == states[i].name ){
			
			std::cout << "Exiting with error.  Added state has the same name as a previously added state.  Names must be unique." << std::endl;
			exit( EXIT_FAILURE );
		}
	}

	if ( x.dist -> identify() == "SilentDistribution" ){
		silentStates.push_back( x );
	}
	else {
		states.push_back( x );
	}
}


//-------------------------------------------------------------------------------------------------------------------------------
void HiddenMarkovModel::add_transition( State &x, State &y, double z ) {
	/* Add a transition between states to the HMM.  The states must have been already added to the HMM, and they must have unique names.
	ARGUMENTS
	---------
	- x: state the transition is going from
	  type: State
	- y: state the transition is going to
	  type: State
	- z: probability of the transition (must be between 0 and 1)
	  type: double */ 

	std::pair< std::string, std::string > key = make_pair( x.name, y.name );
	transitions[ key ] = z;
}


//-------------------------------------------------------------------------------------------------------------------------------
void HiddenMarkovModel::finalise( ){
	/* Run after all states/transitions are added.  Sorts silent states so that transitions to from low to high index. */

	/* add start and end states to the model */
	silentStates.push_back( start );
	silentStates.push_back( end );

	/*run a topological sort on the silent states so that all transitions between silent states go from low index to high index.*/
	HiddenMarkovModel::sort_silentStates( );

	/*index in the states vector where the silent states will begin */
	silentStart = states.size();

	/*append the (sorted) silent states vector to the end of the states vector */
	states.insert( states.end(), silentStates.begin(), silentStates.end() );

	/*Floating point tolerance - make sure that |1 - out_transitions| is less than tol. */
	double tol = 1e-9;

	/*after sorting, find the start and end states (wherever they landed) and set their indices to use later */
	for (unsigned int i = 0; i < states.size(); i++){

		if (states[i].name == STARTNAME){

			startIdx = i;
		}
		else if (states[i].name == ENDNAME){

			endIdx = i;
		}
	}
    
	/*verify that all states have out probabilities that sum to 1
	  while doing so, set the in states and out states for each state i */
	double total;
	std::pair< std::string, std::string > key;

	for ( unsigned int i = 0; i < states.size(); i++ ){

		total = 0.0;

		for ( unsigned int j = 0; j < states.size(); j++ ){

			key = make_pair( states[ i ].name, states[ j ].name );

			if ( transitions.count( key ) > 0 ){

				states[ i ].outTransitions.push_back( j ); //update outTransitions with j
				total += transitions[ key ]; //update total out transitions
			}

			key = make_pair( states[ j ].name, states[ i ].name );

			if ( transitions.count( key ) > 0 ){

				states[ i ].inTransitions.push_back( j ); //update outTransitions with j
			}
		}

		/*throw an error if the total out transition probabilities of any state (except the end state) don't sum to 1.*/
		if ( i != endIdx ){

			if (std::abs( 1 - total ) > tol){

				std::cout << "Exiting with error.  Out transitions for " << states[i].name << " do not sum to 1." << std::endl;
				exit( EXIT_FAILURE );
			}
		}
	}
    
	/*Go through the emitting states. If they have the same underlying distribution object, then they are tethered and the distributions should update together. */
	Distribution *d;
	for ( unsigned int i = 0; i < silentStart; i++){
    
		d = states[i].dist;
		tethers[ d ].push_back( i );
	}

	finalised = true;
}


//-------------------------------------------------------------------------------------------------------------------------------
void HiddenMarkovModel::sort_silentStates( ) {
	/*uses Kahn's algorithm to perform a topological sort on the silent states */

	bool pushS;
	std::pair< std::string, std::string > key;

	std::map< std::pair< std::string, std::string >, double > transitionsTemp = transitions;

	/*final vector of sorted silent states */
	std::vector< State > L;
	
	/*INITIALISATION: build a vector of nodes with no incoming edges. */
	std::vector< State > S;

	for ( unsigned int i = 0; i < silentStates.size(); i++ ) {

		pushS = true;

		/*and each silent state j that state i could transition to... */
		for ( unsigned int j = 0; j < silentStates.size(); j++ ) {

			key = make_pair( silentStates[ j ].name, silentStates[ i ].name );
				
				/*If node i has an incoming edge, we don't want this one - pass to the next loop iteration. */
				if ( transitionsTemp.count( key ) > 0 ){
					pushS = false;
					continue;
				}
		}

        /*If we never found a state j that has an incoming transition to state i, then we never flipped the boolean pushS.  So push state i onto S. */
		if (pushS == true ) S.push_back( silentStates[i] );
	}

	/*While there are still states in vector S with no incomming transitions */
	while ( S.size() > 0 ) {

		State n = S[ 0 ];
		L.push_back( n );
		S.erase( S.begin() );

		for ( unsigned int i = 0; i < silentStates.size(); i++ ) {
			
			key = make_pair( n.name, silentStates[ i ].name );

			if ( transitionsTemp.count( key ) > 0 ) {

				transitionsTemp.erase( key );
				int edgeCount = 0;

				for ( unsigned int j = 0; j < silentStates.size(); j++ ) {
					
					key = make_pair( silentStates[ j ].name, silentStates[ i ].name );

					if (transitionsTemp.count( key ) > 0) {
						edgeCount++;
					}
				}
				if ( edgeCount == 0 ) {
					S.push_back( silentStates[ i ] );
				}				
			}
		}
	}

	#if TEST_SORTING
	for ( unsigned int i = 0; i < L.size(); i++){

		for (unsigned int j = 0; j < L.size(); j++){

			key = make_pair( L[ i ].name, L[ j ].name );
			if (transitions.count(key) > 0) assert(i<j);
		}
	}

	#endif

	silentStates = L;
}


//-------------------------------------------------------------------------------------------------------------------------------
void HiddenMarkovModel::summarise(){
	/* returns human readable results of training */

	std::cout << "state\tinfo\toriMu\ttrMu\toriSig\ttrSig" << std::endl;

	for (unsigned int i = 0; i < silentStart; i++ ){

		if ( states[ i ].dist -> identify() == "NormalDistribution" ){
		
			if ( states[ i ].name.substr( states[ i ].name.length() - 2 ) == "M2" ) continue;

			std::cout << states[ i ].name << "\t" << states[ i ].meta << "\t" << states[ i ].dist -> param1 << "\t" << states[ i ].dist -> trained_param1 << "\t" << states[ i ].dist -> param2 << "\t" << states[ i ].dist -> trained_param2 << std::endl;
		}
	}
}


//-------------------------------------------------------------------------------------------------------------------------------
double HiddenMarkovModel::sequenceProbability( std::vector <double> &observations ) {
	/* Calculates the probability of a sequence of observations
	ARGUMENTS
	---------
	- observations: the observations on which to train the HMM
	  type: std::vector
	RETURNS
	-------
	- sequence probability
	  type: double */

	if ( not finalised ) throw NotFinalised();

	std::vector< std::vector< double > > forward( states.size(), std::vector< double >( observations.size() + 1, NAN ) );
	std::pair< std::string, std::string > key;
	double logProbability, pdfEval;
	unsigned int colIdx = 0;
	unsigned int k;


	/*-----------INITIALISATION----------- */
	/*starting state: set that equal to 1 and everything else to 0 */
	forward[ startIdx ][ 0 ] = 0.0;

	/*transitions between silent states before we emit the first symbol */
	for ( unsigned int i = silentStart; i < states.size(); i++ ){

		if ( i != startIdx ){
			logProbability = NAN;

			/*for all silent states... */
			for ( unsigned int j = 0; j < states[ i ].inTransitions.size(); j++ ){

				k = states[ i ].inTransitions[ j ];
				if ( k < silentStart ) continue; //only doing this for silent states

				key = make_pair( states[ k ].name, states[ i ].name );

				logProbability = lnSum( logProbability,
                                           		lnProd( forward[ k ][ 0 ], \
                                                   		eln (transitions.at( key ) ) ) );
			}
			forward[ i ][ 0 ] = logProbability;
		}
	}

	/*-----------RECURSION----------- */
	for ( unsigned int t = 0; t < observations.size(); t++ ){

		colIdx = t + 1;

		/*emitting states */
		for ( unsigned int i = 0; i < silentStart; i++ ){ 

			pdfEval = eln( states[ i ].dist -> pdf( observations[ t ] ) );

			logProbability = NAN;

			for ( unsigned int j = 0; j < states[ i ].inTransitions.size(); j++ ){

				k = states[ i ].inTransitions[ j ];

				key = make_pair( states[ k ].name, states[ i ].name );

				logProbability = lnSum( logProbability, \
                                            		lnProd( forward[ k ][ colIdx - 1 ], \
                                                    		lnProd( eln ( transitions.at( key ) ), \
                                                            			pdfEval ) ) );
			}
			forward[ i ][ colIdx ] = logProbability;
		}
		/*silent states */
		for ( unsigned int i = silentStart; i < states.size(); i++ ){
			logProbability = NAN;

			/*for all states j that transition into silent state i */
			for ( unsigned int j = 0; j < states[ i ].inTransitions.size(); j++ ){
				
				k = states[ i ].inTransitions[ j ];

				key = make_pair( states[ k ].name, states[ i ].name );

				logProbability = lnSum( logProbability, \
                                            		lnProd( forward[ k ][ colIdx ], \
                                                   		eln( transitions.at( key ) ) ) );
			}
			forward[ i ][ colIdx ] = logProbability;
		}
	}

	/*-----------TERMINATION----------- */
	double forwardProb = NAN;

	for ( unsigned int j = 0; j < states[ endIdx ].inTransitions.size(); j++ ){

		k = states[ endIdx ].inTransitions[ j ];

		key = make_pair( states[ k ].name ,ENDNAME );	

		forwardProb = lnSum( forwardProb, \
				     lnProd( forward[ k ][ colIdx ], \
					     eln( transitions.at( key ) ) ) );
	}
	return forwardProb;
}


//-------------------------------------------------------------------------------------------------------------------------------
std::pair< double, std::vector< std::vector< double > > > HiddenMarkovModel::forward( std::vector <double> &observations ) {
	/* The forward algorithm for HMM learning.
	ARGUMENTS
	---------
	- observations: the observations on which to train the HMM
	  type: std::vector
	RETURNS
	-------
	- a tuple that gives the forward probability and the forward algorithm matrix
	  type: std::pair */

	if ( not finalised ) throw NotFinalised();

	std::vector< std::vector< double > > forward( states.size(), std::vector< double >( observations.size() + 1, NAN ) );
	std::pair< std::string, std::string > key;
	double logProbability, pdfEval;
	unsigned int colIdx = 0;
	unsigned int k;

	/*-----------INITIALISATION----------- */
	/*starting state: set that equal to 1 and everything else to 0 */
	forward[ startIdx ][ 0 ] = 0.0;

	/*transitions between silent states before we emit the first symbol */
	for ( unsigned int i = silentStart; i < states.size(); i++ ){

		if ( i != startIdx ){
			logProbability = NAN;

			/*for all silent states... */
			for ( unsigned int j = 0; j < states[ i ].inTransitions.size(); j++ ){

				k = states[ i ].inTransitions[ j ];
				if ( k < silentStart or k >= i ) continue;

				key = make_pair( states[ k ].name, states[ i ].name );

				logProbability = lnSum( logProbability,
                                           		lnProd( forward[ k ][ 0 ], \
                                                   		eln (transitions.at( key ) ) ) );
			}
			forward[ i ][ 0 ] = logProbability;
		}
	}

	/*-----------RECURSION----------- */
	for ( unsigned int t = 0; t < observations.size(); t++ ){

		colIdx = t + 1;

		/*emitting states */
		for ( unsigned int i = 0; i < silentStart; i++ ){ 

			pdfEval = eln( states[ i ].dist -> pdf( observations[ t ] ) );

			logProbability = NAN;
						
			for ( unsigned int j = 0; j < states[ i ].inTransitions.size(); j++ ){

				k = states[ i ].inTransitions[ j ];

				key = make_pair( states[ k ].name, states[ i ].name );

				logProbability = lnSum( logProbability, \
                                            		lnProd( forward[ k ][ colIdx - 1 ], \
                                                    		eln ( transitions.at( key ) ) ) );
			}
			forward[ i ][ colIdx ] = lnProd(logProbability, pdfEval);
		}

		/*silent states - first pass*/
		for ( unsigned int i = silentStart; i < states.size(); i++ ){
			logProbability = NAN;

			/*for all states j that transition into silent state i */
			for ( unsigned int j = 0; j < states[ i ].inTransitions.size(); j++ ){
				
				k = states[ i ].inTransitions[ j ];
				if ( k >= silentStart ) continue;

				key = make_pair( states[ k ].name, states[ i ].name );

				logProbability = lnSum( logProbability, \
                                            		lnProd( forward[ k ][ colIdx ], \
                                                   		eln( transitions.at( key ) ) ) );
			}
			forward[ i ][ colIdx ] = logProbability;
		}

		/*silent states - second pass */
		for ( unsigned int i = silentStart; i < states.size(); i++ ){
			logProbability = NAN;

			/*for all states j that transition into silent state i */
			for ( unsigned int j = 0; j < states[ i ].inTransitions.size(); j++ ){
				
				k = states[ i ].inTransitions[ j ];
				if ( k < silentStart or k >= i ) continue;

				key = make_pair( states[ k ].name, states[ i ].name );

				logProbability = lnSum( logProbability, \
                                            		lnProd( forward[ k ][ colIdx ], \
                                                   		eln( transitions.at( key ) ) ) );
			}
			forward[ i ][ colIdx ] = lnSum(forward[ i ][ colIdx ], logProbability);
		}			

	}

	/*-----------TERMINATION----------- */
	double forwardProb = NAN;

	for ( unsigned int j = 0; j < states[ endIdx ].inTransitions.size(); j++ ){

		k = states[ endIdx ].inTransitions[ j ];

		key = make_pair( states[ k ].name ,ENDNAME );	

		forwardProb = lnSum( forwardProb, \
				     lnProd( forward[ k ][ colIdx ], \
					     eln( transitions.at( key ) ) ) );
	}

	#if TEST_FORWARD
	for (int j = 0; j < observations.size() + 1; j++){
		for (int i = 0; i < states.size(); i++){
			std::cout << forward[i][j] << " ";
		}
		std::cout << "]"<< std::endl;
	}
	#endif

	std::pair< double, std::vector< std::vector< double > > > forwardPair = std::make_pair( forwardProb, forward );

	return forwardPair;
}


//-------------------------------------------------------------------------------------------------------------------------------
std::pair< double, std::vector< std::vector< double > > > HiddenMarkovModel::backward( std::vector <double> &observations ) {
	/* The backward algorithm for HMM learning.
	ARGUMENTS
	---------
	- observations: the observations on which to train the HMM
	  type: std::vector
	RETURNS
	-------
	- a tuple that gives the backward probability and the backward algorithm matrix
	  type: std::pair */

	if ( not finalised ) throw NotFinalised();

	std::vector< std::vector< double > > backward( states.size(), std::vector< double >( observations.size() + 1, NAN ) );
	std::pair< std::string, std::string > key;
	double logProbability, pdfEval;
	unsigned int k;

	/*-----------INITIALISATION----------- */
	/*we must start in the end state */
	backward[ endIdx ][ observations.size() ] = 0.0;

	/*transitions between silent states between the last observation and the end */
	for ( int i = states.size() - 1; i >= 0; i-- ){

		if (i != endIdx ){
			logProbability = NAN;

			/*for all silent states... */
			for ( unsigned int j = 0; j < states[ i ].outTransitions.size(); j++ ){
			
				k = states[ i ].outTransitions[ j ];

				if ( k < silentStart or k >= i  ) continue;

				key = make_pair( states[ i ].name, states[ k ].name );

				logProbability = lnSum( logProbability, \
							lnProd( backward[ k ][ observations.size() ], \
								eln( transitions.at( key ) ) ) );
			}
			backward[ i ][ observations.size() ] = logProbability;
		}
	}

	/*-----------RECURSION----------- */
	for ( int t = observations.size() - 1; t >= 0; t-- ){
		
		/*silent states */
		for ( unsigned int i = states.size() - 1; i >= silentStart; i-- ){

			logProbability = NAN;

			/*for all emitting states... */
			for ( unsigned int j = 0; j < states[ i ].outTransitions.size(); j++ ){

				k = states[ i ].outTransitions[ j ];

				if ( k >= silentStart) continue; //only for emitting states

				key = make_pair( states[ i ].name, states[ k ].name );

				pdfEval = eln( states[ k ].dist -> pdf( observations[ t ] ) );

				logProbability = lnSum( logProbability, \
							lnProd( backward[ k ][ t + 1 ], \
								lnProd( eln( transitions.at( key ) ), \
									pdfEval ) ) );
			}
			/*for all silent states... */
			for ( unsigned int j = 0; j < states[ i ].outTransitions.size(); j++ ){

				k = states[ i ].outTransitions[ j ];

				if ( k < silentStart ) continue; //only for silent states
			
				key = make_pair( states[ i ].name, states[ k ].name );

				logProbability = lnSum( logProbability, \
                                            		lnProd( backward[ k ][ t ], \
								eln( transitions.at( key ) ) ) );

			}
			backward[ i ][ t ] = logProbability;
		}

		/*emitting states */
		for ( unsigned int i = 0; i < silentStart; i++ ){ 
			
			logProbability = NAN;
		
			/*for all emitting states... */
			for ( unsigned int j = 0; j < states[ i ].outTransitions.size(); j++ ){

				k = states[ i ].outTransitions[ j ];

				if ( k >= silentStart ) continue; //do only for emitting states

				key = make_pair( states[ i ].name, states[ k ].name );
				
				pdfEval = eln( states[ k ].dist -> pdf( observations[ t ] ) );

				logProbability = lnSum( logProbability, \
							lnProd( backward[ k ][ t + 1 ], \
								lnProd( eln( transitions.at( key ) ), \
									pdfEval ) ) );
			}
			/*for all silent states... */
			for ( unsigned int j = 0; j < states[ i ].outTransitions.size(); j++ ){

				k = states[ i ].outTransitions[ j ];

				if ( k < silentStart ) continue; //do only for silent states

				key = make_pair( states[ i ].name, states[ k ].name );

				logProbability = lnSum( logProbability, \
							lnProd( backward[ k ][ t ], \
								eln( transitions.at( key ) ) ) );
			}
			backward[ i ][ t ] = logProbability;
		}
	}

	/*-----------TERMINATION----------- */
	double backwardProb = NAN;
	for ( unsigned int j = 0; j < states[ startIdx ].outTransitions.size(); j++ ){

		k = states[ startIdx ].outTransitions[ j ];

		key = make_pair( STARTNAME, states[ k ].name );	

		backwardProb = lnSum( backwardProb, \
					lnProd( backward[ k ][ 0 ], \
						eln( transitions.at( key ) ) ) );
	}

	#if TEST_BACKWARD
	for (int i = 0; i < backward.size(); i++){
		for (int j = 0; j < backward[i].size(); j++){
			std::cout << backward[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << backwardProb << std::endl;
	#endif

	std::pair< double, std::vector< std::vector< double > > > backwardPair = std::make_pair( backwardProb, backward );

	return backwardPair;
}


//-------------------------------------------------------------------------------------------------------------------------------
std::pair< double, std::vector< std::string > > HiddenMarkovModel::viterbi( std::vector< double > &observations ){

	if ( not finalised ) throw NotFinalised();

	std::vector< std::vector< double > > viterbiLattice( states.size(), std::vector< double >( observations.size() + 1, NAN ) );
	std::vector< std::vector< unsigned int > > backtraceS( states.size(), std::vector< unsigned int >( observations.size() + 1 ) ); /*stores state indices for the Viterbi backtrace */
	std::vector< std::vector< unsigned int > > backtraceT( states.size(), std::vector< unsigned int >( observations.size() + 1 ) ); /*stores observation indices for the Viterbi backtrace */

	std::vector< std::string > viterbiSeq; /*final vector of the state names, in order, that most likely produced the observation */

	std::pair< std::string, std::string > key;
	double logProbability, pdfEval;
	unsigned int colIdx;
	unsigned int k;

	/*-----------INITIALISATION----------- */
	viterbiLattice[ startIdx ][ 0 ] = 0.0; //log probability of the start state is ln(1) = 0

	/*transitions between silent states before we emit the first symbol */
	for ( unsigned int i = silentStart; i < states.size(); i++ ){

		if ( i != startIdx ){

			/*for all silent states... */
			for ( unsigned int j = 0; j < states[ i ].inTransitions.size(); j++ ){

				k = states[ i ].inTransitions[ j ];
				if (k < silentStart ) continue; //do only for silent states

				key = make_pair( states[ k ].name, states[ i ].name );

				logProbability = lnProd( viterbiLattice[ k ][ 0 ], \
							 eln( transitions.at( key ) ) );

				/*if logProbability is larger than the current value of viterbiLattice[i][0], then make this the new value and update the viterbi backtrace */
				if ( lnGreaterThan( logProbability, viterbiLattice[ i ][ 0 ] ) ) {
					viterbiLattice[ i ][ 0 ] = logProbability;
					backtraceS[ i ][ 0 ] = k;
					backtraceT[ i ][ 0 ] = 0;
				}
			}
		}
	}

	/*-----------RECURSION----------- */
	for ( unsigned int t = 0; t < observations.size(); t++ ){

		colIdx = t + 1;

		/*emitting states */
		for ( unsigned int i = 0; i < silentStart; i++ ){ 
			
			pdfEval = eln( states[ i ].dist -> pdf( observations[ t ] ) );

			for ( unsigned int j = 0; j < states[ i ].inTransitions.size(); j++ ) {

				k = states[ i ].inTransitions[ j ];

				key = make_pair( states[ k ].name, states[ i ].name );
				

				logProbability = lnProd( viterbiLattice[ k ][ colIdx - 1 ], \
							 lnProd( eln( transitions.at( key ) ), \
								 pdfEval ) );

				/*if logProbability is larger than the current value of viterbiLattice[i][0], then make this the new value and update the viterbi backtrace */
				if ( lnGreaterThan( logProbability, viterbiLattice[ i ][ colIdx ] ) ) {
					viterbiLattice[ i ][ colIdx ] = logProbability;
					backtraceS[ i ][ colIdx ] = k;
					backtraceT[ i ][ colIdx ] = t;
				}
			}
		}
		/*silent states */
		for ( unsigned int i = silentStart; i < states.size(); i++ ){

			for ( unsigned int j = 0; j < states[ i ].inTransitions.size(); j++ ){

				k = states[ i ].inTransitions[ j ];

				key = make_pair( states[ k ].name, states[ i ].name );

				logProbability = lnProd( viterbiLattice[ k ][ colIdx ], \
                                        		 eln( transitions.at( key ) ) );

				/*if logProbability is larger than the current value of viterbiLattice[i][0], then make this the new value and update the viterbi backtrace */
				if ( lnGreaterThan( logProbability, viterbiLattice[ i ][ colIdx ] ) ) {
					viterbiLattice[ i ][ colIdx ] = logProbability;
					backtraceS[ i ][ colIdx ] = k;
					backtraceT[ i ][ colIdx ] = colIdx;
				}
			}
		}
	}

	/*-----------TERMINATION----------- */
	/*score the viterbi sequence */
	double viterbiScore = viterbiLattice[ endIdx ][ colIdx ];

	/*if this path is impossible, don't traceback and terminate */
	if ( std::isnan( viterbiScore ) ) return std::make_pair( viterbiScore, viterbiSeq );

	/*-----------TRACEBACK----------- */
	/*start at the end state */
	viterbiSeq.push_back( states[ endIdx ].name );
	unsigned int tracebackStates_old = endIdx;
	unsigned int tracebackStates_new;
	unsigned int tracebackT = colIdx;

	/*iterate through the traceback until we reach the start state */
	while (tracebackStates_old != startIdx){

		tracebackStates_new = backtraceS[ tracebackStates_old ][ tracebackT ];
		tracebackT = backtraceT[ tracebackStates_old ][ tracebackT ];

		viterbiSeq.push_back( states[ tracebackStates_new ].name );
		tracebackStates_old = tracebackStates_new;

	}

	/*we added states in the reverse order, so reverse the list to get the right way round*/
	std::reverse( viterbiSeq.begin(), viterbiSeq.end() );

	return std::make_pair( viterbiScore, viterbiSeq );
}


//-------------------------------------------------------------------------------------------------------------------------------
void HiddenMarkovModel::BaumWelch( std::vector< std::vector< double > > &observations, double tolerance, int maxIter, double emissionInertia, bool trainTransitions, int threads, bool verbose ){

	if ( not finalised ) throw NotFinalised();

	if ( verbose ){
		std::cout << "Starting Baum-Welch iterations..." << std::endl;
		std::cout << "\tTraining on " << observations.size() << " observations..." << std::endl;
	}

	std::pair< std::string, std::string > key;
	int tetheredStateIdx;
	double muRunning, sigmaRunning, normaliserRunning;
	double improvement = std::numeric_limits< double >::max();
	int iter = 0;
	double logProbabilityOld;
	
	while ( improvement > tolerance && iter < maxIter ){

		/*vectors for sums over observations */
		std::vector< double > muNumerator( states.size(), NAN ), sigmaNumerator( states.size(), NAN ), forwardSum( states.size(), NAN ), totalTransCount( states.size(), NAN );

		std::vector< std::vector< double > > softTrans( states.size(), std::vector< double >( states.size(), NAN ) );

		/*Create a boolean vector of the length of the emitting states.  We'll use this to keep track of which states we've visited while computing new emissions for tethered states */
		double logProbabilityNew = 0.0;
		bool abort = false;

		/*parallel loop over observations */
		#pragma omp parallel for default(none) shared(muNumerator,sigmaNumerator,forwardSum,observations,trainTransitions,softTrans,abort) private(key) reduction(+:logProbabilityNew) num_threads(threads)
		for ( unsigned int obs = 0; obs < observations.size(); obs++ ){
			try{
			/*private variables for openmp */
			double runningProbability;

			/*calculate forward and backward probabilities and matrices for this observation */
			std::pair< double, std::vector< std::vector< double > > > forwardPair = forward( observations[ obs ] );
			double sequenceProb = forwardPair.first;
			if ( std::isnan( sequenceProb ) ) continue;
			logProbabilityNew += sequenceProb;			
			std::vector< std::vector< double > > alpha = forwardPair.second;
			std::pair< double, std::vector< std::vector< double > > > backwardPair = backward( observations[ obs ] );
			std::vector< std::vector< double > > beta = backwardPair.second;

			/*-----------TRAIN EMISSIONS----------- */
			for ( unsigned int i = 0; i < silentStart; i++ ){
				/*only train normal distributions */
				if ( states[ i ].dist -> identify() == "NormalDistribution" ){
					
					/*for each state, update the mu, sigma values for each normal emitting state */
					for ( unsigned int t = 0; t < observations[ obs ].size(); t++){

						muNumerator[ i ] = lnSum( muNumerator[ i ], \
						                          lnProd( eln( observations[ obs ][ t ] ), \
						                                  lnProd( alpha[ i ][ t + 1 ], \
						                                          lnQuot( beta[ i ][ t + 1 ], \
						                                                  sequenceProb ) ) ) );
						sigmaNumerator[ i ] = lnSum( sigmaNumerator[ i ], \
						                             lnProd( eln( pow( observations[ obs ][ t ] - states[ i ].dist -> trained_param1, 2.0 ) ), \
						                                     lnProd( alpha[ i ][ t + 1 ], \
						                                            lnQuot( beta[ i ][ t + 1 ], \
						                                                    sequenceProb ) ) ) );
						forwardSum[ i ] = lnSum( forwardSum[ i ], \
						                         lnProd( alpha[ i ][ t + 1 ], \
						                                 lnQuot( beta[ i ][ t + 1 ], \
						                                         sequenceProb ) ) );
					}
				}
			}	

			/*-----------TRAIN TRANSITIONS----------- */
			if ( trainTransitions == true ){
				/*loop over ALL STATES */
				//#pragma omp critical
				for ( unsigned int i = 0; i < states.size(); i++ ){

					/*loop over EMITTING STATES */
					/*update the transitions by looping over all states again, so we get all {states x states} combinations */
					for ( unsigned int j = 0; j < silentStart; j++ ){

						key = make_pair( states[ i ].name, states[ j ].name );
						
						double currentTransProb;
						if ( transitions.count( key ) < 1 ) currentTransProb = 0.0;
						else currentTransProb = transitions.at( key );

						runningProbability = NAN;

						for ( unsigned int t = 0; t < observations[ obs ].size(); t++ ){

							runningProbability = lnSum( runningProbability, \
								                        lnProd( alpha[ i ][ t ], \
								                                lnProd( eln( currentTransProb ), \
								                                        lnProd( beta[ j ][ t + 1 ], \
								                                                eln( states[ j ].dist -> pdf( observations[ obs ][ t ] ) ) ) ) ) );
						}
						/*divide softTrans by the sequence probability */
						softTrans[ i ][ j ] = lnSum( softTrans[ i ][ j ], \
                                                			     lnQuot( runningProbability, \
                                                				     sequenceProb ) );
					}

					/*loop over SILENT STATES */
					for ( unsigned int j = silentStart; j < states.size(); j++ ){

						key = make_pair( states[ i ].name, states[ j ].name );

						double currentTransProb;
						if ( transitions.count( key ) < 1 ) currentTransProb = 0.0;
						else currentTransProb = transitions.at( key );

						runningProbability = NAN;

						for ( unsigned int t = 0; t < observations[ obs ].size() + 1; t++ ){

							runningProbability = lnSum( runningProbability, \
						                                    lnProd( alpha[ i ][ t ], \
						                                	    lnProd( eln( currentTransProb ), \
						                                		    beta[ j ][ t ] ) ) );
						}
						/*divide softTrans by the sequence probability */
						softTrans[ i ][ j ] = lnSum( softTrans[ i ][ j ], \
                                                			     lnQuot( runningProbability, \
                                                				     sequenceProb ) );
					}
				}
			}
		}
		catch ( NegativeLog &nl ){
			abort = true;
			obs = observations.size();
		}
		}

		if (abort) throw NumericalInstability();

		/*-----------UPDATE EMISSIONS----------- */
		/*update the gaussian distribution parameters */
        
		/*reset the visited vector to zero from the last round of training */
		std::vector< bool > visited( silentStart, false );
        
		/*for each emitting state... */
		for ( unsigned int i = 0; i < silentStart; i++ ){
            
			/*Check if it's a normal distribution and that we haven't been here before. */
			if ( states[ i ].dist -> identify() == "NormalDistribution" && visited[ i ] == false ){
		        
				/*Reset the running totals to zero from the last round of training */
				muRunning = NAN;
				sigmaRunning = NAN;
				normaliserRunning = NAN;
				    
				/*For each state i, if it's untethered, then tethers[ i ] will just return that state's index - a vector of length 1.  Otherwise, it will return a vector of length > 1 with all of the indices of the states that are tied to state i.  Iterate through these and keep a running total of the mu and sigma numerators, as well as the normalising factor for them */
				for ( unsigned int j = 0; j < tethers[ states[ i ].dist ].size(); j++ ){

					tetheredStateIdx = tethers[ states[ i ].dist ][ j ];

					muRunning = lnSum( muRunning, \
							   muNumerator[ tetheredStateIdx ] );
					sigmaRunning = lnSum( sigmaRunning, \
							      sigmaNumerator[ tetheredStateIdx ] );
					normaliserRunning = lnSum( normaliserRunning, \
								   forwardSum[ tetheredStateIdx ] );
				        
					/*note that we've visited this state so we don't retrain */
					visited[ tetheredStateIdx ] =  true;
				}
		
				states[ i ].dist -> trained_param1 = emissionInertia*(states[ i ].dist -> trained_param1) + (1.0-emissionInertia)*(eexp( lnQuot( muRunning, normaliserRunning ) ));

				/*we've computed the variance, so take sqrt to get the standard deviation */
				states[ i ].dist -> trained_param2 = emissionInertia*(states[ i ].dist -> trained_param2) + (1.0-emissionInertia)*(sqrt( eexp( lnQuot( sigmaRunning, normaliserRunning ) ) ));
			}
		}

		/*-----------UPDATE TRANSITIONS----------- */
		if ( trainTransitions == true ){

			/*calculate normalising factor totalTransCount for each state */
			for ( unsigned int i = 0; i < states.size(); i++ ){

				for ( unsigned int j = 0; j < states.size(); j++ ){
					
					totalTransCount[ i ] = lnSum( totalTransCount[ i ], 
								      softTrans[ i ][ j ] );
				}
			}

			/*normalise and update the transitions */
			for ( unsigned int i = 0; i < states.size(); i++ ){

				/*OUT: normalise the transition probabilities from state i */
				if ( lnGreaterThan( totalTransCount[ i ], NAN ) ){

					for ( unsigned int j = 0; j < states.size(); j++ ){
				
						/*update soft transition matrix */
						key = make_pair( states[ i ].name, states[ j ].name );

						transitions[ key ] = eexp( lnQuot( softTrans[ i ][ j ], \
										   totalTransCount[ i ] ) );
					}
				}

				/*IN: normalise the transition probabilities going into state i */
				for ( unsigned int j = 0; j < states.size(); j++ ){

					if ( lnGreaterThan( totalTransCount[ j ], NAN ) ){
				
						/*update soft transition matrix */
						key = make_pair( states[ j ].name, states[ i ].name );

						transitions[ key ] = eexp( lnQuot( softTrans[ j ][ i ], \
										   totalTransCount[ j ] ) );
					}
				}
			}
		}

		#if TEST_BAUMWELCH
		/* check probabilities sum to 1 */
		double total;
		double tol = 1e-8;

		for ( unsigned int i = 0; i < states.size(); i++ ){

			total = 0.0;

			for ( unsigned int j = 0; j < states.size(); j++ ){

				key = make_pair( states[ i ].name, states[ j ].name );

				if ( transitions.count( key ) > 0 ){

					total += transitions[ key ];

				}
			}

			if ( i != endIdx ){

				assert( std::abs( 1 - total ) < tol );

			}

		}
		#endif
        
		/*How did we do?  Calculate the training improvement by taking the difference between the total (over all observations) log probability of the previous training and the log 			probability of the current training.  Update the improvement so that the main while loop breaks if we reach the tolerance. */
		if ( iter > 0 ){
			improvement = logProbabilityNew - logProbabilityOld;
			if ( verbose ){
				std::cout << "\tImprovement: " << improvement << std::endl;
			}	
		}
		logProbabilityOld = logProbabilityNew;
		iter++;
	}
    
	/*Display a final message stating whether we managed to converge or not. */
	if ( improvement < tolerance ){
		if ( verbose ){
			std::cout << "\tConverged to tolerance " << tolerance << std::endl;
			std::cout << "Done." << std::endl;
		}
	}
	else{
		std::cout << "\tWARNING: congergence failed." << std::endl;
	}
}
