#include <iostream>

/* local header files */
#include "../hmm.h"


int main(){

	/* TEST FOR SORTING SILENT STATES */

	/*make a hmm object and a states vector to put the states in */
	HiddenMarkovModel hmm();
	std::vector< State > states;

	/*make a few distributions - one of each type */
	SilentDistribution sd(0,0);
	SilentDistribution *sp = &sd;
	NormalDistribution nd(0.0,1.0);
	NormalDistribution *np = &nd;
	UniformDistribution ud(0.0,1.0);
	UniformDistribution *up = &ud;

	/*push the state objects into the states vector, then add them to the model */
	states.push_back(State(np,"first_added",1.0));
	states.push_back(State(sp,"second_added",1.0));
	states.push_back(State(up,"third_added",1.0));
	states.push_back(State(sp,"fourth_added",1.0));
	hmm.add_state(states[0]);
	hmm.add_state(states[1]);
	hmm.add_state(states[2]);
	hmm.add_state(states[3]);

	/*add a few test transitions */
	hmm.add_transition(states[1], states[0], 1.0);
	hmm.add_transition(states[2], states[1], 1.0);
	hmm.add_transition(states[3], states[2], 1.0);

	/*print out the names of the states.  They will be in the order that we added them to the model */
	std::cout << "Printing the non-silent states:" << std::endl;
	for (int i = 0; i < hmm.states.size(); i++){
		std::cout << hmm.states[i].name << std::endl;
	}

	std::cout << "Finalise the HMM (quicksort states):" << std::endl;
	/*now finalise the model and test whether the quicksort on silent states worked */
	hmm.finalise();
	for (int i = 0; i < hmm.states.size(); i++){
		std::cout << hmm.states[i].name << std::endl;
	}


	/* TEST FOR PDF POLYMORPHISM 
	NormalDistribution nd(0.0,1.0);
	NormalDistribution *np = &nd;
	State state(np,"name",1.0);
	std::cout << state.dist -> pdf(0.0) << std::endl; 
	*/

	
	return 0;
}

