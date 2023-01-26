# caHMM
A light-weight, self-contained C++ library for hidden Markov models.  Features include:
- silent states,
- linked emitting states where the distributions will train together,
- explicitly defined logarithm computations for numerical stability,
- no dependencies (just STL).

The algorithms that the library provides are:
- forward,
- backward,
- Viterbi,
- Baum-Welch.

### Installation and Testing
After cloning the repository, nagivate to the caHMM directory and build the library by running:
```shell
make
```
Test that the library is running properly by running one of the tests:
```shell
g++ -L. -std=c++11 -o tests/simpleHMM tests/simpleHMM.cpp -lcaHMM -fopenmp
./tests/simpleHMM
```

### Basic Usage
The best way to get to grips with basic usage is to look through the simpleHMM.cpp example file in the examples directory.  This file captures all of the basic usage and features of the library.
Object constructors:
```c++
NormalDistribution nd( double mean, double standardDeviation );
UniformDistribution ud( double a, double b );
SilentDistribution sd( 0, 0 );
State st( Distribution *distPointer, std::string stateName, std::string stateTether, double weight );
HiddenMarkovModel hmm();
```


We can do the following with a HiddenMarkovModel object:
```c++
HiddenMarkovModel hmm();
hmm.add_state( State myState ); //adds a state object to the HMM object
hmm.add_transition( State myFirstState, State mySecondState, double p ); //state myFirstState can transition to mySecondState with probability p
hmm.finalise(); //must be called before any HMM algorithms (forward, backward, Viterbi, Baum-Welch) are called


std::pair< double, std::vector< std::vector< double > > > forward = hmm.forward( std::vector <double> observations );
forward.first; //the forward probability of observations
forward.second; //vector of vectors that gives the forward lattice of observations

std::pair< double, std::vector< std::vector< double > > > backward = hmm.backward( std::vector <double> observations );
backward.first; //the backward probability of observations
backward.second; //vector of vectors that gives the backward lattice of observations

std::pair< double, std::vector< std::string > > viterbi = hmm.viterbi( std::vector< double > observations );
viterbi.first; //the Viterbi probability of observations
viterbi.second; //the Viterbi path of observations

hmm.BaumWelch( std::vector< std::vector< double > > observations, double tolerance, int maxItererationsAllowed, bool trainTransitions, int threads );

std::stringstream ss = hmm.summarise(); //returns a human readable outcome of the HMM training
```

### Parallel Processing
For the Baum-Welch algorithm, parallel processing is done with OpenMP.

### Credit and Thanks
Inspiration was drawn from Jared Simpson's Nanopolish (https://github.com/jts/nanopolish) and Jacob Schreiber's Pomegranate (https://github.com/jmschrei/pomegranate).  Some of the functions for numerical stability were inspired by Tobias Mann's article (http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf).  All due credit and thanks to these developers.
