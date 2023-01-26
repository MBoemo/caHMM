from pomegranate import *
import numpy
numpy.set_printoptions(threshold=sys.maxsize)


#pomegranate==0.4.0
#numpy==1.8.0


def import_poreModel(filename):
#	takes the filename of an ONT pore model file and returns a map from kmer (string) to [mean,std] (list of floats)
#	ARGUMENTS
#       ---------
#	- filename: path to an ONT model file
#	  type: string
#	OUTPUTS
#       -------
#	- kmer2MeanStd: a map, keyed by a kmer, that returns the model mean and standard deviation signal for that kmer
#	  type: dictionary

	f = open(filename,'r')
	g = f.readlines()
	f.close()

	kmer2MeanStd = {}
	for line in g:
		if line[0] != '#' and line[0:4] != 'kmer': #ignore the header
			splitLine = line.split('\t')
			kmer2MeanStd[ splitLine[0] ] = [ float(splitLine[1]), float(splitLine[2]) ]
	g = None

	return kmer2MeanStd


def build_TrainingHMM(refSequence, thymidineModel, scalings):

	hmm = HiddenMarkovModel()

	windowSize = 10
	brduStart = windowSize - 5
	brduEnd = windowSize + refSequence[windowSize:windowSize+6].rfind('T')

	refLength = len(refSequence)

	#new HMM transition parameters
	internalM12I = 0.3475
	internalI2I = 0.5
	internalM12M1 = 0.4

	externalD2D = 0.3
	externalD2M1 = 0.7
	externalI2M1 = 0.5
	externalM12D = 0.0025
	externalM12M1 = 0.25

	###########################################################################################################################
	# Add States to Model
	#Create the HMM states.  Iterate through the reference sequence, and make the repeating HMM module for each position in the sequence.

	#two-dimensional array for the states, where the columns are the positions in the reference
	states = [[0 for x in range(refLength)] for y in range(3)]

	#the first base
	i = 0
	emissions = thymidineModel[refSequence[0:6]]
	level_mu = scalings[0] + scalings[1] * emissions[0]
	level_sig = scalings[2] * emissions[1]

	#print(refSequence[i:i+6])
	#print(level_mu,level_sig)

	states[0][i] = State( UniformDistribution(0, 250, frozen=True), name='Insertion_'+str(i) )    
	states[1][i] = State( NormalDistribution(level_mu, level_sig), name='Match_'+str(i) )    
	states[2][i] = State( None, name='Deletion_'+str(i) )  
	for j in range(3):
		hmm.add_state(states[j][i])

	#make an insertion state before the first base
	firstI = State( UniformDistribution(0, 250, frozen=True), name='Insertion_Pre' )
	hmm.add_state(firstI)

	#insertion state before the first base
	hmm.add_transition(hmm.start, firstI, 0.25) #start to the first insertion
	hmm.add_transition(firstI, firstI, 0.25) #self loop

	#to the base 1 insertion
	hmm.add_transition(states[0][i], states[0][i], internalI2I , group='internal_I-to-I')
	hmm.add_transition(states[1][i], states[0][i], internalM12I , group='internal_M-to-I')

	#to the base 1 match
	hmm.add_transition(firstI, states[1][i], 0.5) #first insertion to first match
	hmm.add_transition(states[1][i], states[1][i], internalM12M1 , group='internal_M-to-M')
	hmm.add_transition(hmm.start, states[1][0], 0.5)  #start to M
	
	#to the base 1 deletion
	hmm.add_transition(firstI, states[2][i], 0.25) #first insertion to the first deletion
	hmm.add_transition(hmm.start, states[2][0], 0.25) #start to D

	#the rest of the sequence
	for i, char in enumerate(refSequence[1:-5]):

		i += 1

		emissions = thymidineModel[refSequence[i:i+6]]

		#correct for shift/scale/var
		level_mu = scalings[0] + scalings[1] * emissions[0]
		level_sig = scalings[2] * emissions[1]

		#print(refSequence[i:i+6])
		#print(level_mu,level_sig)
		
		#create states for this nucleotide
		states[0][i] = State( UniformDistribution(0, 250, frozen=True), name='Insertion_'+str(i) )    
		states[1][i] = State( NormalDistribution(level_mu, level_sig), name='Match_'+str(i) )    
		states[2][i] = State( None, name='Deletion_'+str(i) )    

		#add the state objects to the hmm model object
		for j in range(3):
			hmm.add_state(states[j][i])

		#internal transitions for this nucleotide
		hmm.add_transition(states[1][i], states[1][i], internalM12M1 , group='internal_M-to-M')
		hmm.add_transition(states[1][i], states[0][i], internalM12I , group='internal_M-to-I')
		hmm.add_transition(states[0][i], states[0][i], internalI2I , group='internal_I-to-I')

	#this is really just a safety thing... get the last index in the iterator
	for last_i,x in enumerate(refSequence[:-5]): 
		pass	

	###########################################################################################################################
	# Handle Transitions Between Modules, Handle Analogue Branch

	#We have to reach forward to the next state (state i+1) so it's easier to just do this in a separate loop from the internal transitions one
	for i, char in enumerate(refSequence[:-5]): 

		#Don't execute this if we're at the end, because there's no i+1 to reach forward to.
		if i != last_i:

			hmm.add_transition(states[1][i], states[2][i+1], externalM12D, group='external_M-to-D')
			hmm.add_transition(states[2][i], states[2][i+1], externalD2D, group='external_D-to-D')
			hmm.add_transition(states[2][i], states[1][i+1], externalD2M1, group='external_D-to-M')
			hmm.add_transition(states[1][i], states[1][i+1], externalM12M1, group='external_M-to-M')
			hmm.add_transition(states[0][i], states[1][i+1], externalI2M1, group='external_I-to-M')

	###########################################################################################################################
	# Handle Start and End

	#handle end states
	hmm.add_transition(states[0][last_i], hmm.end, externalI2M1 )
	hmm.add_transition(states[1][last_i], hmm.end, externalM12M1 + externalM12D)
	hmm.add_transition(states[2][last_i], hmm.end, 1.0)

	#bake the model
	#hmm.bake(merge='all',verbose=True)
	hmm.bake()
	return hmm

###########################################################################################################################
# MAIN


#make some simple test data - uncomment one or the other and check against simpleHMM.cpp

sequnceSnippet = 'CTCCACTCGTTACCCTGTCCCATTCA'
events = [98.5804, 94.1981, 85.6842, 87.9788, 91.598, 109.494, 86.048, 72.0098, 84.7888, 97.4931, 95.5623, 100.11, 91.7845, 86.4678, 92.8339, 88.0668, 96.8873, 95.0026, 94.0932, 97.2692, 92.5841, 101.299, 83.0631, 98.8736, 109.134, 105.496, 95.7022, 97.3252, 94.303, 73.0359, 70.9371, 86.8595, 89.7138, 82.4568, 87.727, 93.1837, 87.3073]
scalings = [4.55653, 0.956433, 1.75943]

'''
sequnceSnippet = 'CACACACACATCCTAACACTACCCTA'
events = [96.6349, 95.5902, 94.0932, 102.325, 99.7597, 102.418, 95.3664, 92.0644, 94.1911, 95.2824, 97.1946, 101.299, 104.563, 103.607, 100.627, 101.858, 97.8476, 103.178, 95.8621, 102.884, 96.6816, 99.0602, 95.8701, 100.646, 97.5677, 82.3502, 79.0056, 81.4868, 116.06, 107.175, 96.7282, 96.8215, 112.492, 110.487, 114.031, 105.776, 91.7845, 88.893, 93.9113, 92.0644, 100.179, 92.3442, 92.624, 95.2824, 96.3551, 101.802, 112.492, 108.994, 93.8633, 93.6034, 103.864, 99.7877, 97.5677, 100.515, 102.511]
scalings = [4.55653, 0.956433, 1.75943]
'''

poreModel = import_poreModel('template_median68pA.model')
hmm = build_TrainingHMM(sequnceSnippet, poreModel, scalings)

print('Log Probability:',hmm.log_probability(events))


