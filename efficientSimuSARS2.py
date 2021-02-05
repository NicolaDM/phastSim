# import sys
# import os
# import math
import numpy as np
# import os.path
# from os import path
# import argparse
from ete3 import Tree
import time
import phastSim
from importlib import reload
reload(phastSim)

"""
Script that simulates sequence evolution along a given input phylogeny.
The algorithm (based on Gillespie approach) is fast for trees with short branches,
as for example in genomic epidemiology.
It can be instead be slower than traditional approaches when longer branch lengths are considered.

example run:
python3 efficientSimuSARS2.py --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 --scale 2.0 --hyperMutProbs 0.01
 --hyperMutRates  100.0 --codon --treeFile rob-11-7-20.newick --createFasta

Possible future extensions: 
add wrapper to simulate whole sars-cov-2 genome evolution ?
CONSIDER TO EXTEND THE RANGE OF ALLOWED MODELS to allow easier specification og e.g. HKY, JC, etc?
ALLOW discretized gamma?
ALLOW INDELS ?
ALLOW TREE GENERATED ON THE GO, WITH MUTATIONS AFFECTING FITNESS and therefore birth rate of a lineage ?
allow mixtures of different models for different parts of the genome (e.g. coding and non-coding at the same time)
"""

# setup the argument parser and read the arguments from command line
parser = phastSim.setup_argument_parser()
args = parser.parse_args()

# instantiate a phastSim run. This class holds all arguments and constants, which can be easily called as e.g.
# sim_run.args.path or sim_run.const.alleles
sim_run = phastSim.phastSim_run(args=args)
self=sim_run

pathSimu=args.path
reference=args.reference
rootGenomeLength=args.rootGenomeLength
rootGenomeFrequencies=args.rootGenomeFrequencies
np.random.seed(args.seed)
seed=args.seed
scale=args.scale
alpha=args.alpha
invariable=args.invariable
mutationRates=args.mutationRates
categoryProbs=args.categoryProbs
categoryRates=args.categoryRates
hyperMutProbs=args.hyperMutProbs
hyperMutRates=args.hyperMutRates
verbose=args.verbose
outputFile=args.outputFile
createNewick=args.createNewick
createFasta=args.createFasta
createPhylip=args.createPhylip
treeFile=args.treeFile
hierarchy=not args.noHierarchy

codon=args.codon
omegaAlpha=args.omegaAlpha
omegaCategoryProbs=args.omegaCategoryProbs
omegaCategoryRates=args.omegaCategoryRates

alleles = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3, "u": 3, "U": 3}
allelesList = ["A", "C", "G", "T"]
nAlleles = 4


# initialise the root genome. Reads either from file or creates a genome in codon or nucleotide mode
ref, refList = sim_run.init_rootGenome()

# set up substitution rates
mutMatrix = sim_run.init_substitution_rates()

# set up gamma rates
gammaRates = sim_run.init_gamma_rates()

# set up hypermutation rates
hyperCategories = sim_run.init_hypermutation_rates()

# set up codon substitution model
if sim_run.args.codon:
	omegas = sim_run.init_codon_substitution_model()
	gammaRates, omegas = phastSim.check_start_stop_codons(ref=ref, gammaRates=gammaRates, omegas=omegas)
else:
	omegas = None


#SARS-CoV-2 genome annotation - not used yet but will be useful when simulating under a codon model.
#geneEnds=[[266,13468],[13468,21555],[21563,25384],[25393,26220],[26245,26472],[26523,27191],[27202,27387],[27394,27759],[27894,28259],[28274,29533],[29558,29674]]


# Loads a tree structure from a newick string in ETE2. The returned variable t is the root node for the tree.
start = time.time()
t = Tree(pathSimu+treeFile)
time1 = time.time() - start
print("Time for reading tree with ETE3: "+str(time1))
	


#save information about the categories of each site on a file
file=open(pathSimu+outputFile+".info","w")
file.write("pos\t"+"cat\t"+"hyperCat\t"+"hyperAlleleFrom\t"+"hyperAlleleTo\n")


#Hierarchical approach (DIVIDE ET IMPERA ALONG THE GENOME),
# WHEN THE RATES ARE UPDATED, UPDATE ONLY RATE OF THE SITE AND OF ALL THE NODES ON TOP OF THE HYRARCHY.
# THAT IS, DEFINE A tree STRUCTURE, WITH TERMINAL NODES BEING GENOME LOCI AND WITH INTERNAL NODES BING MERGING GROUPS OF LOCI CONTAINING INFORMATION ABOUT THEIR CUMULATIVE RATES.
# THIS WAY UPDATING A MUTATION EVENT REQUIRES COST LOGARITHMIC IN GENOME SIZE EVEN IF EVERY SITE HAS A DIFFERENT RATE. 
if hierarchy:
	#CODONS: maybe don't create all the matrices from the start (might have too large an memory and time preparation cost).
	#instead, initialize only the rates from the reference allele (only 9 rates are needed), and store them in a dictionary at level 0 terminal nodes, and when new codons at a position 
	#are reached, extend the dictionary and calculate these new rates. Most positions will have only a few codons explored.

	# instantiate a GenomeTree with all needed rates and categories
	genome_tree = phastSim.GenomeTree_hierarchical(
		nCodons=sim_run.nCodons,
		codon=sim_run.args.codon,
		ref=ref,
		gammaRates=gammaRates,
		omegas=omegas,
		mutMatrix=mutMatrix,
		hyperCategories=hyperCategories,
		hyperMutRates=sim_run.args.hyperMutRates,
		file=file,
		verbose=sim_run.args.verbose)

	# populate the GenomeTree
	genome_tree.populateGenomeTree(node=genome_tree.genomeRoot)

	# normalize all rates
	genome_tree.normalize_rates(scale=sim_run.args.scale)

	time2 = time.time() - start
	print("Total time after preparing for simulations: " + str(time2))

	# NOW DO THE ACTUAL SIMULATIONS. DEFINE TEMPORARY STRUCTURE ON TOP OF THE CONSTANT REFERENCE GENOME TREE.
	# define a multi-layered tree; we start the simulations with a genome tree.
	# as we move down the phylogenetic tree, new layers are added below the starting tree.
	# Nodes to layers below link to nodes above, or to nodes on the same layer, but never to nodes in the layer below.
	# while traversing the tree, as we move up gain from a node back to its parent
	# (so that we can move to siblings etc), the nodes in layers below the current one are simply "forgotten"
	# (in C they could be de-allocated, but the task here is left to python automation).
	genome_tree.mutateBranchETEhierarchy(t, genome_tree.genomeRoot, 1, sim_run.args.createNewick)




#use simpler approach that collates same rates along the genome - less efficient with more complex models.
else:

	# instantiate a genome tree for the non hierarchical case
	genome_tree = phastSim.GenomeTree_simpler(
		nCat=sim_run.nCat,
		ref=ref,
		mutMatrix=mutMatrix,
		categories=sim_run.categories,
		categoryRates=sim_run.args.categoryRates,
		hyperMutRates=sim_run.args.hyperMutRates,
		hyperCategories=hyperCategories,
		file=file,
		verbose=sim_run.args.verbose)

	# prepare the associated lists and mutation rate matrices
	genome_tree.prepare_genome()

	# normalize the mutation rates
	genome_tree.normalize_rates(scale=sim_run.args.scale)

	time2 = time.time() - start
	print("Total time after preparing for simulations: " + str(time2))

	# Run sequence evolution simulation along tree
	genome_tree.mutateBranchETE(t, genome_tree.muts, genome_tree.totAlleles, genome_tree.totMut, genome_tree.extras,
								sim_run.args.createNewick)




time2 = time.time() - start
print("Total time after simulating sequence evolution along tree with Gillespie approach: " + str(time2))

# depending on the type of genome_tree, this automatically uses the correct version
genome_tree.write_genome_short(tree=t, output_path=args.path, output_file=args.outputFile)

time3 = time.time() - start
print("Total time after writing short file: " + str(time3))


# If requested, create a newick output
if createNewick:
	file = open(args.path + args.outputFile + ".tree", "w")
	newickTree = phastSim.writeGenomeNewick(t) + ";\n"
	file.write(newickTree)
	file.close()
	
	time3 = time.time() - start
	print("Total time after writing newick file: "+str(time3))


# If requested, create a fasta output.
# Depending on which type of simulation (hierarchical or non-hierarchical), the correct function will be used
if createFasta:
	genome_tree.write_genome(tree=t, output_path=args.path, output_file=args.outputFile, refList=refList)
	
	time3 = time.time() - start
	print("Total time after writing fasta file: "+str(time3))
	


#If requested, create a phylip output
if createPhylip:
	file=open(pathSimu+outputFile+".phy","w")
	file.write("\t"+str(len(t))+"\t"+str(len(ref))+"\n")
	
	if hierarchy:
		#function to write a complete sequence output iteratively
		def writeGenomePhylip(node,file,nRefList):
			#update list
			for m in node.mutations:
				nRefList[m[0]]=allelesList[m[2]]
			#print leaf entry to file
			if node.is_leaf():
				file.write(node.name+"\t"+(''.join(nRefList))+"\n")
			else:
				for c in node.children:
					writeGenomePhylip(c,file,nRefList)
			#de-update the list so it can be used by siblings etc. 
			for n in range(len(node.mutations)):
				m=node.mutations[len(node.mutations)-(n+1)]
				nRefList[m[0]]=allelesList[m[1]]
		writeGenomePhylip(t,file,refList)
	else:
		#function to write a phylip output iteratively
		def writeGenomePhylip(node,file):
			if node.is_leaf():
				seq=genomeSeq(node.mutations)
				file.write(node.name+"\t"+seq+"\n")
			for c in node.children:
				writeGenomePhylip(c,file)
		writeGenomePhylip(t,file)
	file.close()
	
	time3 = time.time() - start
	print("Total time after writing phylip file: "+str(time3))

elapsedTime = time.time() - start
print("Overall time: "+str(elapsedTime))










exit()








