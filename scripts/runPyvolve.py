import sys
import os
import math
import numpy as np
import os.path
from os import path
import argparse
from ete3 import Tree
import time

#script to run pyvolve

parser = argparse.ArgumentParser(description='Compare simulators, in particular phastSim, seq-gen, indelible and pyvolve.')
parser.add_argument('--path',default="", help='Path where to run simulations.')
parser.add_argument('--reference',default="MN908947.3.fasta", help='File containing the reference genome to be used as root genome. To be found in the folder specified with --path. By default the reference is generated randomly.')
parser.add_argument('--treeFile',default="exampleTree.tree", help='Name of file containing the tree used to simulate sequences (assumed within the --path and in newick format).')
parser.add_argument('--scale',default=1.0,type=float, help='Scale the simulation tree by this amount (default 1.0). Branch lengths are assumed in terms of expected substitutions per site (more or less, as frequencies changes through time, total mutation rate might also  change).')
parser.add_argument("--seed", help="Seed for random simulator in genes simulations", type=int, default=1)
parser.add_argument("--mutationRates", help="Mutation rates, by default using the neutral rates estimated from SARS-CoV-2; so far only exactly 12 input values allowed (r_AC, r_AG, r_AT, r_CA, etc...) corresponding to an UNREST nucleotide substitution model.", type=float, nargs='+',  default=[0.039, 0.310, 0.123, 0.140, 0.022, 3.028, 0.747, 0.113, 2.953, 0.056, 0.261, 0.036])
parser.add_argument("--categoryProbs", help="Probabilities of rate categories. They are supposed to sum up to one, but if they don't they are normalized to do so. By default only 1 category is simulated.", type=float, nargs='+',  default=[1.0])
parser.add_argument("--categoryRates", help="Rates of site categories. The overall mutation rate is renormalized so that the expected number of bstitutions per branch length unit per site is 1. By default only 1 category of rate 1.0 is simulated. The number of rates has to be the same as the number of categgory probabilities, otherwise an error is thrown.", type=float, nargs='+',  default=[1.0])
parser.add_argument('--outputFile',default="sars-cov-2_simulation_pyvolve_output", help='Output file name containing the simulated genomes in succint format. The file will be created within the folder specified with --path.')
args = parser.parse_args()

pathSimu=args.path
reference=args.reference
np.random.seed(args.seed)
seed=args.seed
scale=args.scale
mutationRates=args.mutationRates
categoryProbs=args.categoryProbs
categoryRates=args.categoryRates
treeFile=args.treeFile
outputFile=args.outputFile

#collect reference
ref=""
if reference!="":
    file=open(reference)
    line=file.readline()
    while line!="":
        line=file.readline()
        ref+=line.replace("\n","")
    file.close()
    print("\n Finished reading reference genome at "+reference+" with "+str(len(ref))+" bases.")
    refList=list(ref)
else:
    print("random ancestor not implemente yet!")
    exit()

#SARS-CoV-2 genome annotation - not used yet but will be useful when simulating under a codon model.
geneEnds=[[266,13468],[13468,21555],[21563,25384],[25393,26220],[26245,26472],[26523,27191],[27202,27387],[27394,27759],[27894,28259],[28274,29533],[29558,29674]]

#define the mutation matrix
if len(mutationRates)==12:
        #possible alleles
        alleles={"A":0,"C":1,"G":2,"T":3}
        allelesList=["A","C","G","T"]
        nAlleles=4

        #print("\n Assuming UNREST nucleotide substitution model.")
        mutMatrix=np.zeros((4,4),dtype=float)
        index=0
        for i in range(4):
            sum=0.0
            for j in range(4):
                if j!=i:
                    mutMatrix[i][j]=mutationRates[index]
                    sum+=mutationRates[index]
                    index+=1
            mutMatrix[i][i]=-sum
else:
    print("I still have to include this substitution model in pyvolve simulations")
    exit()

nCat=len(categoryProbs)

if nCat!=len(categoryRates):
    print("Issue with number of category probs "+str(len(categoryProbs))+" and number of category rates "+str(len(categoryRates)))
    exit()

    
sum=0.0
for i in categoryProbs:
    sum+=i
for i in range(nCat):
    categoryProbs[i]=categoryProbs[i]/sum
if sum>1.000001 or sum<0.999999:
    print("\n Normalizing probabilities of site categories. New probabilities:")
    print(categoryProbs)


#run pyvolve

print("Starting pyvolve timer")
import pyvolve
start = time.time()
pyvolveTree = pyvolve.read_tree(file = pathSimu+treeFile, scale_tree = args.scale)
#pyvolveTree = pyvolve.read_tree(tree = tString2, scale_tree = args.scale)
nucModel=pyvolve.Model("custom", {"matrix":mutMatrix}, alpha = 0.5, num_categories = len(categoryRates))
partitions=pyvolve.Partition(models = nucModel, root_sequence = ref)
my_evolver = pyvolve.Evolver(tree = pyvolveTree, partitions = partitions)
my_evolver(seqfile=pathSimu+outputFile, algorithm=1) # Algorithm = 1 uses the Gillespie algorithm. 
time2 = time.time() - start
print("Pyvolve timer ended")
print(time2)

exit()
























