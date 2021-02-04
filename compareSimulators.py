import sys
import os
#import math
import numpy as np
import os.path
#from os import path
import argparse
from ete3 import Tree
import time

#NICOLA THINGS TO DO: 
#allow codon models

#ALLOW RANDOM ROOT GENOME in phastSim?
#allow discretised gamma distribution in phastSim like in the others?




#script to compare the running time of phastSim with pyvolve, seq-gen and indelible.

# python2 compareSimulators.py --path /Users/demaio/Desktop/coronavirus/simulations/ --nLeaves 100 --categoryProbs 0.2 0.2 0.2 0.2 0.2 --categoryRates 0.01 0.1 0.2 0.5 1.0 --pyvolveSim --seqgenSim --indelibleSim

parser = argparse.ArgumentParser(description='Compare simulators, in particular phastSim, seq-gen, indelible and pyvolve.')
parser.add_argument('--path',default="", help='Path where to run simulations.')
parser.add_argument('--reference',default="MN908947.3.fasta", help='File containing the reference genome to be used as root genome. To be found in the folder specified with --path. By default the reference is generated randomly.')
parser.add_argument("--nLeaves", help="Ignore tree file and simulate a random tree with the given number of leaves.", type=int, default=0)
parser.add_argument("--replicates", help="Number of replicate simulations to run", type=int, default=1)
parser.add_argument("--length", help="Length of the genome, if not read from file", type=int, default=29903)
parser.add_argument('--scale',default=1.0,type=float, help='Scale the simulation tree by this amount (default 1.0). Branch lengths are assumed in terms of expected substitutions per site (more or less, as frequencies changes through time, total mutation rate might also  change).')
parser.add_argument("--seed", help="Seed for random simulator in genes simulations", type=int, default=1)
parser.add_argument("--alpha", help="Parameter of the gamma distribution for rate variation; each site will then have a separate rate. If specified, continuous rate variation is assumed, otherwise homogeneous rates are used unless the --mutationRates optioon is used.", type=float, default=0.0)
parser.add_argument("--invariable", help="Proportion of invariable sites, that is, sites that have a mutation rate of 0.0 .", type=float, default=0.0)
parser.add_argument("--mutationRates", help="Mutation rates, by default using the neutral rates estimated from SARS-CoV-2; so far only exactly 12 input values allowed (r_AC, r_AG, r_AT, r_CA, etc...) corresponding to an UNREST nucleotide substitution model.", type=float, nargs='+',  default=[0.039, 0.310, 0.123, 0.140, 0.022, 3.028, 0.747, 0.113, 2.953, 0.056, 0.261, 0.036])
parser.add_argument("--codon", help="Run simulations under a codon model, where mutation rates are used to describe nucleotide mutation rates, and omegas are used to describe the effect of selection at the amino acid level. Default is false (uses a nucleotide substitution model).", action="store_true")
parser.add_argument("--omegaAlpha", help="Parameter of the gamma distribution for omega variation; each codon will then have a separate omega. If specified, continuous omega variation is assumed, otherwise homogeneous omegas are used unless the --omegaCategoryRates option is used.", type=float, default=0.0)
parser.add_argument("--omegaCategoryProbs", help="Probabilities of omega categories. They are supposed to sum up to one, but if they don't they are normalized to do so. By default only 1 category is simulated.", type=float, nargs='+',  default=[1.0])
parser.add_argument("--omegaCategoryRates", help="Omegas of different omega categories. The overall evolutionary rate is renormalized so that the expected number of substitutions per branch length unit per site is 1. By default only 1 category of rate 1.0 is simulated. The number of omegas has to be the same as the number of omega category probabilities, otherwise an error is thrown.", type=float, nargs='+',  default=[1.0])
parser.add_argument("--categoryProbs", help="Probabilities of rate categories. They are supposed to sum up to one, but if they don't they are normalized to do so. By default only 1 category is simulated.", type=float, nargs='+',  default=[1.0])
parser.add_argument("--categoryRates", help="Rates of site categories. The overall mutation rate is renormalized so that the expected number of bstitutions per branch length unit per site is 1. By default only 1 category of rate 1.0 is simulated. The number of rates has to be the same as the number of categgory probabilities, otherwise an error is thrown.", type=float, nargs='+',  default=[1.0])
parser.add_argument("--noHierarchy", help="Run without hierarchical algorithm; the latter is faster with more complex models.", action="store_true")
#parser.add_argument("--createFasta", help="Create a fasta file with the simulated genomes (default name sars-cov-2_simulation_output.fasta).", action="store_true")
parser.add_argument("--pyvolveSim", help="run simulations using pyvolve for comparison)", action="store_true")
parser.add_argument("--seqgenSim", help="run simulations using seqgen for comparison)", action="store_true")
parser.add_argument("--indelibleSim", help="run simulations using indelible for comparison)", action="store_true")
parser.add_argument("--indelibleSim2", help="run simulations using indelible method 2 for comparison)", action="store_true")
args = parser.parse_args()

pathSimu=args.path
reference=args.reference
nLeaves=args.nLeaves
nReplicates=args.replicates
np.random.seed(args.seed)
seed=args.seed
scale=args.scale
alpha=args.alpha
invariable=args.invariable
mutationRates=args.mutationRates
codon=args.codon
omegaAlpha=args.omegaAlpha
omegaCategoryProbs=args.omegaCategoryProbs
omegaCategoryRates=args.omegaCategoryRates
categoryProbs=args.categoryProbs
categoryRates=args.categoryRates
length=args.length

#createFasta=args.createFasta
pyvolveSim=args.pyvolveSim
seqgenSim=args.seqgenSim
indelibleSim=args.indelibleSim
indelibleSim2=args.indelibleSim2

hierarchy=not args.noHierarchy

#collect reference
ref=""
if reference!="":
	file=open(pathSimu+reference)
	line=file.readline()
	while line!="":
		line=file.readline()
		ref+=line.replace("\n","")
	file.close()
	print("\n Finished reading reference genome at "+pathSimu+reference+" with "+str(len(ref))+" bases.")
else:
	print("random ancestor not implemente yet!")
	exit()

#SARS-CoV-2 genome annotation - not used yet but will be useful when simulating under a codon model.
geneEnds=[[266,13468],[13468,21555],[21563,25384],[25393,26220],[26245,26472],[26523,27191],[27202,27387],[27394,27759],[27894,28259],[28274,29533],[29558,29674]]
if codon:
	print("Extracting and concatenating CDSs for codon model")
	newRef=""
	for g in geneEnds:
		newRef+=ref[g[0]-1:g[1]]
	ref=newRef
	print("new Ref length: "+str(len(ref)))
	fileRef=open(pathSimu+reference+"_new.fa","w")
	fileRef.write(">reference\n"+ref+"\n")
	fileRef.close()
	reference=reference+"_new.fa"

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

if alpha<0.000000001:
	sum=0.0
	for i in categoryProbs:
		sum+=i
	for i in range(nCat):
		categoryProbs[i]=categoryProbs[i]/sum
	if sum>1.000001 or sum<0.999999:
		print("\n Normalizing probabilities of site categories. New probabilities:")
		print(categoryProbs)
	


def branch_lengths_2_decimals(str_newick_tree):
    """replaces branch lengths in scientific notation with decimals"""
    colon_s = 0
    #comma_back_paren_s = 0
    num = ''
    new_tree = ''
    for count, char in enumerate(str_newick_tree):
        if char == ':':
            colon_s = count
            continue
        if char in (')', ','):
            #comma_back_paren_s = 1
            #num = '%f' % float(num)
            #num=str(float(num))
            num="{:.10f}".format(float(num))
            new_tree += ":" + num
            colon_s = 0
            num = ''
        if colon_s != 0:
            num = num + char
        if colon_s == 0:
            new_tree += char
    new_tree = new_tree.strip('\'').strip('\"').strip('\'') + ";"
    return new_tree
	

times=[]
for r in range(nReplicates):
	times2=[]
	#simulate tree
	start = time.time()
	treeFile="simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+".tree"
	treeFile2="simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_2.tree"
	if nLeaves<=1:
		#use dendropy - too slow
		from dendropy.simulate import treesim
		t = treesim.birth_death_tree(birth_rate=length, death_rate=0.0, num_total_tips=nLeaves)
		tString=t.as_string(schema="newick")
		tString=tString.split()[1]
		tString=tString.replace(":0)",":1.0e-09)").replace(":0,",":1.0e-09,")
		tString2="("+tString.replace(";","")+");"
		tString3=branch_lengths_2_decimals(tString2.replace(";",""))
	elif nLeaves<1:
		#use ngesh - too slow
		import ngesh
		#print(ngesh.gen_tree.__doc__) 
		tree = ngesh.gen_tree(length, 0.0, min_leaves=nLeaves, labels="enum",seed=seed+r)
		tString=tree.write()
		tString=tString.replace(")1",")")
		#print(tString)
		tString2=tString
		tString3=branch_lengths_2_decimals(tString2.replace(";",""))
		#print(tString3)
		#exit()
	else:
		#use custom script - works well even for huge phylogenies
		import random_tree
		tree = random_tree.gen_tree(length, 0.0*float(length), min_leaves=nLeaves, labels="enum",seed=seed+r)
		tString=tree.write()
		tString=tString.replace(")1",")")
		#tString=tString.replace(":0)",":1.0e-09)").replace(":0,",":1.0e-09,")
		print(tString)
		tString2=tString
		tString3=branch_lengths_2_decimals(tString2.replace(";",""))

	file=open(pathSimu+treeFile,"w")
	file.write(tString+"\n")
	file.close()
	file=open(pathSimu+treeFile2,"w")
	file.write(tString2+"\n")
	file.close()
	#t.write(path=pathSimu+treeFile, schema="newick")
	time1 = time.time() - start
	times2.append(time1)
	print("Time for generating tree: "+str(time1))
	

	
	#run phastSim
	print("Running phastSim")
	start = time.time()
	
	stringRun="python3 /Users/demaio/Desktop/coronavirus/simulations/scripts/efficientSimuSARS2.py --path /Users/demaio/Desktop/coronavirus/simulations/  --seed "+str(seed+r)+" --scale "+str(scale)+" --treeFile "+treeFile+" --reference "+reference+" --outputFile simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+".txt"
	#if createFasta:
	#	stringRun+=" --createFasta "
	if alpha>=0.000000001:
		stringRun+=(" --alpha "+str(alpha))
	elif len(categoryProbs)>1:
		stringRun+=" --categoryProbs "
		for c in range(len(categoryProbs)):
			stringRun+=" "+str(categoryProbs[c])
		stringRun+=" --categoryRates "
		for c in range(len(categoryProbs)):
			stringRun+=" "+str(categoryRates[c])
	if invariable>=0.000000001:
		stringRun+=(" --invariable "+str(invariable))
	if codon:
		stringRun+=" --codon "
		if omegaAlpha>=0.000000001:
			stringRun+=(" --omegaAlpha "+str(omegaAlpha))
		elif len(omegaCategoryProbs)>1:
			stringRun+=" --omegaCategoryProbs "
			for c in range(len(omegaCategoryProbs)):
				stringRun+=" "+str(omegaCategoryProbs[c])
			stringRun+=" --omegaCategoryRates "
			for c in range(len(omegaCategoryRates)):
				stringRun+=" "+str(omegaCategoryRates[c])
	os.system(stringRun)
	
	time2 = time.time() - start
	times2.append(time2)
	print("Total time after simulating sequence evolution along tree with phastSim: "+str(time2))
	
	start = time.time()
	os.system(stringRun+" --createFasta ")
	time2 = time.time() - start
	times2.append(time2)
	print("Total time after simulating sequence evolution along tree with phastSim and writing fasta file: "+str(time2))
	
	#run pyvolve
	#Run simulations using pyvolve if requested
	if pyvolveSim:
		print("Running pyvolve")
		start = time.time()
		
		#import pyvolve
		#pyvolveTree = pyvolve.read_tree(file = pathSimu+treeFile, scale_tree = args.scale)
		#pyvolveTree = pyvolve.read_tree(tree = tString2, scale_tree = args.scale)
		#nucModel=pyvolve.Model("custom", {"matrix":mutMatrix}, alpha = 0.5, num_categories = len(categoryRates))
		#partitions=pyvolve.Partition(models = nucModel, root_sequence = ref)
		#my_evolver = pyvolve.Evolver(tree = pyvolveTree, partitions = partitions)
		#my_evolver(seqfile=pathSimu+"simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_pyvolveOutput.txt")
		
		stringRun="python2 /Users/demaio/Desktop/coronavirus/simulations/scripts/runPyvolve.py --path /Users/demaio/Desktop/coronavirus/simulations/  --seed "+str(seed+r)+" --scale "+str(scale)+" --treeFile "+treeFile2+" --outputFile simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_pyvolve.txt"
		if alpha>=0.000000001 or invariable>=0.000000001:
			print("I am not allowing continuoous rate variation in pyvolve!")
		elif len(categoryProbs)>1:
			stringRun+=" --categoryProbs "
			for c in range(len(categoryProbs)):
				stringRun+=" "+str(categoryProbs[c])
			stringRun+=" --categoryRates "
			for c in range(len(categoryProbs)):
				stringRun+=" "+str(categoryRates[c])
			os.system(stringRun)
			elapsedTime = time.time() - start
			times2.append(elapsedTime)
			print("Time for simulating with pyvolve: "+str(elapsedTime))



	#run simulations with seqgen
	if seqgenSim:
		print("Running seq-gen")
		start = time.time()
		if alpha>=0.000000001:
			categoryString="-a "+str(alpha)+" "
		elif len(categoryProbs)>1:
			categoryString="-g "+str(len(categoryRates))+" "
		else:
			categoryString=" "
		if invariable>=0.000000001:
			categoryString+=(" -i "+str(invariable)+" ")
		stringRun="/Applications/Seq-Gen-1.3.4/source/seq-gen -l "+str(length)+" -s "+str(scale)+" -m GTR -a 0.5 -f 0.3 0.2 0.2 0.3 -r "+str(mutMatrix[0][1]/mutMatrix[2][3])+" "+str(mutMatrix[0][2]/mutMatrix[2][3])+" "+str(mutMatrix[0][3]/mutMatrix[2][3])+" "+str(mutMatrix[1][2]/mutMatrix[2][3])+" "+str(mutMatrix[1][3]/mutMatrix[2][3])+" 1.0 "+categoryString+" -z "+str(seed+r)+" "+pathSimu+treeFile+" > "+pathSimu+"simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_seqgenOutput.txt"
		os.system(stringRun)
		time2 = time.time() - start
		times2.append(time2)
		print("Total time after simulating sequence evolution along tree with seqgen: "+str(time2))
	
	
	
	#run simulations with indelible
	if indelibleSim:
		print("Running indelible method 1")
		start = time.time()
		
		file=open(pathSimu+"control.txt","w")
		file.write("[TYPE] NUCLEOTIDE 1\n")
		file.write("[MODEL]    modelname\n")
		file.write("  [submodel]     UNREST "+str(mutMatrix[3][1]/mutMatrix[2][0])+" "+str(mutMatrix[3][0]/mutMatrix[2][0])+" "+str(mutMatrix[3][2]/mutMatrix[2][0])+" "+str(mutMatrix[1][3]/mutMatrix[2][0])+" "+str(mutMatrix[1][0]/mutMatrix[2][0])+" "+str(mutMatrix[1][2]/mutMatrix[2][0])+" "+str(mutMatrix[0][3]/mutMatrix[2][0])+" "+str(mutMatrix[0][1]/mutMatrix[2][0])+" "+str(mutMatrix[0][2]/mutMatrix[2][0])+" "+str(mutMatrix[2][3]/mutMatrix[2][0])+" "+str(mutMatrix[2][1]/mutMatrix[2][0])+"  \n")
		file.write("  [statefreq] 0.3 0.2 0.2 0.3 \n")
		if alpha>=0.000000001:
			file.write("  [rates] "+str(invariable)+" "+str(alpha)+" 0 \n")
		elif len(categoryProbs)>1:
			file.write("  [rates] "+str(invariable)+" 1 "+str(len(categoryProbs))+" \n")
		file.write("[TREE] treename  "+tString3+" \n")
		file.write("[PARTITIONS] partitionname [treename modelname "+str(length)+"]  \n")
		file.write("[EVOLVE] partitionname 1 simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_indelibleOutput.txt \n")
		file.close()

		os.system("cd "+pathSimu+" ; "+"/Users/demaio/Desktop/INDELibleV1.03/bin/indelible ")
		
		time2 = time.time() - start
		times2.append(time2)
		print("Total time after simulating sequence evolution along tree with indelible method 1: "+str(time2))
		
	if indelibleSim2:
		#now version with full Gillespie
		print("Running indelible method 2")
		start = time.time()
		
		file=open(pathSimu+"control.txt","w")
		file.write("[TYPE] NUCLEOTIDE 2\n")
		file.write("[MODEL]    modelname\n")
		file.write("  [submodel]     UNREST "+str(mutMatrix[3][1]/mutMatrix[2][0])+" "+str(mutMatrix[3][0]/mutMatrix[2][0])+" "+str(mutMatrix[3][2]/mutMatrix[2][0])+" "+str(mutMatrix[1][3]/mutMatrix[2][0])+" "+str(mutMatrix[1][0]/mutMatrix[2][0])+" "+str(mutMatrix[1][2]/mutMatrix[2][0])+" "+str(mutMatrix[0][3]/mutMatrix[2][0])+" "+str(mutMatrix[0][1]/mutMatrix[2][0])+" "+str(mutMatrix[0][2]/mutMatrix[2][0])+" "+str(mutMatrix[2][3]/mutMatrix[2][0])+" "+str(mutMatrix[2][1]/mutMatrix[2][0])+"  \n")
		file.write("  [statefreq] 0.3 0.2 0.2 0.3 \n")
		if alpha>=0.000000001:
			file.write("  [rates] "+str(invariable)+" "+str(alpha)+" 0 \n")
		elif len(categoryProbs)>1:
			file.write("  [rates] "+str(invariable)+" 1 "+str(len(categoryProbs))+" \n")
		file.write("[TREE] treename  "+tString3+" \n")
		file.write("[PARTITIONS] partitionname [treename modelname "+str(length)+"]  \n")
		file.write("[EVOLVE] partitionname 1 simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_indelible2Output.txt \n")
		file.close()

		os.system("cd "+pathSimu+" ; "+"/Users/demaio/Desktop/INDELibleV1.03/bin/indelible ")
		
		time2 = time.time() - start
		times2.append(time2)
		print("Total time after simulating sequence evolution along tree with indelible method 2: "+str(time2))
		
	
	times.append(times2)

print("Total matrix of running times:")
print(times)

exit()
























