#import sys
import os
#import math
import numpy as np
import os.path
#from os import path
import argparse
#from ete3 import Tree
import time

#script to compare the running time of phastSim with pyvolve, seq-gen and indelible.

# python2 compareSimulators.py --nLeaves 100 --createFasta --pyvolveSim --seqgenSim --indelibleSim --indelibleSim2

parser = argparse.ArgumentParser(description='Compare simulators, in particular phastSim, seq-gen, indelible and pyvolve.')
#parser.add_argument('--path',default="/Users/demaio/Desktop/coronavirus/simulations/", help='Path where to write simulation output.')
#parser.add_argument('--reference',default="/Users/demaio/Documents/GitHub/phastSim/phastSim/example/MN908947.3.fasta", help='File containing the reference genome to be used as root genome.')
parser.add_argument('--path',default="/home/will/Desktop/projects/embl/phastSim/simulation_output_6/", help='Path where to write simulation output.')
parser.add_argument('--reference',default="/home/will/Desktop/projects/embl/phastSim/phastSim/example/MN908947.3.fasta", help='File containing the reference genome to be used as root genome.')
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
parser.add_argument("--pyvolveSim", help="run simulations using pyvolve for comparison", action="store_true")
parser.add_argument("--seqgenSim", help="run simulations using seqgen for comparison", action="store_true")
parser.add_argument("--indelibleSim", help="run simulations using indelible for comparison", action="store_true")
parser.add_argument("--indelibleSim2", help="run simulations using indelible method 2 for comparison", action="store_true")
parser.add_argument("--createFasta", help="also run phastSim with fasta file generation.", action="store_true")
parser.add_argument("--phastSim", help="Run standard phastSim with hierarchical approach.", action="store_true")
parser.add_argument("--generatePlots", help="Generate plots from previous simulation runs", action="store_true")
parser.add_argument("--useIndels", help="Run using indels with 0.1 insertion/deletion rates and a geometric(0.5) distribution.", action="store_true")
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
use_indels=args.useIndels

#createFasta=args.createFasta
pyvolveSim=args.pyvolveSim
seqgenSim=args.seqgenSim
indelibleSim=args.indelibleSim
indelibleSim2=args.indelibleSim2
createFasta=args.createFasta
generatePlots=args.generatePlots
phastSim=args.phastSim

noHierarchy=args.noHierarchy
#hierarchy=not args.noHierarchy

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
else:
	print("random ancestor not implemente yet!")
	exit()

#SARS-CoV-2 genome annotation - not used yet but will be useful when simulating under a codon model.
geneEnds=[[266,13468],[13468,21555],[21563,25384],[25393,26220],[26245,26472],[26523,27191],[27202,27387],[27394,27759],[27894,28259],[28274,29533],[29558,29674]]
if codon:
	print("Extracting and concatenating CDSs for codon model")
	newRef=""
	for g in geneEnds:
		newRef+=ref[g[0]+2:g[1]-3]
	ref=newRef
	print("new Ref length: "+str(len(ref)))
	fileRef=open(reference+"_new.fa","w")
	fileRef.write(">reference\n"+ref+"\n")
	fileRef.close()
	reference=reference+"_new.fa"
	
length=len(ref)
	

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
	
def rescaleTree(node,scale):
	node.dist*=scale
	for c in node.get_children():
		rescaleTree(c,scale)

times=[[],[],[],[],[],[],[],[]]
for r in range(nReplicates):
	times2=[]
	#simulate tree
	start = time.time()
	treeFile="simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+".tree"
	treeFile2="simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_2.tree"

	#use custom script - works well even for huge phylogenies
	import random_tree
	tree = random_tree.gen_tree(length, 0.0*float(length), min_leaves=nLeaves, labels="enum", seed=seed+r)
	rescaleTree(tree,scale)
	tString=tree.write()
	tString=tString.replace(")1",")")
	#tString=tString.replace(":0)",":1.0e-09)").replace(":0,",":1.0e-09,")
	#print(tString)
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
	times[0].append(time1)
	#times2.append(time1)
	print("Time for generating tree: "+str(time1))
	

	
	#run phastSim
	print("Running phastSim")
	
	
	stringRun="phastSim --outpath "+pathSimu+"  --seed "+str(seed+r)+" --treeFile "+pathSimu+treeFile+" --reference "+reference+" --outputFile simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+".txt" # --scale "+str(scale)+"
	#if createFasta:
	#	stringRun+=" --createFasta "
	if use_indels:
		stringRun += " --indels --insertionRate CONSTANT 0.1 --deletionRate CONSTANT 0.1 --insertionLength GEOMETRIC 0.5 --deletionLength GEOMETRIC 0.5 " 
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
	if phastSim:
		start = time.time()
		os.system(stringRun)#+" >/dev/null")
	
		time2 = time.time() - start
		times[1].append(time2)
		#times2.append(time2)
		print("Total time after simulating sequence evolution along tree with phastSim: "+str(time2))
	else:
		times[1].append("NaN")
	
	if createFasta:
		start = time.time()
		os.system(stringRun+" --createFasta "+" >/dev/null")
		time2 = time.time() - start
		#times2.append(time2)
		times[2].append(time2)
		print("Total time after simulating sequence evolution along tree with phastSim and writing fasta file: "+str(time2))
	else:
		times[2].append("NaN")
		
	if noHierarchy:
		start = time.time()
		os.system(stringRun+" --noHierarchy "+" >/dev/null")
		time2 = time.time() - start
		#times2.append(time2)
		times[3].append(time2)
		print("Total time after simulating sequence evolution along tree with phastSim and no multilayer genome tree: "+str(time2))
	else:
		times[3].append("NaN")
	
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
		
		stringRun="python /home/will/Desktop/projects/embl/phastSim/scripts/runPyvolve.py --path /home/will/Desktop/projects/embl/phastSim/simulation_output_6/  --seed "+str(seed+r)+" --scale "+str(scale)+" --treeFile "+treeFile2+" --outputFile simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_pyvolve.txt"

		if alpha>=0.000000001 or invariable>=0.000000001:
			print("I am not allowing continuous rate variation in pyvolve!")
			times[3].append("NaN")
		else:
			if len(categoryProbs)>1:
				stringRun+=" --categoryProbs "
				for c in range(len(categoryProbs)):
					stringRun+=" "+str(categoryProbs[c])
				stringRun+=" --categoryRates "
				for c in range(len(categoryProbs)):
					stringRun+=" "+str(categoryRates[c])
			os.system(stringRun)
			time2 = time.time() - start
			#times2.append(elapsedTime)
			times[4].append(time2)
			print("Time for simulating with pyvolve: "+str(time2))
	else:
		#times2.append("NaN")
		times[4].append("NaN")



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
		#stringRun="/Applications/Seq-Gen-1.3.4/source/seq-gen -l "+str(length)+" -m GTR -f 0.3 0.2 0.2 0.3 -r "+str(mutMatrix[0][1]/mutMatrix[2][3])+" "+str(mutMatrix[0][2]/mutMatrix[2][3])+" "+str(mutMatrix[0][3]/mutMatrix[2][3])+" "+str(mutMatrix[1][2]/mutMatrix[2][3])+" "+str(mutMatrix[1][3]/mutMatrix[2][3])+" 1.0 "+categoryString+" -z "+str(seed+r)+" "+pathSimu+treeFile+" > "+pathSimu+"simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_seqgenOutput.txt" # -s "+str(scale)+"
		stringRun="seq-gen -l "+str(length)+" -m GTR -f 0.3 0.2 0.2 0.3 -r "+str(mutMatrix[0][1]/mutMatrix[2][3])+" "+str(mutMatrix[0][2]/mutMatrix[2][3])+" "+str(mutMatrix[0][3]/mutMatrix[2][3])+" "+str(mutMatrix[1][2]/mutMatrix[2][3])+" "+str(mutMatrix[1][3]/mutMatrix[2][3])+" 1.0 "+categoryString+" -z "+str(seed+r)+" "+pathSimu+treeFile+" > "+pathSimu+"simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_seqgenOutput.txt" # -s "+str(scale)+"
		os.system(stringRun+" ") #&>/dev/null
		time2 = time.time() - start
		#times2.append(time2)
		times[5].append(time2)
		print("Total time after simulating sequence evolution along tree with seqgen: "+str(time2))
	else:
		#times2.append("NaN")
		times[5].append("NaN")
	
	
	#run simulations with indelible
	if indelibleSim:
		print("Running indelible method 1")
		start = time.time()
		
		file=open(pathSimu+"control.txt","w")
		if codon:
			file.write("[TYPE] CODON 1\n")
			file.write("[MODEL]    modelname\n")
			file.write("  [submodel]  2.0 ")
			for c in range(len(omegaCategoryProbs)-1):
				file.write(" "+str(omegaCategoryProbs[c]))
			for c in range(len(omegaCategoryRates)):
				file.write(" "+str(omegaCategoryRates[c]))
			file.write("  \n")
			if use_indels:
				file.write(" [insertmodel]	NB	0.5	1 \n")
				file.write(" [deletemodel]	NB	0.5	1 \n")
				file.write(" [insertrate]    0.1 \n")
				file.write(" [deleterate]    0.1 \n")
			
			file.write("[TREE] treename  "+tString3+" \n")
			file.write("[PARTITIONS] partitionname [treename modelname "+str(int(length/3))+"]  \n")
			file.write("[EVOLVE] partitionname 1 simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_codon_indelibleOutput.txt \n")
		else:
			file.write("[TYPE] NUCLEOTIDE 1\n")
			file.write("[MODEL]    modelname\n")
			file.write("  [submodel]     UNREST "+str(mutMatrix[3][1]/mutMatrix[2][0])+" "+str(mutMatrix[3][0]/mutMatrix[2][0])+" "+str(mutMatrix[3][2]/mutMatrix[2][0])+" "+str(mutMatrix[1][3]/mutMatrix[2][0])+" "+str(mutMatrix[1][0]/mutMatrix[2][0])+" "+str(mutMatrix[1][2]/mutMatrix[2][0])+" "+str(mutMatrix[0][3]/mutMatrix[2][0])+" "+str(mutMatrix[0][1]/mutMatrix[2][0])+" "+str(mutMatrix[0][2]/mutMatrix[2][0])+" "+str(mutMatrix[2][3]/mutMatrix[2][0])+" "+str(mutMatrix[2][1]/mutMatrix[2][0])+"  \n")
			file.write("  [statefreq] 0.3 0.2 0.2 0.3 \n")
			if use_indels:
				file.write(" [insertmodel]	NB	0.5	1 \n")
				file.write(" [deletemodel]	NB	0.5	1 \n")
				file.write(" [insertrate]    0.1 \n")
				file.write(" [deleterate]    0.1 \n")
			if alpha>=0.000000001:
				file.write("  [rates] "+str(invariable)+" "+str(alpha)+" 0 \n")
			elif len(categoryProbs)>1:
				file.write("  [rates] "+str(invariable)+" 1 "+str(len(categoryProbs))+" \n")
			file.write("[TREE] treename  "+tString3+" \n")
			file.write("[PARTITIONS] partitionname [treename modelname "+str(length)+"]  \n")
			file.write("[EVOLVE] partitionname 1 simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_indelibleOutput.txt \n")
		file.close()

		#os.system("cd "+pathSimu+" ; "+"/Users/demaio/Desktop/INDELibleV1.03/bin/indelible ")
		os.system("cd "+pathSimu+" ; "+"/home/will/Downloads/INDELibleV1.03/bin/indelible ")
		
		time2 = time.time() - start
		#times2.append(time2)
		times[6].append(time2)
		print("Total time after simulating sequence evolution along tree with indelible method 1: "+str(time2))
	else:
		#times2.append("NaN")
		times[6].append("NaN")
		
	if indelibleSim2:
		#now version with full Gillespie
		print("Running indelible method 2")
		start = time.time()
		
		file=open(pathSimu+"control.txt","w")
		if codon:
			file.write("[TYPE] CODON 2\n")
			file.write("[MODEL]    modelname\n")
			file.write("  [submodel]  2.0 ")
			for c in range(len(omegaCategoryProbs)-1):
				file.write(" "+str(omegaCategoryProbs[c]))
			for c in range(len(omegaCategoryRates)):
				file.write(" "+str(omegaCategoryRates[c]))
			file.write("  \n")
			if use_indels:
				file.write(" [insertmodel]	NB	0.5	1 \n")
				file.write(" [deletemodel]	NB	0.5	1 \n")
				file.write(" [insertrate]    0.1 \n")
				file.write(" [deleterate]    0.1 \n")
			
			file.write("[TREE] treename  "+tString3+" \n")
			file.write("[PARTITIONS] partitionname [treename modelname "+str(int(length/3))+"]  \n")
			file.write("[EVOLVE] partitionname 1 simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_codon_indelible2Output.txt \n")
		else:
			file.write("[TYPE] NUCLEOTIDE 2\n")
			file.write("[MODEL]    modelname\n")
			file.write("  [submodel]     UNREST "+str(mutMatrix[3][1]/mutMatrix[2][0])+" "+str(mutMatrix[3][0]/mutMatrix[2][0])+" "+str(mutMatrix[3][2]/mutMatrix[2][0])+" "+str(mutMatrix[1][3]/mutMatrix[2][0])+" "+str(mutMatrix[1][0]/mutMatrix[2][0])+" "+str(mutMatrix[1][2]/mutMatrix[2][0])+" "+str(mutMatrix[0][3]/mutMatrix[2][0])+" "+str(mutMatrix[0][1]/mutMatrix[2][0])+" "+str(mutMatrix[0][2]/mutMatrix[2][0])+" "+str(mutMatrix[2][3]/mutMatrix[2][0])+" "+str(mutMatrix[2][1]/mutMatrix[2][0])+"  \n")
			file.write("  [statefreq] 0.3 0.2 0.2 0.3 \n")
			if use_indels:
				file.write(" [insertmodel]	NB	0.5	1 \n")
				file.write(" [deletemodel]	NB	0.5	1 \n")
				file.write(" [insertrate]    0.1 \n")
				file.write(" [deleterate]    0.1 \n")
			if alpha>=0.000000001:
				file.write("  [rates] "+str(invariable)+" "+str(alpha)+" 0 \n")
			elif len(categoryProbs)>1:
				file.write("  [rates] "+str(invariable)+" 1 "+str(len(categoryProbs))+" \n")
			file.write("[TREE] treename  "+tString3+" \n")
			file.write("[PARTITIONS] partitionname [treename modelname "+str(length)+"]  \n")
			file.write("[EVOLVE] partitionname 1 simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_indelible2Output.txt \n")
		file.close()

		#os.system("cd "+pathSimu+" ; "+"/Users/demaio/Desktop/INDELibleV1.03/bin/indelible ")
		os.system("cd "+pathSimu+" ; "+"/home/will/Downloads/INDELibleV1.03/bin/indelible ")
		
		time2 = time.time() - start
		#times2.append(time2)
		times[7].append(time2)
		print("Total time after simulating sequence evolution along tree with indelible method 2: "+str(time2))
	else:
		#times2.append("NaN")
		times[7].append("NaN")
		
	
	#times.append(times2)

print("Total matrix of running times:")
print(times)


if generatePlots:
	
	#import matplotlib.pyplot as plt
	
	import matplotlib.pyplot as plt
	#import numpy as np
	#from matplotlib.patches import Polygon
	
	names=["tree generation","phastSim","phastSim+Fasta","pyvolve","SeqGen","INDELible-m1","INDELible-m2"]
	colors=["blue","red","orange","green","purple","yellow","brown"]
	# the commented out values were used previously (on Nicola's PC).

	times10=[[0.286470890045166, 0.31896352767944336, 0.3150522708892822, 0.35256433486938477, 0.34905338287353516, 0.3493950366973877, 0.3511476516723633, 0.3582336902618408, 0.351032018661499, 0.35051774978637695], [2.0049617290496826, 2.0820443630218506, 1.8900294303894043, 2.2931203842163086, 2.27139949798584, 2.2659008502960205, 2.3040614128112793, 2.22983717918396, 2.2670042514801025, 2.3632993698120117], [1.9599967002868652, 1.9465186595916748, 1.9787395000457764, 2.040527582168579, 2.0918385982513428, 2.0643293857574463, 2.1026904582977295, 2.1438252925872803, 2.1532113552093506, 2.0645127296447754], [6.93420672416687, 7.5393946170806885, 8.075521469116211, 8.473002433776855, 8.589826822280884, 8.422191381454468, 8.41925859451294, 8.49124026298523, 8.545444250106812, 8.691659450531006 ], [0.012601613998413086, 0.01311635971069336, 0.01289224624633789, 0.012771129608154297, 0.013012886047363281, 0.012801170349121094, 0.012769460678100586, 0.012870073318481445, 0.012779474258422852, 0.012956619262695312], [0.03798532485961914, 0.03840351104736328, 0.040375709533691406, 0.04167771339416504, 0.04137134552001953, 0.041043758392333984, 0.042273759841918945, 0.042439937591552734, 0.04085898399353027, 0.04187369346618652],[0.02474522590637207, 0.024831533432006836, 0.02729034423828125, 0.028174638748168945, 0.027623414993286133, 0.028365373611450195, 0.027507305145263672, 0.028249740600585938, 0.028221845626831055, 0.02847886085510254]]
	#[[0.32402634620666504, 0.3214292526245117, 0.3188958168029785, 0.28540635108947754, 0.3108482360839844, 0.31295108795166016, 0.2913055419921875, 0.2925238609313965, 0.2932114601135254, 0.3199167251586914], [1.8718585968017578, 1.862605333328247, 1.9202218055725098, 1.9187707901000977, 1.9224002361297607, 1.9175093173980713, 1.9194765090942383, 1.9422886371612549, 2.1560351848602295, 2.136245012283325], [1.980574131011963, 2.0036814212799072, 1.969557762145996, 1.9862523078918457, 1.914921760559082, 1.999527931213379, 1.9412376880645752, 2.090512275695801, 2.105055809020996, 2.1017000675201416], [0.5465364456176758, 0.5536282062530518, 0.5508022308349609, 0.5633573532104492, 0.5786206722259521, 0.5846483707427979, 0.5441484451293945, 0.5838699340820312, 0.5855481624603271, 0.5608088970184326 ], [0.015013694763183594, 0.014559507369995117, 0.015006542205810547, 0.014104843139648438, 0.014327049255371094, 0.01415395736694336, 0.01450490951538086, 0.014388561248779297, 0.014214515686035156, 0.025898456573486328], [0.0422518253326416, 0.044155120849609375, 0.038712263107299805, 0.04276752471923828, 0.040401458740234375, 0.04051089286804199, 0.04031658172607422, 0.04380369186401367, 0.03952741622924805, 0.07233691215515137],[0.02580547332763672, 0.02463388442993164, 0.023540973663330078, 0.02529001235961914, 0.024139404296875, 0.024132251739501953, 0.026221036911010742, 0.023634672164916992, 0.02323293685913086, 0.047487497329711914]]
	#[[0.0013659000396728516, 0.0012090206146240234, 0.0013148784637451172, 0.0010118484497070312, 0.0013380050659179688, 0.002196073532104492, 0.0014619827270507812, 0.001085042953491211, 0.001199960708618164, 0.0013890266418457031], [1.3747148513793945, 1.3423030376434326, 1.3215079307556152, 1.3170280456542969, 1.3145039081573486, 1.3163139820098877, 1.3323869705200195, 1.3344120979309082, 1.399022102355957, 1.3535301685333252], [1.3641538619995117, 1.3581798076629639, 1.334787130355835, 1.3404159545898438, 1.3326139450073242, 1.3254010677337646, 1.3380720615386963, 1.3444411754608154, 1.3227300643920898, 1.3296558856964111], [10.451174020767212, 10.58082389831543, 10.541486978530884, 10.444769144058228, 10.375016927719116, 10.394562005996704, 10.366996049880981, 10.622323989868164, 10.495352983474731, 10.429331064224243], [0.03925800323486328, 0.032778024673461914, 0.033091068267822266, 0.03274106979370117, 0.036557912826538086, 0.032247066497802734, 0.03708791732788086, 0.041892051696777344, 0.038690805435180664, 0.035842180252075195], [0.21403002738952637, 0.21619105339050293, 0.21575593948364258, 0.21455693244934082, 0.2209010124206543, 0.221451997756958, 0.21570682525634766, 0.2141880989074707, 0.22292304039001465, 0.21918296813964844], [0.12611007690429688, 0.1258530616760254, 0.1236422061920166, 0.12346792221069336, 0.12584614753723145, 0.12546515464782715, 0.12303709983825684, 0.12221598625183105, 0.1220400333404541, 0.12607908248901367]]
	times20=[[0.3538219928741455, 0.35256242752075195, 0.35068202018737793, 0.35391855239868164, 0.35389089584350586, 0.35297393798828125, 0.3555889129638672, 0.35291576385498047, 0.35629749298095703, 0.3530271053314209], [2.2835533618927, 2.3008270263671875, 2.2627522945404053, 2.2767410278320312, 2.2627642154693604, 2.2968878746032715, 2.2799906730651855, 2.2735390663146973, 2.3126041889190674, 2.281782865524292], [2.125321865081787, 2.1067090034484863, 2.1070659160614014, 2.0969557762145996, 2.1393401622772217, 2.1330881118774414, 2.131077289581299, 2.133711814880371, 2.1346688270568848, 2.1495466232299805], [18.068254232406616, 17.70702052116394, 17.693846225738525, 18.588911294937134, 18.20507311820984, 17.854225158691406, 17.68247675895691, 18.212438583374023, 18.805498838424683, 18.89807891845703 ], [0.0244901180267334, 0.023113489151000977, 0.023440837860107422, 0.02306842803955078, 0.02300238609313965, 0.023348331451416016, 0.023157835006713867, 0.0234224796295166, 0.023848533630371094, 0.023596525192260742], [0.07850837707519531, 0.07745718955993652, 0.08003401756286621, 0.07822990417480469, 0.07800412178039551, 0.08019804954528809, 0.07890987396240234, 0.07845377922058105, 0.09174537658691406, 0.0785372257232666],[0.05437588691711426, 0.052553653717041016, 0.05224275588989258, 0.05171704292297363, 0.052687644958496094, 0.05399322509765625, 0.05250859260559082, 0.05260801315307617, 0.053426265716552734, 0.05310535430908203]]
	#[[0.0024929046630859375, 0.0014770030975341797, 0.0023810863494873047, 0.0014848709106445312, 0.001817941665649414, 0.0016529560089111328, 0.0020999908447265625, 0.0015180110931396484, 0.0020749568939208984, 0.0017440319061279297], [1.645247220993042, 1.310811996459961, 1.3231310844421387, 1.3270468711853027, 1.3239338397979736, 1.3379931449890137, 1.3164989948272705, 1.3218259811401367, 1.3193409442901611, 1.5228931903839111], [1.679689884185791, 1.3187479972839355, 1.3616609573364258, 1.329289197921753, 1.346217155456543, 1.3264689445495605, 1.3758349418640137, 1.3410701751708984, 1.3292009830474854, 1.6160540580749512], [22.576815128326416, 22.07684302330017, 22.30865716934204, 21.877304792404175, 21.605812072753906, 21.84848117828369, 22.243227005004883, 22.178505897521973, 22.18565797805786, 23.35649800300598], [0.06273508071899414, 0.05921602249145508, 0.0638580322265625, 0.05792999267578125, 0.05753588676452637, 0.05878496170043945, 0.057292938232421875, 0.06512904167175293, 0.05781197547912598, 0.057282209396362305], [0.44101905822753906, 0.44517993927001953, 0.43025708198547363, 0.4448831081390381, 0.44520998001098633, 0.44223713874816895, 0.43360185623168945, 0.4268670082092285, 0.46094799041748047, 0.4530940055847168], [0.24215412139892578, 0.24007606506347656, 0.24377012252807617, 0.23820710182189941, 0.24094796180725098, 0.2389540672302246, 0.24043488502502441, 0.23656988143920898, 0.29044508934020996, 0.25960588455200195]]
	times50=[[0.3528478145599365, 0.3800690174102783, 0.36669373512268066, 0.3984487056732178, 0.3549675941467285, 0.3582484722137451, 0.35726261138916016, 0.3565552234649658, 0.35793638229370117, 0.357755184173584], [2.3104307651519775, 2.2700436115264893, 2.3246967792510986, 2.6904139518737793, 2.3094537258148193, 2.293480396270752, 2.3130524158477783, 2.3382277488708496, 2.3402392864227295, 2.328073263168335], [2.1109278202056885, 2.1701135635375977, 2.154858350753784, 2.6594481468200684, 2.1045942306518555, 2.0899605751037598, 2.154416561126709, 2.140125274658203, 2.167069911956787, 2.1244428157806396],  [46.586344957351685, 48.117431640625, 53.37920546531677, 53.35781145095825, 49.02303409576416, 47.120819330215454, 46.65156292915344, 46.184895277023315, 46.85412263870239, 45.58630681037903 ], [0.052538156509399414, 0.05588245391845703, 0.05894160270690918, 0.05544638633728027, 0.052907705307006836, 0.0561068058013916, 0.052176475524902344, 0.05303549766540527, 0.053055524826049805, 0.05276012420654297], [0.19448256492614746, 0.21820569038391113, 0.23572468757629395, 0.2046794891357422, 0.19811749458312988, 0.19909405708312988, 0.1802988052368164, 0.18567991256713867, 0.1891040802001953, 0.18378615379333496],[0.13672709465026855, 0.14298367500305176, 0.19878196716308594, 0.1411442756652832, 0.14324355125427246, 0.14290475845336914, 0.13991236686706543, 0.1418323516845703, 0.13859987258911133, 0.13866901397705078]]
	#[[0.003216981887817383, 0.0038487911224365234, 0.0030329227447509766, 0.003144979476928711, 0.004025936126708984, 0.002969980239868164, 0.0034058094024658203, 0.002958059310913086, 0.023120880126953125, 0.02884507179260254], [1.2893130779266357, 1.6244549751281738, 1.486210823059082, 1.4786460399627686, 1.4211831092834473, 1.811263084411621, 1.376521110534668, 1.4128239154815674, 1.5160620212554932, 1.3880500793457031], [1.3925518989562988, 1.4420089721679688, 1.4690890312194824, 1.5120949745178223, 1.4923410415649414, 1.4661610126495361, 1.4662230014801025, 1.4426720142364502, 1.4030230045318604, 1.4014220237731934], [61.306437969207764, 63.28899002075195, 60.32176899909973, 61.62099599838257, 62.32091784477234, 60.96144986152649, 60.207597970962524, 59.269126892089844, 60.00672507286072, 59.49119710922241], [0.13117289543151855, 0.13464593887329102, 0.13071918487548828, 0.1415119171142578, 0.13554096221923828, 0.13779211044311523, 0.1376659870147705, 0.13895797729492188, 0.13930082321166992, 0.1320171356201172], [1.1521000862121582, 1.117422103881836, 1.1462361812591553, 1.1094839572906494, 1.1246919631958008, 1.1122980117797852, 1.109248161315918, 1.1435189247131348, 1.1563069820404053, 1.1218719482421875], [0.6041860580444336, 0.6159520149230957, 0.6443181037902832, 0.6238491535186768, 0.5989768505096436, 0.5900440216064453, 0.609623908996582, 0.5931220054626465, 0.6229169368743896, 0.600348949432373]]
	times100=[[0.3605227470397949, 0.36771273612976074, 0.3620142936706543, 0.3610813617706299, 0.3643789291381836, 0.36153626441955566, 0.36297011375427246, 0.3598511219024658, 0.3606250286102295, 0.36702919006347656], [2.324413776397705, 2.3448030948638916, 2.2683048248291016, 2.355851173400879, 2.302262306213379, 2.614875316619873, 2.308723211288452, 2.2713937759399414, 2.577224016189575, 2.2992186546325684], [2.1528897285461426, 2.1498656272888184, 2.152667760848999, 2.1706840991973877, 2.124654531478882, 2.1586265563964844, 2.1820075511932373, 2.1325128078460693, 2.1518514156341553, 2.1527624130249023], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], [0.10081744194030762, 0.10376214981079102, 0.10110592842102051, 0.10138201713562012, 0.10301685333251953, 0.10176539421081543, 0.10157108306884766, 0.10219907760620117, 0.10280966758728027, 0.1014702320098877], [0.38815999031066895, 0.372699499130249, 0.37894487380981445, 0.36351490020751953, 0.37461280822753906, 0.3726005554199219, 0.369903564453125, 0.3683609962463379, 0.38035154342651367, 0.37946295738220215],[0.2601146697998047, 0.3069155216217041, 0.29543209075927734, 0.29478955268859863, 0.29495882987976074, 0.2943120002746582, 0.3089559078216553, 0.2960972785949707, 0.3049640655517578, 0.31397151947021484]]
	#[[0.009551048278808594, 0.0053958892822265625, 0.005692958831787109, 0.00590205192565918, 0.005836963653564453, 0.0058939456939697266, 0.005366086959838867, 0.005799770355224609, 0.006151914596557617, 0.0052869319915771484], [1.90081787109375, 1.3287749290466309, 1.3533220291137695, 1.4254589080810547, 1.5417671203613281, 1.4508280754089355, 1.4114508628845215, 1.549408197402954, 1.5193719863891602, 1.400057077407837], [1.3768410682678223, 1.385685920715332, 1.412989854812622, 1.3838729858398438, 1.5289270877838135, 1.4544739723205566, 1.4677560329437256, 1.51531982421875, 1.4512879848480225, 1.4553320407867432], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [0.2564408779144287, 0.2659909725189209, 0.26821398735046387, 0.2528820037841797, 0.26355600357055664, 0.2617809772491455, 0.25960516929626465, 0.25537586212158203, 0.24987196922302246, 0.25585007667541504], [2.20414400100708, 2.1862261295318604, 2.231003999710083, 2.1976799964904785, 2.2304821014404297, 2.2315008640289307, 2.573899984359741, 2.246011972427368, 2.2906219959259033, 2.20377516746521], [1.1515510082244873, 1.1431560516357422, 1.1758790016174316, 1.146554946899414, 1.1703548431396484, 1.178049087524414, 1.3542709350585938, 1.2229149341583252, 1.154651165008545, 1.1641769409179688]]
	times200=[[0.3697977066040039, 0.37154245376586914, 0.3704705238342285, 0.3203010559082031, 0.30208826065063477, 0.31093406677246094, 0.30785632133483887, 0.3038299083709717, 0.31195592880249023, 0.30602264404296875], [2.33345365524292, 2.293391466140747, 2.343384265899658, 2.0045413970947266, 2.0226736068725586, 2.0051918029785156, 1.9962029457092285, 2.0674495697021484, 1.9853968620300293, 2.0697474479675293], [2.1809158325195312, 2.202779769897461, 2.2099366188049316, 2.1083574295043945, 2.101482391357422, 2.0778238773345947, 2.071920156478882, 2.155670404434204, 2.0888094902038574, 2.1042933464050293], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], [0.19618630409240723, 0.20128369331359863, 0.22835230827331543, 0.18909144401550293, 0.20229721069335938, 0.2088613510131836, 0.21081781387329102, 0.21525144577026367, 0.2302243709564209, 0.20631861686706543], [0.7894492149353027, 0.8540546894073486, 0.7978489398956299, 0.8405978679656982, 0.6944313049316406, 0.7361471652984619, 0.7282116413116455, 0.7236196994781494, 0.7994036674499512, 0.7387704849243164],[0.5849626064300537, 0.5522007942199707, 0.6045231819152832, 0.5013446807861328, 0.5505819320678711, 0.5447268486022949, 0.5406420230865479, 0.5656862258911133, 0.596081018447876, 0.5491166114807129]]
	#[[0.016273975372314453, 0.011254072189331055, 0.011461973190307617, 0.010603904724121094, 0.013769149780273438, 0.013958930969238281, 0.026243925094604492, 0.014862060546875, 0.011723995208740234, 0.010759115219116211], [1.9284470081329346, 1.7488529682159424, 1.5611369609832764, 1.5293529033660889, 1.6704580783843994, 1.5298941135406494, 1.548914909362793, 1.8947179317474365, 1.6275060176849365, 1.5622169971466064], [1.761125087738037, 1.7953181266784668, 1.5314640998840332, 1.6465811729431152, 1.7108330726623535, 1.6326282024383545, 1.6300709247589111, 1.7247998714447021, 1.568471908569336, 2.247066020965576], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [0.4969940185546875, 0.5056090354919434, 0.4887361526489258, 0.5002598762512207, 0.49825382232666016, 0.49005913734436035, 0.4913029670715332, 0.5112559795379639, 0.49417805671691895, 0.4985392093658447], [5.168694019317627, 5.151010990142822, 5.2190330028533936, 5.296245098114014, 5.013715982437134, 5.4771199226379395, 6.150067090988159, 5.349929094314575, 4.977146148681641, 5.350111961364746], [2.8624649047851562, 2.833103895187378, 2.7423510551452637, 2.6613221168518066, 2.663487195968628, 2.860327959060669, 3.3309261798858643, 2.802933931350708, 2.646231174468994, 2.517848014831543]]
	times500=[[0.30927610397338867, 0.36544275283813477, 0.3450932502746582, 0.30956244468688965, 0.3121151924133301, 0.3526463508605957, 0.3074028491973877, 0.31971120834350586, 0.31170654296875, 0.3200862407684326], [2.0508828163146973, 2.178354024887085, 2.012525796890259, 2.0926146507263184, 2.0488452911376953, 2.1748151779174805, 2.0666205883026123, 2.0955605506896973, 2.0279133319854736, 2.0116307735443115], [2.1624059677124023, 2.1972615718841553, 2.1707282066345215, 2.1957273483276367, 2.122082233428955, 2.1553404331207275, 2.2179973125457764, 2.2125117778778076, 2.262425661087036, 2.249150514602661], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], [0.4888169765472412, 0.4863319396972656, 0.5340523719787598, 0.5313558578491211, 0.47763800621032715, 0.528071403503418, 0.5774078369140625, 0.5229134559631348, 0.5575599670410156, 0.4910142421722412], [1.9245662689208984, 1.9973838329315186, 1.9835631847381592, 2.0207395553588867, 2.0040602684020996, 1.9497828483581543, 2.151925802230835, 1.9248974323272705, 2.159038782119751, 2.052276611328125],[1.3047397136688232, 1.4080557823181152, 1.3182940483093262, 1.3443877696990967, 1.3176193237304688, 1.3390007019042969, 1.270117998123169, 1.3929383754730225, 1.2255518436431885, 1.3066596984863281]]
	#[[0.030160903930664062, 0.02469801902770996, 0.04264998435974121, 0.03104686737060547, 0.03468894958496094, 0.02983403205871582, 0.02934098243713379, 0.024353981018066406, 0.02561211585998535, 0.03130292892456055], [1.777817964553833, 1.9888050556182861, 1.9227111339569092, 1.9003281593322754, 2.0147860050201416, 1.9665369987487793, 1.8426148891448975, 1.4138050079345703, 1.4871950149536133, 1.5098040103912354], [1.5789639949798584, 1.5292949676513672, 1.564305067062378, 1.691133975982666, 1.7239899635314941, 1.560417890548706, 1.5549790859222412, 1.6500108242034912, 1.5650379657745361, 1.653033971786499], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [1.222322940826416, 1.2483329772949219, 1.2385191917419434, 1.2624180316925049, 1.2348461151123047, 1.3084070682525635, 1.2472999095916748, 1.2439031600952148, 1.2595350742340088, 1.2943921089172363], [11.081032991409302, 11.013811111450195, 11.17693018913269, 11.416763067245483, 11.85645604133606, 11.123574018478394, 11.086908102035522, 11.227660894393921, 11.142872095108032, 11.275111198425293], [5.7495949268341064, 5.812233924865723, 6.02468204498291, 5.983781099319458, 6.056072950363159, 5.8660888671875, 5.7771899700164795, 5.78557014465332, 5.831812143325806, 6.089737892150879]]
	times1000=[[0.3621346950531006, 0.38665223121643066, 0.31900787353515625, 0.3224220275878906, 0.32614874839782715, 0.35962581634521484, 0.324965238571167, 0.327625036239624, 0.33636903762817383, 0.3233802318572998], [2.039370536804199, 2.212015151977539, 2.207923650741577, 2.194664239883423, 2.1344387531280518, 2.1163933277130127, 2.0266501903533936, 2.0396029949188232, 2.0094895362854004, 2.040476083755493], [2.405404806137085, 2.44035005569458, 2.4309840202331543, 2.4775707721710205, 2.377474784851074, 2.3623545169830322, 2.3917644023895264, 2.379338502883911, 2.421314001083374, 2.461298942565918], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], [1.0367393493652344, 1.062631368637085, 1.0485081672668457, 1.0177793502807617, 1.046640396118164, 1.045102834701538, 1.0521225929260254, 1.1114771366119385, 1.2765841484069824, 1.1065080165863037], [3.6736600399017334, 3.914623737335205, 3.9768152236938477, 4.0130455493927, 4.228492975234985, 3.9983980655670166, 3.9995248317718506, 4.538363456726074, 3.946955442428589, 4.009037017822266],[2.484671115875244, 2.7940917015075684, 2.73533296585083, 2.740504741668701, 3.156716823577881, 2.6919333934783936, 2.76507306098938, 2.6453404426574707, 2.7502517700195312, 2.656099557876587]]
	#[[0.06780815124511719, 0.06437397003173828, 0.05380582809448242, 0.04997706413269043, 0.05129718780517578, 0.04825186729431152, 0.054646968841552734, 0.04891395568847656, 0.05314779281616211, 0.049278974533081055], [1.4967939853668213, 1.9933040142059326, 1.5371029376983643, 1.5284991264343262, 1.9865779876708984, 1.506382942199707, 1.9183871746063232, 2.0416319370269775, 1.5037190914154053, 1.8662559986114502], [1.849992036819458, 1.8258349895477295, 1.709947109222412, 1.6919009685516357, 1.8142428398132324, 1.6936531066894531, 1.898118019104004, 1.7096149921417236, 1.7108089923858643, 1.7036280632019043], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [2.4352219104766846, 2.4571030139923096, 2.4070231914520264, 2.4291841983795166, 2.451385974884033, 2.439326047897339, 2.4405691623687744, 2.484835147857666, 2.471928834915161, 2.465049982070923], [22.25348210334778, 22.035884857177734, 21.825459957122803, 21.92672896385193, 21.81256103515625, 22.11038112640381, 21.884013175964355, 21.707130908966064, 21.84363889694214, 21.670482873916626], [11.619686126708984, 11.556863069534302, 11.478086948394775, 11.657796144485474, 11.494728088378906, 11.577374935150146, 11.585293054580688, 11.68454885482788, 11.430992126464844, 11.45626187324524]]
	times2000=[[0.3661062717437744, 0.45351290702819824, 0.3601338863372803, 0.35436058044433594, 0.35435938835144043, 0.3600795269012451, 0.36786413192749023, 0.35814714431762695, 0.3976435661315918, 0.3498270511627197], [2.160513162612915, 2.2566988468170166, 2.1803650856018066, 2.1846585273742676, 2.1303911209106445, 2.098470449447632, 2.11624813079834, 2.1036572456359863, 2.1532790660858154, 2.1683804988861084], [2.5164432525634766, 2.90380859375, 2.7777934074401855, 2.9818320274353027, 2.9727835655212402, 2.8844592571258545, 2.8335378170013428, 2.7332332134246826, 2.742670774459839, 2.761235237121582], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], [1.9231503009796143, 2.260104179382324, 2.2015814781188965, 2.1370723247528076, 2.3042519092559814, 2.2724063396453857, 2.228992223739624, 2.184084892272949, 2.1104650497436523, 2.2324934005737305], [7.6945130825042725, 8.943911790847778, 8.64449143409729, 8.120613813400269, 7.753779649734497, 7.903039932250977, 8.951663732528687, 8.138253211975098, 8.089648008346558, 8.02677035331726],[5.074307918548584, 4.969040870666504, 5.2580201625823975, 5.27686333656311, 5.446438789367676, 5.543627977371216, 5.2642505168914795, 5.467250108718872, 5.298097848892212, 5.33390736579895]]
	#[[0.12939810752868652, 0.099761962890625, 0.10983014106750488, 0.137160062789917, 0.11310815811157227, 0.12170195579528809, 0.17563486099243164, 0.10024499893188477, 0.0971071720123291, 0.10719680786132812], [1.6819391250610352, 2.0618398189544678, 2.0402579307556152, 2.0742900371551514, 1.892535924911499, 2.2402780055999756, 1.8665268421173096, 1.5126521587371826, 1.5043201446533203, 1.7153770923614502], [2.1061370372772217, 2.091276168823242, 2.240056037902832, 2.9656341075897217, 2.3872971534729004, 2.5221140384674072, 2.6080589294433594, 2.0910890102386475, 2.0963239669799805, 2.1601290702819824], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [4.886102914810181, 4.879283905029297, 4.931666135787964, 4.81746506690979, 4.882803916931152, 4.861746072769165, 4.885030031204224, 4.975928068161011, 4.911868095397949, 4.934772968292236], [44.09892296791077, 44.005027055740356, 45.88066387176514, 48.627665996551514, 56.57271981239319, 52.096412897109985, 44.491461992263794, 43.75139808654785, 44.21598196029663, 44.35627007484436], [23.274376153945923, 23.85594391822815, 25.263535976409912, 25.600807905197144, 26.603089809417725, 27.35308003425598, 23.306721925735474, 23.075037002563477, 23.742963075637817, 23.250370979309082]]
	times5000=[[0.4319624900817871, 0.578153133392334, 0.4791405200958252, 0.4583749771118164, 0.4581429958343506, 0.4563584327697754, 0.4969463348388672, 0.4782121181488037, 0.4441094398498535, 0.476804256439209], [2.4047598838806152, 2.4544565677642822, 2.3045332431793213, 2.407986640930176, 2.400059223175049, 2.4558863639831543, 2.297607898712158, 2.3934261798858643, 2.3565797805786133, 2.357121706008911], [3.359138250350952, 4.676752805709839, 4.2872374057769775, 4.083299875259399, 4.012995481491089, 4.027643203735352, 3.901313543319702, 4.1031951904296875, 4.5778303146362305, 4.204830884933472], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], [5.241696834564209, 6.00919508934021, 5.414457082748413, 5.3550708293914795, 5.577221870422363, 5.836879253387451, 5.897825002670288, 5.351659297943115, 5.994544744491577, 5.44526481628418], [21.11788034439087, 24.148942947387695, 22.35642910003662, 22.02048683166504, 20.632789611816406, 20.679664373397827, 20.647603273391724, 21.54020667076111, 22.911720037460327, 20.413681030273438],[13.770885229110718, 15.842254638671875, 13.566912412643433, 13.837852239608765, 13.728861808776855, 14.39827275276184, 14.470169067382812, 14.848873376846313, 14.513187646865845, 14.163492679595947]]
	#[[0.31500887870788574, 0.24887490272521973, 0.3096740245819092, 0.2508058547973633, 0.3484029769897461, 0.2693758010864258, 0.27234315872192383, 0.35924696922302246, 0.3092200756072998, 0.42699503898620605], [1.9736239910125732, 2.429014205932617, 1.9887659549713135, 2.0318288803100586, 2.406982183456421, 2.0246009826660156, 1.9706370830535889, 2.4449448585510254, 2.524501085281372, 2.53228497505188], [3.344974994659424, 3.331270933151245, 3.323568820953369, 3.300321102142334, 3.405689001083374, 3.3196661472320557, 3.3752760887145996, 3.3163578510284424, 3.508978843688965, 3.707275867462158], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [12.1112220287323, 12.05104112625122, 12.257277965545654, 12.403314113616943, 12.445502042770386, 12.404994010925293, 12.645019054412842, 12.529222965240479, 12.521138906478882, 12.774863004684448], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [59.62119007110596, 58.707270860672, 58.55328392982483, 59.17231893539429, 63.55593514442444, 59.023497104644775, 58.952008962631226, 59.784512996673584, 63.09826993942261, 60.914510011672974]]
	times10000=[[0.6260547637939453, 0.7074122428894043, 0.6357431411743164, 0.6991472244262695, 0.639310359954834, 0.6259925365447998, 0.710153341293335, 0.7380797863006592, 0.6671605110168457, 0.6532645225524902], [2.7210512161254883, 2.8169045448303223, 2.749976634979248, 2.844238519668579, 2.6944046020507812, 2.797699213027954, 2.9618775844573975, 2.79940128326416, 2.9193601608276367, 2.9056859016418457], [5.002784729003906, 6.24365234375, 6.210233688354492, 6.346195936203003, 6.405362606048584, 6.158026456832886, 6.951162099838257, 6.973858118057251, 6.434947967529297, 6.413838624954224],  ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], [9.8789381980896, 12.122819423675537, 11.20814847946167, 12.080686807632446, 11.188852071762085, 11.808539152145386, 12.18217158317566, 11.679875373840332, 12.776806116104126, 11.668919801712036], [41.407557249069214, 41.701902866363525, 40.675811529159546, 43.204975605010986, 41.2561297416687, 42.822746992111206, 40.41104793548584, 43.9643976688385, 43.87178349494934, 40.83157801628113],[26.497570037841797, 26.034415245056152, 26.1755313873291, 26.574292421340942, 26.59775996208191, 27.335813760757446, 27.705915451049805, 26.716097831726074, 27.45615243911743, 27.301084280014038]]
	#[[0.5868499279022217, 0.580380916595459, 0.6529080867767334, 0.6462969779968262, 0.6883280277252197, 0.6795670986175537, 0.6926059722900391, 0.4921681880950928, 0.7163400650024414, 0.6498241424560547], [2.942147970199585, 2.4787280559539795, 2.444871187210083, 2.4249610900878906, 2.5257160663604736, 2.4837939739227295, 2.465178966522217, 2.3728160858154297, 2.475499153137207, 2.415358066558838], [5.290570974349976, 5.2421088218688965, 5.169349908828735, 5.3277130126953125, 5.383677005767822, 5.364545106887817, 5.328459978103638, 5.3075220584869385, 5.364189147949219, 5.305814981460571], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [23.960747957229614, 23.89430594444275, 24.102982997894287, 24.04936909675598, 24.217010021209717, 24.367241144180298, 24.017180919647217, 24.123044967651367, 24.17894697189331, 23.99029302597046], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	times20000=[[1.0456008911132812, 1.0157508850097656, 0.9977357387542725, 1.0433409214019775, 0.9844036102294922, 0.9949188232421875, 1.0466749668121338, 1.0282049179077148, 1.0322685241699219, 1.0468966960906982], [3.4624462127685547, 3.320486307144165, 3.2932329177856445, 3.3226077556610107, 3.3741281032562256, 3.3315465450286865, 3.404574155807495, 3.3680167198181152, 3.4742331504821777, 3.5063624382019043], [7.9548726081848145, 13.762882471084595, 10.170494794845581, 10.22055435180664, 10.322172403335571, 13.923986673355103, 11.7287437915802, 10.132718563079834, 11.774081707000732, 10.1583092212677],  ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], [20.55828619003296, 27.191519737243652, 22.930957794189453, 22.4307644367218, 22.116400003433228, 25.666804313659668, 23.80967879295349, 21.890169143676758, 22.126861572265625, 22.27155375480652], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'],['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	#[[1.1991140842437744, 1.327165126800537, 1.4817399978637695, 1.4684429168701172, 1.7353639602661133, 1.6934490203857422, 1.5031070709228516, 1.3446290493011475, 1.4934370517730713, 1.8493130207061768], [3.3477210998535156, 3.769974946975708, 3.3797521591186523, 3.4234910011291504, 3.576843023300171, 3.5666489601135254, 3.539125919342041, 3.4996449947357178, 3.4938271045684814, 3.4616119861602783], [9.40938401222229, 9.229989051818848, 9.341056823730469, 9.800310134887695, 9.614714860916138, 9.365905046463013, 9.440688848495483, 9.506255149841309, 9.528976917266846, 9.752545833587646], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [48.57287001609802, 48.325536012649536, 48.499929904937744, 48.55190300941467, 48.723459005355835, 48.555269956588745, 48.50492286682129, 48.813239097595215, 48.42090392112732, 48.982256174087524], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	times50000=[[2.3065271377563477, 2.629556655883789, 2.3764169216156006, 2.29233717918396, 2.3040945529937744, 2.3308117389678955, 2.256986618041992, 2.2236831188201904, 2.3057613372802734, 2.2760653495788574], [5.730430841445923, 5.980284214019775, 5.472230911254883, 5.573661804199219, 5.59115743637085, 5.612621784210205, 5.6475443840026855, 5.602144002914429, 5.522136449813843, 5.645097255706787], [18.477766752243042, 22.74521493911743, 18.818898916244507, 19.26700782775879, 18.477845430374146, 19.159831523895264, 18.33698344230652, 21.301982641220093, 20.87391686439514, 19.02156972885132],  ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'],['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	#[[3.293437957763672, 4.06860089302063, 3.951677083969116, 4.423721075057983, 3.8627521991729736, 4.419393062591553, 4.29544997215271, 4.469686031341553, 4.260180950164795, 4.0369179248809814], [6.553931951522827, 6.685858964920044, 6.578627824783325, 6.476964950561523, 6.540050029754639, 6.447605133056641, 6.578918933868408, 7.231057167053223, 6.6752049922943115, 6.765127897262573], [21.58997893333435, 20.89152193069458, 21.59639310836792, 21.675976037979126, 21.34534788131714, 22.414045095443726, 25.244637966156006, 22.464181900024414, 21.908604860305786, 22.150094032287598], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	times100000=[[4.7787184715271, 4.990846633911133, 4.7829906940460205, 4.697843551635742, 5.73615574836731, 5.246704339981079, 4.7617716789245605, 4.723233938217163, 4.6709020137786865, 4.681627035140991], [8.839689254760742, 9.733030080795288, 9.74084758758545, 9.130806922912598, 8.701268434524536, 9.257972478866577, 8.856314897537231, 8.71550726890564, 9.100442886352539, 9.361463785171509], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'],  ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'],['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	#[[7.647463083267212, 8.008918046951294, 9.324841022491455, 9.470612049102783, 8.557217121124268, 9.787823915481567, 9.107623100280762, 9.38093113899231, 8.35722804069519, 8.31895899772644], [11.86309003829956, 11.798933982849121, 12.203128099441528, 11.762664079666138, 12.412771940231323, 12.195564985275269, 11.996767044067383, 11.904637813568115, 11.997853994369507, 12.214859962463379], [42.35268306732178, 40.621132135391235, 42.493098974227905, 42.65987491607666, 44.19890284538269, 42.4167218208313, 43.31222414970398, 43.10220289230347, 43.47725200653076, 43.2533118724823], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	times200000=[[9.742886304855347, 9.687723398208618, 9.805734395980835, 10.160772800445557, 9.70213007926941, 10.00439453125, 10.242988586425781, 10.556771278381348, 9.81380295753479, 9.89162302017212], [15.708847045898438, 16.888503789901733, 16.527430057525635, 15.489632844924927, 16.673904418945312, 16.86138343811035, 17.193645477294922, 16.379803895950317, 15.804017782211304, 15.780726909637451], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'],['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	#[[15.817861080169678, 16.58775305747986, 19.115757942199707, 16.761270999908447, 19.467989206314087, 17.577764987945557, 17.07235097885132, 19.266084909439087, 16.83725094795227, 18.492412090301514], [21.678630113601685, 21.98533320426941, 22.04712176322937, 22.029582977294922, 22.8069269657135, 21.75388789176941, 21.656525135040283, 21.720540046691895, 21.85643482208252, 21.551968097686768], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	times500000=[[26.701364755630493, 26.59220552444458, 25.53991198539734, 30.91444158554077, 26.915411949157715, 26.484868049621582, 25.6319682598114, 25.93476366996765, 25.34960389137268, 25.988577127456665], [37.68199133872986, 39.63857102394104, 37.82759428024292, 40.06772065162659, 39.95154285430908, 39.504962682724, 36.5359206199646, 38.00857353210449, 38.43317937850952, 38.274829149246216], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'],['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	#[[41.40035009384155, 43.9114990234375, 55.79817795753479, 47.97320508956909, 54.219762086868286, 53.92414689064026, 46.42072582244873, 51.166882038116455, 45.05335998535156, 51.811765909194946], [54.970056772232056, 55.55358910560608, 56.28190588951111, 56.577943086624146, 59.569891929626465, 59.95699906349182, 65.13326716423035, 53.935949087142944, 52.92773199081421, 55.46649694442749], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	#times1000000=[]

	#nLeaves=[10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000]
	nLeaves=["10","20","50","100","200","500","1000","2000","5000","10^4","2x10^4","5x10^4","10^5","2x10^5","5x10^5"]
	times=[times10,times20,times50,times100,times200,times500,times1000,times2000,times5000,times10000,times20000,times50000,times100000,times200000,times500000]
	
	
	#BOXPlot drawing
	def boxplot(valuesLists,axisLabels,plotFileName,labels,colors, xLabel,topPlot,degreeSkew=45):
	
		fig, ax1 = plt.subplots(figsize=(10, 6))
		fig.canvas.set_window_title('Simulation Running times')
		fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
		#num_boxes = len(valuesLists)*len(valuesLists[0])
		
		data=np.zeros((len(valuesLists)*len(valuesLists[0]),len(valuesLists[0][0])))
		print((len(valuesLists)*len(valuesLists[0]),len(valuesLists[0][0])))
		positions=[]
		for i in range(len(valuesLists)):
			for j in range(len(valuesLists[i])):
				c=colors[j]
				positions.append(i+(j-len(valuesLists[i])/2.0)*0.5/len(valuesLists[i]))
				for k in range(len(valuesLists[i][j])):
					#print("")
					#print(valuesLists[i][j][k])
					#print(data[i*len(valuesLists[0])+j][k])
					data[i*len(valuesLists[0])+j][k]=valuesLists[i][j][k]
				position=[i+(j-(len(valuesLists[i])-1)/2.0)*0.5/len(valuesLists[i])]
				#if len(valuesLists[i])%2==0:
				#	position=[i+(j-(len(valuesLists[i])-1)/2.0)*0.5/len(valuesLists[i])]
				#else:
				#	position=[i+(j-(len(valuesLists[i])-1)/2.0)*0.5/len(valuesLists[i])]
				plt.boxplot(data[i*len(valuesLists[0])+j], positions=position, notch=False, patch_artist=True, widths=0.5/len(colors), manage_ticks=False, 
					boxprops=dict(facecolor=c, color=c),
					capprops=dict(color=c),
					whiskerprops=dict(color=c),
					medianprops=dict(color=c),
					flierprops = dict(marker='o', markerfacecolor=c, markersize=1.5,
					  linestyle='none', markeredgecolor=c)
				)
		
		# Add a horizontal grid to the plot, but make it very light in color
		# so we can use it for reading data values but not be distracting
		ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',   alpha=0.5)

		# Hide these grid behind plot objects
		ax1.set_axisbelow(True)
		ax1.set_title('Comparison of Simulation Running times')
		ax1.set_xlabel(xLabel)
		ax1.set_ylabel('Time (seconds)')

		# Set the axes ranges and axes labels
		ax1.set_xlim(-0.7, len(valuesLists) -0.5)
		top = topPlot
		bottom = -1
		ax1.set_ylim(bottom, top)
		#ax1.set_xticklabels(labels, rotation=45, fontsize=8)
		ax1.set_xticks(np.arange(len(valuesLists)))
		ax1.set_xticklabels(labels, rotation=degreeSkew, fontsize=11)

		# Finally, add a basic legend
		for i in range(len(colors)):
			fig.text(0.10, 0.85-0.045*i, axisLabels[i], backgroundcolor=colors[i], color='black', weight='roman', size='medium')
		
		fig.savefig(plotFileName)
		plt.close()

	def errplot(times,labels,plotFileName,n_leaves,colors,topPlot,degreeSkew=45):
		
		mean_times = []
		errors = []
		for t1 in times:
			mean_times.append([])
			errors.append([])
			for t2 in t1:
				mean_times[-1].append(np.mean([float(t) for t in t2]))
				errors[-1].append(np.std([float(t) for t in t2]))

		mean_times = np.array(mean_times).T
		errors = np.array(errors).T    
		

		x = n_leaves
		y = range(100,200)
		fig = plt.figure(figsize=(15, 9))
		ax1 = fig.add_subplot(111)
		# Hide these grid behind plot objects
		ax1.set_axisbelow(True)
		ax1.set_title('Comparison of Simulation Running Times')
		ax1.set_xlabel(topPlot)
		ax1.set_ylabel('Time (seconds)')
		ax1.set_xscale('log')
		#ax1.set_yscale('log')
		ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',   alpha=0.5)
		#ax1.errorbar(x, mean_times[0], yerr=errors[0], marker='x', c='blue', label='tree generation', fmt="o", capsize=2)
		
		reordered_numbers = [0, 1, 2, 4, 6, 5, 3]
		for i in range(len(labels)):
			ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='x', c=colors[i], label=labels[i], capsize=2)
		
		#ax1.errorbar(x, mean_times[1], yerr=errors[1], marker='x', c='red', label='phastSim', capsize=2)
		#ax1.errorbar(x, mean_times[2], yerr=errors[2], marker='x', c='orange', label='phastSim+fasta', capsize=2)
		#ax1.errorbar(x, mean_times[4], yerr=errors[4], marker='x', c='purple', label='Seq-Gen', capsize=2)
		#ax1.errorbar(x, mean_times[6], yerr=errors[6], marker='x', c='brown', label='INDELible-m2', capsize=2)
		#ax1.errorbar(x, mean_times[5], yerr=errors[5], marker='x', c='yellow', label='INDELible-m1', capsize=2)
		#ax1.errorbar(x, mean_times[3], yerr=errors[3], marker='x', c='green', label='pyvolve', capsize=2)
		plt.legend(loc='upper left')
		#plt.show()
		fig.savefig(plotFileName)

	#generate boxplot of general running times
	boxplot(times,names,pathSimu+"boxplot_times_general.pdf",nLeaves,colors,'Number of tips',70)
	nLeaves = [10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000]
	errplot(times,names,pathSimu+"errplot_times_general.pdf",nLeaves,colors,'Number of tips',70)




	#generate boxplot for bacterial simulations
	names=["phastSim","phastSim vanilla","SeqGen"]
	colors=["red","orange","purple"]
	
	timeBac100=[[0.006983041763305664, 0.011324882507324219, 0.009025096893310547, 0.011911869049072266, 0.012537002563476562, 0.0067310333251953125, 0.012540102005004883, 0.005934953689575195, 0.01127314567565918, 0.005726814270019531], [174.96500897407532, 174.36507487297058, 175.25401997566223, 175.28451919555664, 171.74607491493225, 171.94691801071167, 173.5567581653595, 172.69711184501648, 191.58557987213135, 175.2613389492035], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [13.424217939376831, 14.320651054382324, 14.301777839660645, 14.55683422088623, 14.60419511795044, 14.411143064498901, 14.294666051864624, 13.91373586654663, 13.604182958602905, 13.315825939178467], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [34.542505979537964, 37.61758017539978, 37.153862953186035, 37.15251898765564, 36.598602056503296, 36.699177980422974, 36.6919310092926, 36.52304816246033, 36.54633593559265, 34.995906829833984], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	timeBac1000=[[0.05812692642211914, 0.1542520523071289, 0.060903072357177734, 0.06916284561157227, 0.06678295135498047, 0.06496405601501465, 0.05615997314453125, 0.06887292861938477, 0.060105085372924805, 0.057878971099853516], [177.48807787895203, 182.96084809303284, 172.49677300453186, 175.87213397026062, 173.62093901634216, 172.60746812820435, 173.63555788993835, 173.0993959903717, 173.5166938304901, 175.74975085258484], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [13.144881010055542, 14.442451000213623, 13.33045482635498, 14.297406911849976, 12.883018970489502, 12.911453008651733, 12.940874099731445, 12.909605979919434, 12.800656080245972, 12.888622999191284], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [371.2199020385742, 368.92077112197876, 366.3824338912964, 373.98392486572266, 340.9463691711426, 340.8884289264679, 348.04884099960327, 340.19118785858154, 340.25771379470825, 342.55591678619385], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	timeBac10000=[[0.5482289791107178, 0.5563790798187256, 0.6169149875640869, 0.634497880935669, 0.6346278190612793, 0.6516609191894531, 0.6554958820343018, 0.48153090476989746, 0.672562837600708, 0.6116149425506592], [169.9283730983734, 174.08655405044556, 175.7704598903656, 173.69724893569946, 175.6370129585266, 175.27587485313416, 172.99442100524902, 176.0433931350708, 175.20786309242249, 172.27963709831238], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [14.153738021850586, 14.174522876739502, 14.369817018508911, 14.407010078430176, 14.171308040618896, 14.240686178207397, 14.342763900756836, 14.443896055221558, 14.38150691986084, 14.392709016799927], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	timeBac100000=[[7.1389319896698, 7.664942979812622, 8.775321960449219, 7.957852125167847, 8.219660997390747, 8.681435823440552, 7.813961982727051, 9.513015031814575, 8.392756938934326, 8.492002010345459], [183.45800304412842, 181.6222870349884, 179.1631760597229, 186.37382698059082, 180.89060497283936, 179.808434009552, 187.61126399040222, 180.37693405151367, 178.5731999874115, 180.19063591957092], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [27.191282987594604, 26.84230399131775, 27.70342516899109, 29.843647003173828, 28.46603488922119, 28.659116983413696, 29.28634786605835, 28.485222101211548, 29.228497982025146, 28.772555112838745], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	timeBac1000000=[[81.62878894805908, 89.62794303894043, 106.14013600349426, 110.6383581161499, 96.30954599380493, 97.43770599365234, 110.09148812294006, 110.74342083930969, 107.95002603530884, 112.01782989501953], [279.7848780155182, 287.42469906806946, 283.3451671600342, 279.34494400024414, 281.20244097709656, 281.4417998790741, 286.5778498649597, 283.95273303985596, 298.03560304641724, 279.8306269645691], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [174.9824252128601, 172.16432189941406, 170.83757781982422, 168.054692029953, 167.96446204185486, 171.40612292289734, 176.20980095863342, 165.83070993423462, 179.4578800201416, 172.5215358734131], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	
	nLeaves=["100","1000","10^4","10^5","10^6"]
	times=[[timeBac100[1],timeBac100[3],timeBac100[3]],[timeBac1000[1],timeBac1000[3],timeBac1000[5]],[timeBac10000[1],timeBac10000[3],timeBac10000[5]],[timeBac100000[1],timeBac100000[3],timeBac100000[5]],[timeBac1000000[1],timeBac1000000[3],timeBac1000000[5]]]
	
	boxplot(times,names,pathSimu+"boxplot_times_bacteria.pdf",nLeaves,colors,'Number of tips',400)
	nLeaves = [100, 1000, 10000, 100000, 1000000]
	errplot(times,names,pathSimu+"errplot_times_bacteria.pdf",nLeaves,colors,'Number of tips',400)
	
	
	
	
	#generate boxplot for the effects of evolutionary models (rate variation, codon model, etc)
	#Indelible m1 10 tips, indelible m2 100 tips, seq-gen 1000 tips, phastSim 100,000 tips
	#No variation, 10 site classes, continuous distribution, codon model, codon model with 10 site classes, codon model with continuous distributions. 
	INDELible10=[0.2225949764251709, 0.2190871238708496, 0.21699094772338867, 0.22150301933288574, 0.2256309986114502, 0.22266411781311035, 0.22504186630249023, 0.21621108055114746, 0.22166180610656738, 0.22487497329711914]
	INDELible10_10cat=[0.23387789726257324, 0.22591495513916016, 0.22251391410827637, 0.23491692543029785, 0.2268209457397461, 0.22698307037353516, 0.22743487358093262, 0.22669291496276855, 0.223466157913208, 0.22720599174499512]
	INDELible10_alpha=[6.374353885650635, 6.485472917556763, 6.457162141799927, 6.532535076141357, 6.484618902206421, 6.542191028594971, 6.585536956787109, 7.222656965255737, 6.829970121383667, 6.700471878051758]
	INDELible10_codon=[0.20548510551452637, 0.20304203033447266, 0.20437192916870117, 0.21530485153198242, 0.2131330966949463, 0.20339107513427734, 0.2080841064453125, 0.21094799041748047, 0.20357298851013184, 0.2106781005859375]
	INDELible10_codon_10cat=[0.7385458946228027, 0.746147871017456, 0.7383320331573486, 0.7351891994476318, 0.7411129474639893, 0.7182350158691406, 0.7311639785766602, 0.7352509498596191, 0.7086069583892822, 0.7333619594573975]
	
	indelible2_100=[1.178359031677246, 1.162600040435791, 1.1259591579437256, 1.118901014328003, 1.132032871246338, 1.1725618839263916, 1.1374468803405762, 1.1177160739898682, 1.1428818702697754, 1.1853229999542236]
	indelible2_100_10cat=[5.693902969360352, 5.679260969161987, 5.7287890911102295, 5.680241823196411, 5.6857030391693115, 5.7029969692230225, 5.665284872055054, 5.681509971618652, 5.681549072265625, 5.8491740226745605]
	indelible2_100_alpha=[1.1788570880889893, 1.1999309062957764, 1.1450579166412354, 1.139327049255371, 1.1485190391540527, 1.1926980018615723, 1.1499888896942139, 1.1295068264007568, 1.1256649494171143, 1.197829008102417]
	indelible2_100_codon=[0.43164706230163574, 0.4278898239135742, 0.426832914352417, 0.40623998641967773, 0.40640997886657715, 0.42308592796325684, 0.4180629253387451, 0.42006492614746094, 0.4056990146636963, 0.4109060764312744]
	indelible2_100_codon_10cat=[0.43619394302368164, 0.43151402473449707, 0.4308300018310547, 0.44187402725219727, 0.45023298263549805, 0.43655991554260254, 0.4258692264556885, 0.4209752082824707, 0.43906593322753906, 0.4262430667877197]

	seqgen1000=[2.4755289554595947, 2.411590099334717, 2.4759140014648438, 2.3132500648498535, 2.453019142150879, 2.4912800788879395, 2.53383207321167, 2.328788995742798, 2.446610927581787, 2.3403468132019043]
	seqgen1000_10cat=[2.4618189334869385, 2.4114022254943848, 2.459367036819458, 2.3362350463867188, 2.431511878967285, 2.3429830074310303, 2.4409120082855225, 2.3587169647216797, 2.5023059844970703, 2.3628599643707275]
	seqgen1000_alpha=[3.920012950897217, 3.9600830078125, 4.005131006240845, 3.8833510875701904, 3.9653871059417725, 3.9893078804016113, 3.992614984512329, 3.9619109630584717, 3.9531421661376953, 3.959604024887085]
	NaNs=['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']
	
	phastSim20000=[3.8024981021881104, 3.1865439414978027, 3.1835639476776123, 3.0934529304504395, 3.165221929550171, 3.1619739532470703, 3.2406389713287354, 3.072295904159546, 3.0972628593444824, 3.2586851119995117]
	phastSim20000_10cat=[3.226116895675659, 3.116103172302246, 3.2135679721832275, 3.146617889404297, 3.1883559226989746, 3.2167351245880127, 3.1342380046844482, 3.242230176925659, 5.162225008010864, 3.2765140533447266]
	phastSim20000_alpha=[3.089195966720581, 3.113906145095825, 3.1601879596710205, 3.0560100078582764, 3.4628560543060303, 3.225409984588623, 3.119723081588745, 3.0943148136138916, 3.1068480014801025, 3.1178669929504395]
	phastSim20000_codon=[3.292970895767212, 3.2894129753112793, 3.23388409614563, 3.264816999435425, 3.1970441341400146, 3.1810920238494873, 3.1609909534454346, 3.137653112411499, 3.193368911743164, 3.2282588481903076]
	phastSim20000_codon_10cat=[3.237046003341675, 3.310333013534546, 3.2924559116363525, 3.2617690563201904, 3.1413769721984863, 3.296412944793701, 3.244813919067383, 3.250843048095703, 3.157186985015869, 3.184061050415039]
	phastSim20000_codon_alpha=[3.24564790725708, 3.1222660541534424, 3.0887770652770996, 3.0826399326324463, 3.142029047012329, 2.994220018386841, 2.968656063079834, 3.0855040550231934, 3.045428991317749, 3.065749168395996]
	
	names=["phastSim (20,000 tips)","SeqGen (1,000 tips)","INDELible-m1 (10 tips)","INDELible-m2 (100 tips)"]
	colors=["red","purple","yellow","brown"]
	models=["nucleotide","nuc+10cat","nuc+alpha","codon","codon+10cat","codon+alpha"]
	times=[[phastSim20000,seqgen1000,INDELible10,indelible2_100],[phastSim20000_10cat,seqgen1000_10cat,INDELible10_10cat,indelible2_100_10cat],[phastSim20000_alpha,seqgen1000_alpha,INDELible10_alpha,indelible2_100_alpha],[phastSim20000_codon,NaNs,INDELible10_codon,indelible2_100_codon],[phastSim20000_codon_10cat,NaNs,INDELible10_codon_10cat,indelible2_100_codon_10cat],[phastSim20000_codon_alpha,NaNs,NaNs,NaNs]]
	boxplot(times,names,pathSimu+"boxplot_times_models.pdf",models,colors,'Evolutionary model',8)

	#errplot(times,names,pathSimu+"errplot_times_models.pdf",models,colors,'Evolutionary model',8)
	
	
	#checking the effects of branch lengths over running times
	#phastSim 100000 tips, INDELible 1000 tips, Seq-Gen 5000 tips
	phastSim100000=[12.136877059936523, 11.765272855758667, 12.018028974533081, 12.175333023071289, 11.63103699684143, 11.44487714767456, 11.771679878234863, 11.911523818969727, 11.839121103286743, 11.610697984695435]
	phastSim100000_scale01=[6.787471055984497, 6.508502006530762, 6.604642868041992, 6.449320077896118, 6.572448968887329, 6.534075975418091, 6.582338809967041, 6.6191699504852295, 6.599027872085571, 6.590991020202637]
	phastSim100000_scale10=[61.11373996734619, 57.925140142440796, 58.324761152267456, 58.994688987731934, 57.3508939743042, 57.853127002716064, 59.886695861816406, 57.88926315307617, 59.73718214035034, 60.17274498939514]
	
	seqgen5000=[11.902686834335327, 11.977451086044312, 11.844815969467163, 12.060249090194702, 12.185868978500366, 12.155524969100952, 13.859594106674194, 12.145087957382202, 12.049667119979858, 12.08738112449646]
	seqgen5000_scale01=[11.94965410232544, 11.937930822372437, 11.490603923797607, 11.666943073272705, 11.815308809280396, 11.85882306098938, 11.835784196853638, 11.729623079299927, 11.99119520187378, 11.87510895729065]
	seqgen5000_scale10=[11.821383953094482, 11.782198905944824, 12.818541049957275, 12.114749908447266, 12.099608898162842, 11.971168041229248, 11.937159061431885, 12.158576011657715, 11.974681854248047, 12.383517980575562]
	
	indelible1000=[21.826755046844482, 21.699470043182373, 21.36352300643921, 21.440908193588257, 22.4670250415802, 21.947756052017212, 22.008903980255127, 22.392734050750732, 22.055774927139282, 22.210675954818726]
	indelible1000_scale01=[21.488821029663086, 21.76942801475525, 21.745434045791626, 21.717295169830322, 21.533621072769165, 21.723461866378784, 21.697612047195435, 21.479082107543945, 21.48190689086914, 21.40741991996765]
	indelible1000_scale10=[21.98651385307312, 22.482820987701416, 22.185405015945435, 21.92905902862549, 22.216336965560913, 21.83970284461975, 22.006327152252197, 21.74193501472473, 21.808908939361572, 22.07722306251526]
	
	indelible2_1000=[11.328500032424927, 11.311942100524902, 11.490538835525513, 11.408236980438232, 11.429368019104004, 11.483153104782104, 11.532360076904297, 11.781499147415161, 11.747456073760986, 11.791171073913574]
	indelible2_1000_scale01=[11.074132919311523, 12.201006174087524, 12.127370834350586, 12.11881399154663, 11.880905866622925, 11.744970083236694, 11.71963906288147, 11.66642689704895, 11.640264987945557, 11.676317930221558]
	indelible2_1000_scale10=[12.04808497428894, 11.315619945526123, 11.326038122177124, 11.52104902267456, 11.497009992599487, 11.315742015838623, 11.406708002090454, 11.539706945419312, 11.63560700416565, 11.623785018920898]
	
	names=["phastSim (100,000 tips)","SeqGen (5,000 tips)","INDELible-m1 (1000 tips)","INDELible-m2 (1000 tips)"]
	colors=["red","purple","yellow","brown"]
	lengths=["0.1","1","10"]
	times=[[phastSim100000_scale01,seqgen5000_scale01,indelible1000_scale01,indelible2_1000_scale01],[phastSim100000,seqgen5000,indelible1000,indelible2_1000],[phastSim100000_scale10,seqgen5000_scale10,indelible1000_scale10,indelible2_1000_scale10]]
	boxplot(times,names,pathSimu+"boxplot_times_rescale.pdf",lengths,colors,'branch length rescale factor',70,degreeSkew=0)
	
	lengths = [0.1, 1, 10]
	#errplot(times,names,pathSimu+"errplot_times_rescale.pdf",lengths,colors,'branch length rescale factor',70)
	
	# indel simulation comparison - phastSim, indelible-v1, indelible-v2
	# number of tips: 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000
	names=["tree generation","phastSim","phastSim+Fasta","INDELible-m1","INDELible-m2"]
	colors=["blue","red","orange","yellow","brown"]

	times10=[[0.32546424865722656, 0.2624638080596924, 0.2651400566101074, 0.26073694229125977, 0.26463890075683594, 0.2659893035888672, 0.26677727699279785, 0.3185391426086426, 0.31917643547058105, 0.3141472339630127], [2.1556780338287354, 2.0433263778686523, 2.0150911808013916, 2.0070581436157227, 2.0596706867218018, 2.1369547843933105, 2.079869031906128, 2.3088269233703613, 2.3454384803771973, 2.3637101650238037], [2.025832176208496, 2.0558063983917236, 2.0284719467163086, 2.1088755130767822, 2.200312852859497, 2.0535054206848145, 2.2374017238616943, 2.2551021575927734, 2.2804694175720215, 2.2652318477630615], [0.04556894302368164, 0.05311298370361328, 0.04337048530578613, 0.0427241325378418, 0.04316425323486328, 0.04216146469116211, 0.04281044006347656, 0.04262089729309082, 0.0432891845703125, 0.054961442947387695], [0.02457427978515625, 0.029098987579345703, 0.030329465866088867, 0.03525996208190918, 0.031858205795288086, 0.03503155708312988, 0.02840280532836914, 0.036489248275756836, 0.03365349769592285, 0.030407190322875977 ]]
	times20=[[0.31945300102233887, 0.3333101272583008, 0.3152434825897217, 0.3256418704986572, 0.3227427005767822, 0.3118453025817871, 0.322324275970459, 0.3194284439086914, 0.32078099250793457, 0.3193016052246094], [2.3566293716430664, 2.358487844467163, 2.3402223587036133, 2.346405506134033, 2.355820655822754, 2.3659119606018066, 2.375427722930908, 2.359123706817627, 2.3722352981567383, 2.330568552017212], [2.273454427719116, 2.253715991973877, 2.2945330142974854, 2.2657666206359863, 2.286865234375, 2.2691919803619385, 2.2933640480041504, 2.2502806186676025, 2.2781882286071777, 2.3005595207214355], [0.07726049423217773, 0.07605361938476562, 0.07612442970275879, 0.0782313346862793, 0.07828712463378906, 0.07755136489868164, 0.09316277503967285, 0.07691383361816406, 0.07537221908569336, 0.09715390205383301], [0.06331181526184082, 0.06451702117919922, 0.05375504493713379, 0.0571591854095459, 0.08751606941223145, 0.0524594783782959, 0.05438733100891113, 0.054318904876708984, 0.055394649505615234, 0.06630206108093262 ]]
	times50=[[0.3166079521179199, 0.3257465362548828, 0.32345080375671387, 0.31409215927124023, 0.32008934020996094, 0.32721519470214844, 0.3242166042327881, 0.3243870735168457, 0.32608747482299805, 0.3179609775543213], [2.3272178173065186, 2.405362844467163, 2.3365321159362793, 2.3416860103607178, 2.3167293071746826, 2.3644657135009766, 2.357297658920288, 2.3380234241485596, 2.3674230575561523, 2.328831911087036], [2.298748254776001, 2.269660711288452, 2.2928011417388916, 2.292900323867798, 2.2788147926330566, 2.3371686935424805, 2.315415143966675, 2.282987594604492, 2.2783312797546387, 2.259902000427246], [0.17986273765563965, 0.18409228324890137, 0.18003153800964355, 0.18171358108520508, 0.18166804313659668, 0.18123126029968262, 0.18381762504577637, 0.1909928321838379, 0.20838594436645508, 0.18513798713684082], [0.15886974334716797, 0.15141773223876953, 0.24355506896972656, 0.15930700302124023, 0.14786005020141602, 0.14998626708984375, 0.14945697784423828, 0.13439321517944336, 0.15129804611206055, 0.14148926734924316 ]]
	times100=[[0.32665514945983887, 0.3335232734680176, 0.3235647678375244, 0.33245038986206055, 0.32474303245544434, 0.32438158988952637, 0.3254122734069824, 0.3211653232574463, 0.33133578300476074, 0.31354546546936035], [2.370936393737793, 2.390069007873535, 2.3693268299102783, 2.3474628925323486, 2.3586549758911133, 2.353419065475464, 2.3524601459503174, 2.360806703567505, 2.364651918411255, 2.3899550437927246], [2.295365571975708, 2.2972891330718994, 2.3280863761901855, 2.318568706512451, 2.2597038745880127, 2.2971644401550293, 2.4166100025177, 2.2900681495666504, 2.31244158744812, 2.299140214920044], [0.3617417812347412, 0.3728368282318115, 0.36397361755371094, 0.37748003005981445, 0.3655531406402588, 0.38352227210998535, 0.36750197410583496, 0.36978816986083984, 0.3962726593017578, 0.3859894275665283], [0.29948878288269043, 0.3197615146636963, 0.28223109245300293, 0.30565762519836426, 0.3340730667114258, 0.298659086227417, 0.32448458671569824, 0.2979152202606201, 0.32575511932373047, 0.3420422077178955 ]]
	times200=[[0.31435322761535645, 0.31729888916015625, 0.33643007278442383, 0.32003211975097656, 0.29729342460632324, 0.30275464057922363, 0.3011605739593506, 0.3002660274505615, 0.3044776916503906, 0.2939119338989258], [2.4505326747894287, 2.3537113666534424, 2.3443708419799805, 2.358020067214966, 2.3700480461120605, 2.3751583099365234, 2.329376459121704, 2.430290699005127, 2.3356897830963135, 2.3660731315612793], [2.3395607471466064, 2.3806312084198, 2.331080913543701, 2.370893716812134, 2.316906690597534, 2.4248690605163574, 2.3881564140319824, 2.372387409210205, 2.3621037006378174, 2.3346774578094482], [0.8451361656188965, 0.81907057762146, 0.7356603145599365, 0.8783695697784424, 0.7894411087036133, 0.755021333694458, 0.7861297130584717, 0.7856440544128418, 0.7581984996795654, 0.7601583003997803], [0.5781629085540771, 0.5203299522399902, 0.6084508895874023, 0.5540952682495117, 0.6588683128356934, 0.6407825946807861, 0.649280309677124, 0.6574001312255859, 0.8320307731628418, 0.6658432483673096 ]]
	times500=[[0.3241617679595947, 0.30945801734924316, 0.30818819999694824, 0.3087890148162842, 0.3054022789001465, 0.32735371589660645, 0.301347017288208, 0.3094305992126465, 0.3248178958892822, 0.3079488277435303], [2.4548282623291016, 2.3391318321228027, 2.297805070877075, 2.387371063232422, 2.2436463832855225, 2.649502754211426, 2.399074077606201, 2.421020746231079, 2.4090769290924072, 2.3196232318878174], [2.4686121940612793, 2.460338592529297, 2.374307632446289, 2.5748894214630127, 2.5234196186065674, 2.577294111251831, 2.446418046951294, 2.4995436668395996, 2.495103120803833, 2.5030055046081543], [2.233119487762451, 2.2040374279022217, 2.210883617401123, 2.304008960723877, 2.047483205795288, 2.13580060005188, 2.1504695415496826, 2.0436949729919434, 2.2270612716674805, 2.1424834728240967], [1.566122055053711, 1.5635943412780762, 1.4824836254119873, 1.3837003707885742, 1.404653787612915, 1.3873927593231201, 1.3557803630828857, 1.3557384014129639, 1.456052541732788, 1.6567285060882568 ]]
	times1000=[[0.323681116104126, 0.3277304172515869, 0.3154482841491699, 0.3198223114013672, 0.32904601097106934, 0.3184659481048584, 0.32576870918273926, 0.3365957736968994, 0.3242917060852051, 0.3187093734741211], [2.371051073074341, 2.254274845123291, 2.351447343826294, 2.8683762550354004, 2.44075083732605, 2.374760389328003, 2.265185832977295, 2.3820109367370605, 2.2783946990966797, 2.364931344985962], [2.881812334060669, 2.8860726356506348, 2.555190324783325, 2.6442477703094482, 2.692300319671631, 2.6801698207855225, 2.6771442890167236, 2.6435608863830566, 2.655925750732422, 2.6343815326690674], [4.412667512893677, 4.473804712295532, 4.408661603927612, 4.612804174423218, 4.409607172012329, 4.496532917022705, 4.276963710784912, 4.577553033828735, 4.600027322769165, 4.296308279037476], [2.999229669570923, 3.2141993045806885, 2.8491432666778564, 2.756218671798706, 3.016702175140381, 2.8072690963745117, 3.433241605758667, 3.1156091690063477, 2.9944984912872314, 3.1689796447753906 ]]
	times2000=[[0.3514270782470703, 0.35779261589050293, 0.35210466384887695, 0.3484079837799072, 0.36373448371887207, 0.36445188522338867, 0.3786911964416504, 0.3601212501525879, 0.35206127166748047, 0.36800146102905273], [2.3589372634887695, 3.038125514984131, 2.4358866214752197, 2.4174444675445557, 2.434441328048706, 2.346853256225586, 2.5356507301330566, 2.384561777114868, 2.4205384254455566, 2.4122507572174072], [2.861046314239502, 3.133699655532837, 3.1299421787261963, 3.402874708175659, 3.2692577838897705, 3.223698616027832, 3.3882076740264893, 2.953687906265259, 3.1378424167633057, 3.1192641258239746], [8.86726999282837, 9.359963417053223, 11.262237071990967, 9.429405212402344, 9.390536308288574, 8.65852427482605, 8.90850830078125, 9.446406126022339, 9.125349521636963, 9.307991981506348], [5.787947177886963, 5.444770812988281, 5.99629020690918, 6.386816024780273, 6.040544271469116, 5.958037853240967, 6.130545139312744, 5.8382158279418945, 5.549508571624756, 5.9274632930755615 ]]
	times5000=[[0.4435710906982422, 0.4505586624145508, 0.46629881858825684, 0.47099757194519043, 0.47298765182495117, 0.4531590938568115, 0.47820043563842773, 0.46857595443725586, 0.4473426342010498, 0.47362518310546875], [2.7430531978607178, 2.8311283588409424, 2.8181192874908447, 2.6468775272369385, 2.714951276779175, 2.8614823818206787, 2.9344606399536133, 2.8177590370178223, 2.663815975189209, 2.8534817695617676], [5.345435857772827, 5.129849433898926, 4.8792314529418945, 4.964365005493164, 4.670850992202759, 5.54222297668457, 6.245495796203613, 5.196593999862671, 4.57506251335144, 5.526819229125977], [22.235101461410522, 23.393762350082397, 23.297765731811523, 21.836959838867188, 21.698474168777466, 22.42943286895752, 23.31689691543579, 23.26853609085083, 22.301843404769897, 23.24411177635193], [15.466925144195557, 14.184302568435669, 14.439735412597656, 14.578090190887451, 14.50532341003418, 15.01891303062439, 15.582456111907959, 14.817336559295654, 14.58914566040039, 17.945835828781128 ]]
	times10000=[[0.6393427848815918, 0.8643717765808105, 0.8417501449584961, 0.8623208999633789, 0.8620443344116211, 0.7940115928649902, 0.7983183860778809, 0.8101792335510254, 0.820847749710083, 0.8228857517242432], [3.37762188911438, 3.5640625953674316, 3.488410472869873, 3.501112937927246, 3.525834321975708, 3.3079264163970947, 3.293438196182251, 3.414933681488037, 3.3546743392944336, 3.5214638710021973], [8.738857507705688, 8.462504386901855, 9.897857904434204, 8.927416324615479, 11.436707258224487, 9.176197290420532, 7.158588886260986, 6.7803122997283936, 8.938203811645508, 8.948585987091064], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [ 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ]]
	times20000=[[1.0342249870300293, 1.0626230239868164, 1.026402473449707, 1.0000288486480713, 1.035221815109253, 1.085937738418579, 1.071202039718628, 1.0741698741912842, 0.9978606700897217, 2.2949447631835938], [4.18130898475647, 4.177823543548584, 4.231494188308716, 4.272785186767578, 4.3068318367004395, 4.1911780834198, 4.394071102142334, 4.267273187637329, 4.4542577266693115, 7.634886980056763], [13.645295143127441, 12.847994089126587, 12.240689039230347, 11.498322248458862, 12.562551975250244, 14.182104110717773, 12.418807029724121, 13.753870248794556, 12.3284170627594, 19.558689832687378], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ]]
	times50000=[[2.3604001998901367, 2.474050521850586, 2.2536497116088867, 2.3335251808166504, 2.3143179416656494, 2.323009967803955, 2.289350748062134, 2.356456756591797, 2.2813494205474854, 4.916504144668579], [7.5116636753082275, 8.104324579238892, 7.911961555480957, 7.674963474273682, 7.6164631843566895, 8.099236011505127, 7.426865577697754, 7.491829872131348, 7.959242343902588, 12.753919839859009], [24.294604778289795, 23.12367296218872, 25.443199634552002, 27.09094762802124, 24.62966513633728, 25.876599073410034, 21.594374895095825, 22.818297147750854, 22.23173689842224, 39.72312140464783], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ]]
	times100000=[[5.205708980560303, 4.865880966186523, 4.8200905323028564, 4.796455144882202, 4.789062261581421, 4.758198499679565, 4.938420534133911, 5.058539390563965, 4.974266290664673, 9.818556308746338], [13.955179452896118, 13.301479816436768, 13.86931824684143, 13.275344133377075, 13.166597127914429, 13.78001856803894, 12.810855388641357, 13.321237087249756, 13.175382137298584, 23.704246044158936], [43.47269606590271, 44.65556740760803, 41.99473738670349, 40.41061878204346, 41.905001640319824, 40.23197841644287, 40.65919733047485, 42.2672119140625, 40.82607626914978, 76.90111446380615], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ]]
	times200000=[[10.31586766242981, 10.294050216674805, 11.406583070755005, 10.03695273399353, 10.388937950134277, 10.030796527862549, 9.85703420639038, 9.926616430282593, 10.26825475692749, 26.83154296875], [23.811527252197266, 25.46506667137146, 25.27086591720581, 23.648561716079712, 24.2943594455719, 23.874767303466797, 24.403047800064087, 25.116068601608276, 27.48458504676819, 60.54904747009277], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ]]
	times500000=[[26.597529888153076, 26.255186319351196, 25.45920157432556, 26.60333824157715, 25.42155122756958, 26.99942398071289, 27.04834246635437, 26.52974557876587, 25.9949049949646], [55.8808376789093, 55.209341526031494, 58.32216477394104, 55.4736225605011, 58.93256855010986, 55.726035594940186, 55.22979021072388, 59.64866805076599, 55.40555381774902],['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN' ]]
	nLeaves=["10","20","50","100","200","500","1000","2000","5000","10^4", "2x10^4", "5x10^4","10^5", "2x10^5", "5x10^5"]
	times=[times10,times20,times50,times100,times200,times500,times1000,times2000,times5000,times10000, times20000, times50000, times100000, times200000, times500000]
	boxplot(times, names, pathSimu+"boxplot_times_indels.pdf", nLeaves, colors, "", 40)
	nLeaves = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]
	errplot(times, names, pathSimu+"errplot_times_indels.pdf", nLeaves, colors, "Number of tips", 40)

exit()











