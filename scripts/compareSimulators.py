import sys
import os
#import math
import numpy as np
import os.path
#from os import path
import argparse
from ete3 import Tree
import time

#script to compare the running time of phastSim with pyvolve, seq-gen and indelible.

# python2 compareSimulators.py --nLeaves 100 --createFasta --pyvolveSim --seqgenSim --indelibleSim --indelibleSim2

parser = argparse.ArgumentParser(description='Compare simulators, in particular phastSim, seq-gen, indelible and pyvolve.')
parser.add_argument('--path',default="/Users/demaio/Desktop/coronavirus/simulations/", help='Path where to write simulation output.')
parser.add_argument('--reference',default="/Users/demaio/Documents/GitHub/phastSim/phastSim/example/MN908947.3.fasta", help='File containing the reference genome to be used as root genome.')
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
createFasta=args.createFasta
generatePlots=args.generatePlots
phastSim=args.phastSim

hierarchy=not args.noHierarchy

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
	fileRef=open(pathSimu+reference+"_new.fa","w")
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
	

times=[[],[],[],[],[],[],[],[]]
for r in range(nReplicates):
	times2=[]
	#simulate tree
	start = time.time()
	treeFile="simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+".tree"
	treeFile2="simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_2.tree"

	#use custom script - works well even for huge phylogenies
	import random_tree
	tree = random_tree.gen_tree(length, 0.0*float(length), min_leaves=nLeaves, labels="enum",seed=seed+r)
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
	start = time.time()
	
	stringRun="phastSimulate --outpath "+pathSimu+"  --seed "+str(seed+r)+" --scale "+str(scale)+" --treeFile "+pathSimu+treeFile+" --reference "+reference+" --outputFile simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+".txt"
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
	if phastSim:
		os.system(stringRun+" >/dev/null")
	
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
		
	if not hierarchy:
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
		
		stringRun="python2 /Users/demaio/Desktop/coronavirus/simulations/scripts/runPyvolve.py --path /Users/demaio/Desktop/coronavirus/simulations/  --seed "+str(seed+r)+" --scale "+str(scale)+" --treeFile "+treeFile2+" --outputFile simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_pyvolve.txt"
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
		stringRun="/Applications/Seq-Gen-1.3.4/source/seq-gen -l "+str(length)+" -s "+str(scale)+" -m GTR -f 0.3 0.2 0.2 0.3 -r "+str(mutMatrix[0][1]/mutMatrix[2][3])+" "+str(mutMatrix[0][2]/mutMatrix[2][3])+" "+str(mutMatrix[0][3]/mutMatrix[2][3])+" "+str(mutMatrix[1][2]/mutMatrix[2][3])+" "+str(mutMatrix[1][3]/mutMatrix[2][3])+" 1.0 "+categoryString+" -z "+str(seed+r)+" "+pathSimu+treeFile+" > "+pathSimu+"simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_seqgenOutput.txt"
		os.system(stringRun+" &>/dev/null")
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
	from matplotlib.patches import Polygon
	
	names=["tree generation","phastSim","phastSim+Fasta","pyvolve","SeqGen","INDELible-m1","INDELible-m2"]
	colors=["blue","red","orange","green","purple","yellow","brown"]
	times10=[[0.0013659000396728516, 0.0012090206146240234, 0.0013148784637451172, 0.0010118484497070312, 0.0013380050659179688, 0.002196073532104492, 0.0014619827270507812, 0.001085042953491211, 0.001199960708618164, 0.0013890266418457031], [1.3747148513793945, 1.3423030376434326, 1.3215079307556152, 1.3170280456542969, 1.3145039081573486, 1.3163139820098877, 1.3323869705200195, 1.3344120979309082, 1.399022102355957, 1.3535301685333252], [1.3641538619995117, 1.3581798076629639, 1.334787130355835, 1.3404159545898438, 1.3326139450073242, 1.3254010677337646, 1.3380720615386963, 1.3444411754608154, 1.3227300643920898, 1.3296558856964111], [10.451174020767212, 10.58082389831543, 10.541486978530884, 10.444769144058228, 10.375016927719116, 10.394562005996704, 10.366996049880981, 10.622323989868164, 10.495352983474731, 10.429331064224243], [0.03925800323486328, 0.032778024673461914, 0.033091068267822266, 0.03274106979370117, 0.036557912826538086, 0.032247066497802734, 0.03708791732788086, 0.041892051696777344, 0.038690805435180664, 0.035842180252075195], [0.21403002738952637, 0.21619105339050293, 0.21575593948364258, 0.21455693244934082, 0.2209010124206543, 0.221451997756958, 0.21570682525634766, 0.2141880989074707, 0.22292304039001465, 0.21918296813964844], [0.12611007690429688, 0.1258530616760254, 0.1236422061920166, 0.12346792221069336, 0.12584614753723145, 0.12546515464782715, 0.12303709983825684, 0.12221598625183105, 0.1220400333404541, 0.12607908248901367]]
	times20=[[0.0024929046630859375, 0.0014770030975341797, 0.0023810863494873047, 0.0014848709106445312, 0.001817941665649414, 0.0016529560089111328, 0.0020999908447265625, 0.0015180110931396484, 0.0020749568939208984, 0.0017440319061279297], [1.645247220993042, 1.310811996459961, 1.3231310844421387, 1.3270468711853027, 1.3239338397979736, 1.3379931449890137, 1.3164989948272705, 1.3218259811401367, 1.3193409442901611, 1.5228931903839111], [1.679689884185791, 1.3187479972839355, 1.3616609573364258, 1.329289197921753, 1.346217155456543, 1.3264689445495605, 1.3758349418640137, 1.3410701751708984, 1.3292009830474854, 1.6160540580749512], [22.576815128326416, 22.07684302330017, 22.30865716934204, 21.877304792404175, 21.605812072753906, 21.84848117828369, 22.243227005004883, 22.178505897521973, 22.18565797805786, 23.35649800300598], [0.06273508071899414, 0.05921602249145508, 0.0638580322265625, 0.05792999267578125, 0.05753588676452637, 0.05878496170043945, 0.057292938232421875, 0.06512904167175293, 0.05781197547912598, 0.057282209396362305], [0.44101905822753906, 0.44517993927001953, 0.43025708198547363, 0.4448831081390381, 0.44520998001098633, 0.44223713874816895, 0.43360185623168945, 0.4268670082092285, 0.46094799041748047, 0.4530940055847168], [0.24215412139892578, 0.24007606506347656, 0.24377012252807617, 0.23820710182189941, 0.24094796180725098, 0.2389540672302246, 0.24043488502502441, 0.23656988143920898, 0.29044508934020996, 0.25960588455200195]]
	times50=[[0.003216981887817383, 0.0038487911224365234, 0.0030329227447509766, 0.003144979476928711, 0.004025936126708984, 0.002969980239868164, 0.0034058094024658203, 0.002958059310913086, 0.023120880126953125, 0.02884507179260254], [1.2893130779266357, 1.6244549751281738, 1.486210823059082, 1.4786460399627686, 1.4211831092834473, 1.811263084411621, 1.376521110534668, 1.4128239154815674, 1.5160620212554932, 1.3880500793457031], [1.3925518989562988, 1.4420089721679688, 1.4690890312194824, 1.5120949745178223, 1.4923410415649414, 1.4661610126495361, 1.4662230014801025, 1.4426720142364502, 1.4030230045318604, 1.4014220237731934], [61.306437969207764, 63.28899002075195, 60.32176899909973, 61.62099599838257, 62.32091784477234, 60.96144986152649, 60.207597970962524, 59.269126892089844, 60.00672507286072, 59.49119710922241], [0.13117289543151855, 0.13464593887329102, 0.13071918487548828, 0.1415119171142578, 0.13554096221923828, 0.13779211044311523, 0.1376659870147705, 0.13895797729492188, 0.13930082321166992, 0.1320171356201172], [1.1521000862121582, 1.117422103881836, 1.1462361812591553, 1.1094839572906494, 1.1246919631958008, 1.1122980117797852, 1.109248161315918, 1.1435189247131348, 1.1563069820404053, 1.1218719482421875], [0.6041860580444336, 0.6159520149230957, 0.6443181037902832, 0.6238491535186768, 0.5989768505096436, 0.5900440216064453, 0.609623908996582, 0.5931220054626465, 0.6229169368743896, 0.600348949432373]]
	times100=[[0.009551048278808594, 0.0053958892822265625, 0.005692958831787109, 0.00590205192565918, 0.005836963653564453, 0.0058939456939697266, 0.005366086959838867, 0.005799770355224609, 0.006151914596557617, 0.0052869319915771484], [1.90081787109375, 1.3287749290466309, 1.3533220291137695, 1.4254589080810547, 1.5417671203613281, 1.4508280754089355, 1.4114508628845215, 1.549408197402954, 1.5193719863891602, 1.400057077407837], [1.3768410682678223, 1.385685920715332, 1.412989854812622, 1.3838729858398438, 1.5289270877838135, 1.4544739723205566, 1.4677560329437256, 1.51531982421875, 1.4512879848480225, 1.4553320407867432], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [0.2564408779144287, 0.2659909725189209, 0.26821398735046387, 0.2528820037841797, 0.26355600357055664, 0.2617809772491455, 0.25960516929626465, 0.25537586212158203, 0.24987196922302246, 0.25585007667541504], [2.20414400100708, 2.1862261295318604, 2.231003999710083, 2.1976799964904785, 2.2304821014404297, 2.2315008640289307, 2.573899984359741, 2.246011972427368, 2.2906219959259033, 2.20377516746521], [1.1515510082244873, 1.1431560516357422, 1.1758790016174316, 1.146554946899414, 1.1703548431396484, 1.178049087524414, 1.3542709350585938, 1.2229149341583252, 1.154651165008545, 1.1641769409179688]]
	times200=[[0.016273975372314453, 0.011254072189331055, 0.011461973190307617, 0.010603904724121094, 0.013769149780273438, 0.013958930969238281, 0.026243925094604492, 0.014862060546875, 0.011723995208740234, 0.010759115219116211], [1.9284470081329346, 1.7488529682159424, 1.5611369609832764, 1.5293529033660889, 1.6704580783843994, 1.5298941135406494, 1.548914909362793, 1.8947179317474365, 1.6275060176849365, 1.5622169971466064], [1.761125087738037, 1.7953181266784668, 1.5314640998840332, 1.6465811729431152, 1.7108330726623535, 1.6326282024383545, 1.6300709247589111, 1.7247998714447021, 1.568471908569336, 2.247066020965576], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [0.4969940185546875, 0.5056090354919434, 0.4887361526489258, 0.5002598762512207, 0.49825382232666016, 0.49005913734436035, 0.4913029670715332, 0.5112559795379639, 0.49417805671691895, 0.4985392093658447], [5.168694019317627, 5.151010990142822, 5.2190330028533936, 5.296245098114014, 5.013715982437134, 5.4771199226379395, 6.150067090988159, 5.349929094314575, 4.977146148681641, 5.350111961364746], [2.8624649047851562, 2.833103895187378, 2.7423510551452637, 2.6613221168518066, 2.663487195968628, 2.860327959060669, 3.3309261798858643, 2.802933931350708, 2.646231174468994, 2.517848014831543]]
	times500=[[0.030160903930664062, 0.02469801902770996, 0.04264998435974121, 0.03104686737060547, 0.03468894958496094, 0.02983403205871582, 0.02934098243713379, 0.024353981018066406, 0.02561211585998535, 0.03130292892456055], [1.777817964553833, 1.9888050556182861, 1.9227111339569092, 1.9003281593322754, 2.0147860050201416, 1.9665369987487793, 1.8426148891448975, 1.4138050079345703, 1.4871950149536133, 1.5098040103912354], [1.5789639949798584, 1.5292949676513672, 1.564305067062378, 1.691133975982666, 1.7239899635314941, 1.560417890548706, 1.5549790859222412, 1.6500108242034912, 1.5650379657745361, 1.653033971786499], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [1.222322940826416, 1.2483329772949219, 1.2385191917419434, 1.2624180316925049, 1.2348461151123047, 1.3084070682525635, 1.2472999095916748, 1.2439031600952148, 1.2595350742340088, 1.2943921089172363], [11.081032991409302, 11.013811111450195, 11.17693018913269, 11.416763067245483, 11.85645604133606, 11.123574018478394, 11.086908102035522, 11.227660894393921, 11.142872095108032, 11.275111198425293], [5.7495949268341064, 5.812233924865723, 6.02468204498291, 5.983781099319458, 6.056072950363159, 5.8660888671875, 5.7771899700164795, 5.78557014465332, 5.831812143325806, 6.089737892150879]]
	times1000=[[0.06780815124511719, 0.06437397003173828, 0.05380582809448242, 0.04997706413269043, 0.05129718780517578, 0.04825186729431152, 0.054646968841552734, 0.04891395568847656, 0.05314779281616211, 0.049278974533081055], [1.4967939853668213, 1.9933040142059326, 1.5371029376983643, 1.5284991264343262, 1.9865779876708984, 1.506382942199707, 1.9183871746063232, 2.0416319370269775, 1.5037190914154053, 1.8662559986114502], [1.849992036819458, 1.8258349895477295, 1.709947109222412, 1.6919009685516357, 1.8142428398132324, 1.6936531066894531, 1.898118019104004, 1.7096149921417236, 1.7108089923858643, 1.7036280632019043], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [2.4352219104766846, 2.4571030139923096, 2.4070231914520264, 2.4291841983795166, 2.451385974884033, 2.439326047897339, 2.4405691623687744, 2.484835147857666, 2.471928834915161, 2.465049982070923], [22.25348210334778, 22.035884857177734, 21.825459957122803, 21.92672896385193, 21.81256103515625, 22.11038112640381, 21.884013175964355, 21.707130908966064, 21.84363889694214, 21.670482873916626], [11.619686126708984, 11.556863069534302, 11.478086948394775, 11.657796144485474, 11.494728088378906, 11.577374935150146, 11.585293054580688, 11.68454885482788, 11.430992126464844, 11.45626187324524]]
	times2000=[[0.12939810752868652, 0.099761962890625, 0.10983014106750488, 0.137160062789917, 0.11310815811157227, 0.12170195579528809, 0.17563486099243164, 0.10024499893188477, 0.0971071720123291, 0.10719680786132812], [1.6819391250610352, 2.0618398189544678, 2.0402579307556152, 2.0742900371551514, 1.892535924911499, 2.2402780055999756, 1.8665268421173096, 1.5126521587371826, 1.5043201446533203, 1.7153770923614502], [2.1061370372772217, 2.091276168823242, 2.240056037902832, 2.9656341075897217, 2.3872971534729004, 2.5221140384674072, 2.6080589294433594, 2.0910890102386475, 2.0963239669799805, 2.1601290702819824], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [4.886102914810181, 4.879283905029297, 4.931666135787964, 4.81746506690979, 4.882803916931152, 4.861746072769165, 4.885030031204224, 4.975928068161011, 4.911868095397949, 4.934772968292236], [44.09892296791077, 44.005027055740356, 45.88066387176514, 48.627665996551514, 56.57271981239319, 52.096412897109985, 44.491461992263794, 43.75139808654785, 44.21598196029663, 44.35627007484436], [23.274376153945923, 23.85594391822815, 25.263535976409912, 25.600807905197144, 26.603089809417725, 27.35308003425598, 23.306721925735474, 23.075037002563477, 23.742963075637817, 23.250370979309082]]
	times5000=[[0.31500887870788574, 0.24887490272521973, 0.3096740245819092, 0.2508058547973633, 0.3484029769897461, 0.2693758010864258, 0.27234315872192383, 0.35924696922302246, 0.3092200756072998, 0.42699503898620605], [1.9736239910125732, 2.429014205932617, 1.9887659549713135, 2.0318288803100586, 2.406982183456421, 2.0246009826660156, 1.9706370830535889, 2.4449448585510254, 2.524501085281372, 2.53228497505188], [3.344974994659424, 3.331270933151245, 3.323568820953369, 3.300321102142334, 3.405689001083374, 3.3196661472320557, 3.3752760887145996, 3.3163578510284424, 3.508978843688965, 3.707275867462158], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [12.1112220287323, 12.05104112625122, 12.257277965545654, 12.403314113616943, 12.445502042770386, 12.404994010925293, 12.645019054412842, 12.529222965240479, 12.521138906478882, 12.774863004684448], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [59.62119007110596, 58.707270860672, 58.55328392982483, 59.17231893539429, 63.55593514442444, 59.023497104644775, 58.952008962631226, 59.784512996673584, 63.09826993942261, 60.914510011672974]]
	times10000=[[0.5868499279022217, 0.580380916595459, 0.6529080867767334, 0.6462969779968262, 0.6883280277252197, 0.6795670986175537, 0.6926059722900391, 0.4921681880950928, 0.7163400650024414, 0.6498241424560547], [2.942147970199585, 2.4787280559539795, 2.444871187210083, 2.4249610900878906, 2.5257160663604736, 2.4837939739227295, 2.465178966522217, 2.3728160858154297, 2.475499153137207, 2.415358066558838], [5.290570974349976, 5.2421088218688965, 5.169349908828735, 5.3277130126953125, 5.383677005767822, 5.364545106887817, 5.328459978103638, 5.3075220584869385, 5.364189147949219, 5.305814981460571], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [23.960747957229614, 23.89430594444275, 24.102982997894287, 24.04936909675598, 24.217010021209717, 24.367241144180298, 24.017180919647217, 24.123044967651367, 24.17894697189331, 23.99029302597046], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	times20000=[[1.1991140842437744, 1.327165126800537, 1.4817399978637695, 1.4684429168701172, 1.7353639602661133, 1.6934490203857422, 1.5031070709228516, 1.3446290493011475, 1.4934370517730713, 1.8493130207061768], [3.3477210998535156, 3.769974946975708, 3.3797521591186523, 3.4234910011291504, 3.576843023300171, 3.5666489601135254, 3.539125919342041, 3.4996449947357178, 3.4938271045684814, 3.4616119861602783], [9.40938401222229, 9.229989051818848, 9.341056823730469, 9.800310134887695, 9.614714860916138, 9.365905046463013, 9.440688848495483, 9.506255149841309, 9.528976917266846, 9.752545833587646], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], [48.57287001609802, 48.325536012649536, 48.499929904937744, 48.55190300941467, 48.723459005355835, 48.555269956588745, 48.50492286682129, 48.813239097595215, 48.42090392112732, 48.982256174087524], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	times50000=[[3.293437957763672, 4.06860089302063, 3.951677083969116, 4.423721075057983, 3.8627521991729736, 4.419393062591553, 4.29544997215271, 4.469686031341553, 4.260180950164795, 4.0369179248809814], [6.553931951522827, 6.685858964920044, 6.578627824783325, 6.476964950561523, 6.540050029754639, 6.447605133056641, 6.578918933868408, 7.231057167053223, 6.6752049922943115, 6.765127897262573], [21.58997893333435, 20.89152193069458, 21.59639310836792, 21.675976037979126, 21.34534788131714, 22.414045095443726, 25.244637966156006, 22.464181900024414, 21.908604860305786, 22.150094032287598], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	times100000=[[7.647463083267212, 8.008918046951294, 9.324841022491455, 9.470612049102783, 8.557217121124268, 9.787823915481567, 9.107623100280762, 9.38093113899231, 8.35722804069519, 8.31895899772644], [11.86309003829956, 11.798933982849121, 12.203128099441528, 11.762664079666138, 12.412771940231323, 12.195564985275269, 11.996767044067383, 11.904637813568115, 11.997853994369507, 12.214859962463379], [42.35268306732178, 40.621132135391235, 42.493098974227905, 42.65987491607666, 44.19890284538269, 42.4167218208313, 43.31222414970398, 43.10220289230347, 43.47725200653076, 43.2533118724823], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	times200000=[[15.817861080169678, 16.58775305747986, 19.115757942199707, 16.761270999908447, 19.467989206314087, 17.577764987945557, 17.07235097885132, 19.266084909439087, 16.83725094795227, 18.492412090301514], [21.678630113601685, 21.98533320426941, 22.04712176322937, 22.029582977294922, 22.8069269657135, 21.75388789176941, 21.656525135040283, 21.720540046691895, 21.85643482208252, 21.551968097686768], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	times500000=[[41.40035009384155, 43.9114990234375, 55.79817795753479, 47.97320508956909, 54.219762086868286, 53.92414689064026, 46.42072582244873, 51.166882038116455, 45.05335998535156, 51.811765909194946], [54.970056772232056, 55.55358910560608, 56.28190588951111, 56.577943086624146, 59.569891929626465, 59.95699906349182, 65.13326716423035, 53.935949087142944, 52.92773199081421, 55.46649694442749], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'], ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN']]
	#times1000000=[]

	#nLeaves=[10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000]
	nLeaves=["10","20","50","100","200","500","1000","2000","5000","10^4","2x10^4","5x10^4","10^5","2x10^5","5x10^5"]
	times=[times10,times20,times50,times100,times200,times500,times1000,times2000,times5000,times10000,times20000,times50000,times100000,times200000,times500000]
	
	
	#BOXPlot drawing
	def boxplot(valuesLists,axisLabels,plotFileName,labels,colors, xLabel,topPlot):
	
		fig, ax1 = plt.subplots(figsize=(10, 6))
		fig.canvas.set_window_title('Simulation Running times')
		fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
		num_boxes = len(valuesLists)*len(valuesLists[0])
		
		data=np.zeros((len(valuesLists)*len(valuesLists[0]),len(valuesLists[0][0])))
		print((len(valuesLists)*len(valuesLists[0]),len(valuesLists[0][0])))
		positions=[]
		for i in range(len(valuesLists)):
			for j in range(len(valuesLists[i])):
				c=colors[j]
				positions.append(i+(j-len(valuesLists[i])/2.0)*0.5/len(valuesLists[i]))
				for k in range(len(valuesLists[i][j])):
					data[i*len(valuesLists[0])+j][k]=valuesLists[i][j][k]
				position=[i+(j-(len(valuesLists[i])-1)/2.0)*0.5/len(valuesLists[i])]
				#if len(valuesLists[i])%2==0:
				#	position=[i+(j-(len(valuesLists[i])-1)/2.0)*0.5/len(valuesLists[i])]
				#else:
				#	position=[i+(j-(len(valuesLists[i])-1)/2.0)*0.5/len(valuesLists[i])]
				plt.boxplot(data[i*len(valuesLists[0])+j], positions=position, notch=False, patch_artist=True, widths=0.5/len(colors), manage_xticks=False, 
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
		ax1.set_xticklabels(labels, rotation=45, fontsize=11)

		# Finally, add a basic legend
		for i in range(len(colors)):
			fig.text(0.10, 0.85-0.045*i, axisLabels[i], backgroundcolor=colors[i], color='black', weight='roman', size='medium')
		
		fig.savefig(plotFileName)
		plt.close()

	#generate boxplot of general running times
	boxplot(times,names,pathSimu+"boxplot_times_general.pdf",nLeaves,colors,'Number of tips',70)




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




exit()











