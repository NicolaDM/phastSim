import sys
import os
import math
import numpy as np
import os.path
from os import path
import argparse
from ete3 import Tree
import time

parser = argparse.ArgumentParser(description='Run simulations about the effective alignment size (EASi).')
parser.add_argument('--path',default="", help='path where to run simulations.')
parser.add_argument('--reference',default="MN908947.3.fasta", help='file containing the reference genome to be used as root genome. To be found in the folder specified with --path.')
parser.add_argument('--treeFile',default="exampleTree.tree", help='Name of file containing the tree used to simulate sequences (assumed within the --path and in newick format).')
parser.add_argument('--scale',default=1.0,type=float, help='scale the simulation tree by this amount (default 1.0). Branch lengths are assumed in terms of expected substitutions per site (more or less, as frequencies changes through time, total mutation rate might also  change).')
parser.add_argument("--seed", help="seed for random simulator in genes simulations", type=int, default=1)
parser.add_argument("--mutationRates", help="mutation rates, by default using the neutral rates estimated from SARS-CoV-2; so far only 12 input values allowed (r_AC, r_AG, r_AT, r_CA, etc...).", type=float, nargs='+',  default=[0.039, 0.310, 0.123, 0.140, 0.022, 3.028, 0.747, 0.113, 2.953, 0.056, 0.261, 0.036])
parser.add_argument("--verbose", help="Turns on verbose mode.", action="store_true")
parser.add_argument('--outputFile',default="sars-cov-2_simulation_output.txt", help='Output file name containing the simulated genomes in succint format. The file will be created within the folder specified with --path.')
parser.add_argument("--createNewick", help="Create a newick file annotated with the simulated mutation events (default name sars-cov-2_simulation_output.tree).", action="store_true")
parser.add_argument('--newickO',default="sars-cov-2_simulation_output.tree", help='Output file name containing the tree annotated with simulated mutation events. The file will be created within the folder specified with --path.')
parser.add_argument("--createFasta", help="Create a fasta file with the simulated genomes (default name sars-cov-2_simulation_output.fasta).", action="store_true")
parser.add_argument('--fastaO',default="sars-cov-2_simulation_output.fasta", help='Output file name containing the simulated genomes. The file will be created within the folder specified with --path.')
parser.add_argument("--createPhylip", help="Create a phylip file with the simulated genomes (default name sars-cov-2_simulation_output.phy).", action="store_true")
parser.add_argument('--phylipO',default="sars-cov-2_simulation_output.phy", help='Output phylip file name containing the simulated genomes. The file will be created within the folder specified with --path.')
parser.add_argument("--pyvolveSim", help="run simulations using pyvolve (not the default, it's only to allow easy comparison)", action="store_true")
parser.add_argument('--pyvolveFastaO',default="pyvolve_simulation_output.fasta", help='Output file name containing the genomes simulated by pyvolve, if pyvolve is used. The file will be created within the folder specified with --path.')
parser.add_argument("--treeFormat", help="remove branch support values from tree file to use in pyvolve", action="store_true")
args = parser.parse_args()

pathSimu=args.path
reference=args.reference
np.random.seed(args.seed)
seed=args.seed
mutationRates=args.mutationRates
verbose=args.verbose
outputFile=args.outputFile
createNewick=args.createNewick
newickO=args.newickO
createFasta=args.createFasta
fastaO=args.fastaO
createPhylip=args.createPhylip
phylipO=args.phylipO
pyvolveSim=args.pyvolveSim
pyvolveFastaO=args.pyvolveFastaO
treeFormat=args.treeFormat
treeFile=args.treeFile


#possible alleles
alleles={"A":0,"C":1,"G":2,"T":3}
allelesList=["A","C","G","T"]

#collect reference
file=open(pathSimu+reference)
line=file.readline()
ref=""
while line!="":
	line=file.readline()
	ref+=line.replace("\n","")
file.close()
print("\n Finished reading reference genome at "+pathSimu+reference+" with "+str(len(ref))+" bases.")
refList=list(ref)

#SARS-CoV-2 genome annotation - not used yet but will be useful when simulating under a codon model.
geneEnds=[[266,13468],[13468,21555],[21563,25384],[25393,26220],[26245,26472],[26523,27191],[27202,27387],[27394,27759],[27894,28259],[28274,29533],[29558,29674]]

#substitution rates
if len(mutationRates)==12:
	print("\n Assuming UNREST nucleotide substitution model.")
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
	print("\n Number of mutation rates "+str(len(mutationRates))+", model not implemented yet.")
	print(mutationRates)
	exit()
print("\n Mutation rate matrix:")
print(mutMatrix)

#List of dictionaries that let you know at wich position (starting from 0) is the 1st A of the genome, the second A of the genome, etc (well, actually the 0th, the 1st, etc..).
#here I could use arrays instead of dictionaries, it might be more efficient.
positions=[{},{},{},{}]
#tot number of A's, C's, G's and T's in the reference genome
totAlleles=[0,0,0,0]
for pos in range(len(ref)):
	a=alleles[ref[pos]]
	positions[a][totAlleles[a]]=pos
	totAlleles[a]+=1
print("\n Number of each nucleotide in the genome:")
print(totAlleles)
		
#print information about a tree node
def printNode(node):
	print("\n Node name " + node.name+" , sampled="+str(node.is_leaf()))
	print(str(len(node.children))+ " children, distance from parent "+str(node.dist))
	print("Mutations:")
	for i in range(4):
		print(node.mutations[i])
			
#generate the genome sequence of a sample using node mutations and the reference
#useful, for example, for generating a fasta file.
def genomeSeq(mutations):
	seqList=list(refList)
	for i in range(4):
		for m in mutations[i]:
			seqList[positions[i][m[0]]]=allelesList[m[1]]
	return ''.join(seqList)

	
#Total mutation rates cumulatively across the genome
totMutMatrix=np.zeros((4,4),dtype=float)
totMut=0.0
for i in range(4):
	for j in range(4):
		totMutMatrix[i][j]=mutMatrix[i][j]*totAlleles[i]
	totMut=totMut-totMutMatrix[i][i]
print("\n Total genome-wide cumulative mutation rate: "+str(totMut))
if verbose:
	print("Base-wise:")
	print(totMutMatrix)
	


#Function to simulate evolution on one branch,using ETE tree structure;
# given details of the parent node, it generates details of the child node.
#To simulate whole tree, it needs to be called on root.
def mutateBranchETE(childNode,parentMuts,parentTotAlleles,parentRate):
	bLen=childNode.dist
	currTime=0.0
	#if newick output is requested, prepare format
	if createNewick:
		childNode.mutAnnotation=[]
	#Initialize child rate and allele numbers with parent ones
	rate=parentRate
	childTotAlleles=list(parentTotAlleles)
	childMutations=[[],[],[],[]]
	#Initialize child mutation list with parent one
	for i in range(4):
		for k in range(len(parentMuts[i])):
			childMutations[i].append(list(parentMuts[i][k]))
	
	#Sample new mutation event with Gillespie algorithm
	currTime+=np.random.exponential(scale=1.0/rate)
	if verbose:
		print("\n Node "+childNode.name+" BLen: "+str(bLen)+". Mutations at start:")
		print(parentMuts)
		print("num alleles at start:")
		print(parentTotAlleles)
		print("tot rate at start:"+str(parentRate))
		print("First sampled time: "+str(currTime)+" ; mutation rate: "+str(rate))
	while currTime<bLen:
		#Now, sample which type of mutation event it is (from which nucleotide to which nucleotide)
		rand=np.random.random()*rate
		tot=0.0
		found=False
		#print("Random value: "+str(rand))
		for i in range(4):
			for j in range(4):
				if j!=i:
					tot+=childTotAlleles[i]*mutMatrix[i][j]
					#print("i "+str(i)+" j "+str(j)+" tot "+str(tot))
					if rand<(tot):
						found=True
						break
			if found:
				break
		if not found:
			print("Error in selecting mutation type")
			exit()
		#Now, sample the specific position of the genome (among those with the mutated alleles) that mutates.
		mutatedBasePos=np.random.randint(childTotAlleles[i])
		if verbose:
			print("mutation position "+str(mutatedBasePos)+" allele "+str(i)+" to "+str(j))
		
		#So far I am not keeping track of the mutation events, but only of the ancestral and descendant genomes.
		#If the set of simulated mutation events is needed, I might need a structure to track them and write them to file in some format.
		
		#in this case, we are mutating a position that was not mutated before. We need to add one entry to the mutation list without removing any old one.
		if mutatedBasePos<totAlleles[i]-len(childMutations[i]):
			#print("New mutation")
			newMutPos=mutatedBasePos
			m=0
			while m<len(childMutations[i]):
			#for m in range(len(childMutations[i])):
				if childMutations[i][m][0]<=newMutPos:
					newMutPos+=1
					m+=1
				else:
					break
			#print(m)
			childMutations[i].insert(m,[newMutPos,j])
			if createNewick:
				childNode.mutAnnotation.append(allelesList[i]+str(positions[i][newMutPos]+1)+allelesList[j])
		#in this case, we are mutating a position that was already mutated.
		#this means that one item from the mutation list needs to be removed or modified, and no item needs to be added.
		else:
			newMutatedBasePos=mutatedBasePos-(totAlleles[i]-len(childMutations[i]))
			if verbose:
				print("Modifying pre-existing mutation")
				print(newMutatedBasePos)
			#exit()
			added=False
			for i2 in range(4):
				if i2!=i:
					for m in range(len(childMutations[i2])):
						if childMutations[i2][m][1]==i:
							if newMutatedBasePos==0:
								if createNewick:
									childNode.mutAnnotation.append(allelesList[i]+str(positions[i][childMutations[i2][m][0]]+1)+allelesList[j])
								if j==i2:
									if verbose:
										print("Deleted mutation \n\n\n")
									del childMutations[i2][m]
								else:
									childMutations[i2][m][1]=j
								added=True
								break
							newMutatedBasePos-=1
				if added:
					break
		
		childTotAlleles[i]-=1
		childTotAlleles[j]+=1
		rate+=mutMatrix[i][i]
		rate-=mutMatrix[j][j]
				
		currTime+=np.random.exponential(scale=1.0/rate)
		if verbose:
			print("new time "+str(currTime)+" mutation events:")
			print(childMutations)
		
	if verbose:
		print("mutations at the end:")
		print(childMutations)
	childNode.mutations=childMutations
	#now mutate children of the current node, calling this function recursively on the node children.
	for c in childNode.children:
		mutateBranchETE(c,childMutations,childTotAlleles,rate)





start = time.time()

# Loads a tree structure from a newick string in ETE2. The returned variable t is the root node for the tree.
t = Tree(pathSimu+treeFile)
scale=args.scale
for node in t.traverse("postorder"):
	node.dist=node.dist * scale
time1 = time.time() - start
print("Time for reading tree with ETE3: "+str(time1))

#Run sequence evolution simulation along tree
mutateBranchETE(t,[[],[],[],[]],totAlleles,totMut)
time2 = time.time() - start
print("Total time after simulating sequence evolution along tree with Gillespie approach: "+str(time2))






#function to write a succint output iteratively
def writeGenomeShort(node,file):
	if node.is_leaf():
		file.write(">"+node.name+"\n")
		mutDict={}
		for i in range(4):
			for m in node.mutations[i]:
				mutDict[positions[i][m[0]]+1]=allelesList[m[1]]
		mutList=list(mutDict.keys())
		mutList.sort()
		for m in mutList:
			file.write(str(m)+" "+mutDict[m]+"\n")
	for c in node.children:
		writeGenomeShort(c,file)

#Create a succint output
file=open(pathSimu+outputFile,"w")
writeGenomeShort(t,file)
file.close()
time3 = time.time() - start
print("Total time after writing short file: "+str(time3))

#function to write a newick output iteratively
def writeGenomeNewick(node):
	if node.is_leaf():
		outString=node.name+'['
	else:
		outString='('
		for c in range(len(node.children)):
			outString+=writeGenomeNewick(node.children[c])
			if c<len(node.children)-1:
				outString+=','
		outString+=')['
	stringToAdd=''
	for i in range(len(node.mutAnnotation)):
		stringToAdd+=node.mutAnnotation[i]
		if i<len(node.mutAnnotation)-1:
			stringToAdd+=','
	stringToAdd+=(']:'+str(node.dist))
	return outString+stringToAdd
	
#If requested, create a newick output
if createNewick:
	file=open(pathSimu+newickO,"w")
	newickTree=writeGenomeNewick(t)+";\n"
	file.write(newickTree)
	file.close()
	
	time3 = time.time() - start
	print("Total time after writing newick file: "+str(time3))

#function to write a fasta output iteratively
def writeGenome(node,file):
	if node.is_leaf():
		seq=genomeSeq(node.mutations)
		file.write(">"+node.name+"\n"+seq+"\n")
	for c in node.children:
		writeGenome(c,file)

#If requested, create a fasta output
if createFasta:
	file=open(pathSimu+fastaO,"w")
	writeGenome(t,file)
	file.close()
	
	time3 = time.time() - start
	print("Total time after writing fasta file: "+str(time3))
	
#function to write a phylip output iteratively
def writeGenomePhylip(node,file):
	if node.is_leaf():
		seq=genomeSeq(node.mutations)
		file.write(node.name+"\t"+seq+"\n")
	for c in node.children:
		writeGenomePhylip(c,file)

#If requested, create a phylip output
if createPhylip:
	file=open(pathSimu+phylipO,"w")
	file.write("\t"+str(len(t))+"\t"+str(len(ref))+"\n")
	writeGenomePhylip(t,file)
	file.close()
	
	time3 = time.time() - start
	print("Total time after writing phylip file: "+str(time3))

elapsedTime = time.time() - start
print("Overall time for new simulation approach: "+str(elapsedTime))








#If requested, modify tree format to remove branch supports - pyvolve can be picky wrt newick format
if treeFormat:
	file=open(pathSimu+treeFile)
	line=file.readline()
	file.close()
	i=0
	treeList=[]
	while line[i]!=";":
		treeList.append(line[i])
		if line[i]==")":
			i+=1
			while line[i]!=":" and line[i]!=";":
				i+=1
		else:
			i+=1
	file=open(pathSimu+treeFile+"_newFormat.tree","w")
	file.write(''.join(treeList)+"\n")
	file.close()
		
		
#Run simulations using pyvolve if requested
if pyvolveSim:
	start = time.time()

	import pyvolve
	
	if treeFormat:
		pyvolveTree = pyvolve.read_tree(file = pathSimu+treeFile+"_newFormat.tree", scale_tree = args.scale)
	else:
		pyvolveTree = pyvolve.read_tree(file = pathSimu+treeFile, scale_tree = args.scale)
		
	#nucModel=pyvolve.Model("nucleotide", {"matrix":mutMatrix})
	nucModel=pyvolve.Model("custom", {"matrix":mutMatrix})
	
	partitions=pyvolve.Partition(models = nucModel, root_sequence = ref)
	
	my_evolver = pyvolve.Evolver(tree = pyvolveTree, partitions = partitions)
	
	my_evolver(seqfile=pathSimu+pyvolveFastaO)
	
	elapsedTime = time.time() - start
	print("Time for simulating with pyvolve: "+str(elapsedTime))
	
	
	

exit()






