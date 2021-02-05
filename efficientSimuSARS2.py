import sys
import os
#import math
import numpy as np
import os.path
#from os import path
import argparse
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

	#	codonRates= np.zeros(nCodons,9)
	#CODONS: maybe don't create all the matrices from the start (might have too large an memory and time preparation cost).
	#instead, initialize only the rates from the reference allele (only 9 rates are needed), and store them in a dictionary at level 0 terminal nodes, and when new codons at a position 
	#are reached, extend the dictionary and calculate these new rates. Most positions will have only a few codons explored.
	
	genome_tree = phastSim.GenomeTree(
		nCodons=sim_run.nCodons,
		codon=sim_run.args.codon,
		ref=ref,
		gammaRates=gammaRates,
		omegas=omegas,
		mutMatrix=mutMatrix,
		hyperCategories=hyperCategories,
		hyperMutRates=sim_run.args.hyperMutRates,
		file=file)


	
	genome_tree.populateGenomeTree(node=genome_tree.genomeRoot)
	
	#I am assuming the branch lengths are in number of substitutions per nucleotide, even though we might be simulating a codon model.
	norm=genome_tree.genomeRoot.rate/len(ref)
	if verbose:
		print("\n Total cumulative mutation rate per site before normalization: "+str(norm))
	#We rescale by the input normalization factor, this is the same as rescaling all the branch lengths by this rescaling factor
	norm=norm/scale
	
	if codon:
		#function to iteratively normalize all rates
		def normalize(node,norm):
			node.rate/=norm
			if node.isTerminal:
				node.rates[node.allele]/=norm
			else:
				normalize(node.belowNodes[0],norm)
				normalize(node.belowNodes[1],norm)
		normalize(genome_tree.genomeRoot,norm)
	else:
		#function to iteratively normalize all rates
		def normalize(node,norm):
			node.rate/=norm
			if node.isTerminal:
				node.rates/=norm
			else:
				normalize(node.belowNodes[0],norm)
				normalize(node.belowNodes[1],norm)
		normalize(genome_tree.genomeRoot,norm)
	
	alRange=range(nAlleles)
	range9=range(9)
	
	#sample a new allele to mutate to, given the mutation rates at a current node
	#NICOLA: THIS MIGHT BE MADE FASTER, ESPECIALLY WITH LARGE STATE SPACES (E.G. CODON MODELS). 
	#INSTEAD OF ITERATING OVER ALL STATES, IN FACT, BELOW ONE COULD USE A DIVIDE AND CONQUER APPROACH OVER THE STATE SPACE, FOR EXAMPLE DEFINING ANOTHER CONSTANT TREE STRUCTURE OVER ALLELE SPACE.
	def sampleMutationCodon(rates,rand):
		for j in range9:
				if rand<rates[j]:
					return j
				else:
					rand-=rates[j]
		print("You should not have got here - there was a bug somewhere or some unlucky sampling in the machine error area")
		exit()
	def sampleMutation(allele,rates,rand):
		for j in alRange:
			if j!=allele:
				if rand<rates[j]:
					return j
				else:
					rand-=rates[j]
		print("You should not have got here - there was a bug somewhere or some unlucky sampling in the machine error area")
		exit()
			
	#NOW DO THE ACTUAL SIMULATIONS. DEFINE TEMPORARY STRUCTURE ON TOP OF THE CONSTANT REFERENCE GENOME TREE.
	#define a multi-layered tree; we start the simulations with a genome tree.
	#as we move down the phylogenetic tree, new layers are added below the starting tree. Nodes to layers below link to nodes above, or to nodes on the same layer, but never to nodes in the layer below.
	#while traversing the tree, as we move up gain from a node back to its parent (so that we can move to siblings etc), the nodes in layers below the current one are simply "forgotten" (in C they could be de-allocated, but the task here is left to python automation).
	
	#find position to mutate along the genome, and update temporary genome tree structure as you go
	def findPos(rand,parentGenomeNode,level):
		if parentGenomeNode.isTerminal:
			#reached a terminal node, now sample the mutation event at the position and update all rates
			node=parentGenomeNode.refNode
			a=parentGenomeNode.allele
			if codon:
				j=sampleMutationCodon(node.rates[a],rand)
				indeces=genome_tree.codonIndices[a]
				i2=int(j/3)
				i3=j%3
				newIndeces=list(indeces)
				newIndeces[i2]=(newIndeces[i2]+i3+1)%4
				parentGenomeNode.allele=genome_tree.codonIndices2[newIndeces[0],newIndeces[1],newIndeces[2]]
				mutEvent=[node.genomePos[0]*3+i2,indeces[i2],newIndeces[i2]]
				if verbose:
					print("Mutation from "+str(a)+" "+genome_tree.codonAllelesList[a]+" to "+str(parentGenomeNode.allele)+" "+genome_tree.codonAllelesList[parentGenomeNode.allele]+" , position "+str(mutEvent[0])+" category rate "+str(gammaRates[mutEvent[0]])+" hyperCat "+str(hyperCategories[mutEvent[0]])+" omega "+str(omegas[node.genomePos[0]])+" old rate "+str(node.rates[a][9])+" old rates:")
					print(node.rates[a])
				if not( parentGenomeNode.allele in node.rates ):
					node.rates[parentGenomeNode.allele]=np.zeros(10)
					indeces=genome_tree.codonIndices[parentGenomeNode.allele]
					parentGenomeNode.rate=0.0
					for i2 in range(3):
						pos2=node.genomePos[0]*3+i2
						nuc1=indeces[i2]
						for i3 in range(3):
							nuc2=(nuc1+i3+1)%4
							if genome_tree.isNonsynom[parentGenomeNode.allele,i2,i3]:
								if genome_tree.isIntoStop[parentGenomeNode.allele,i2,i3]:
									node.rates[parentGenomeNode.allele][i2*3+i3]=0.0
								else:
									node.rates[parentGenomeNode.allele][i2*3+i3]=omegas[node.genomePos[0]]*mutMatrix[nuc1][nuc2]*gammaRates[pos2]/norm
							else:
								node.rates[parentGenomeNode.allele][i2*3+i3]=mutMatrix[nuc1][nuc2]*gammaRates[pos2]/norm
						if hyperCategories[pos2]>0:
							if node.hyper[i2][0]==nuc1:
								node.rates[parentGenomeNode.allele][i2*3+node.hyper[i2][1]]*=hyperMutRates[hyperCategories[pos2]-1]
						for i3 in range(3):
							node.rates[parentGenomeNode.allele][9]+=node.rates[parentGenomeNode.allele][i2*3+i3]
				parentGenomeNode.rate=node.rates[parentGenomeNode.allele][9]
				if verbose:
					print(" new rate "+str(parentGenomeNode.rate)+" all rates:")
					print(node.rates[parentGenomeNode.allele])		
			else:
				j=sampleMutation(a,node.rates[a],rand)
				mutEvent=[node.genomePos[0],a,j]
				if verbose:
					print("Mutation from "+str(a)+" to "+str(j)+" , position "+str(node.genomePos[0])+" category rate "+str(gammaRates[node.genomePos[0]])+" hyperCat "+str(hyperCategories[node.genomePos[0]])+" old rate "+str(parentGenomeNode.rate)+" old rates:")
					print(node.rates)
				parentGenomeNode.rate=-node.rates[j,j]
				if verbose:
					print(" new rate "+str(parentGenomeNode.rate)+" all rates:")
					print(node.rates)
				parentGenomeNode.allele=j
			return mutEvent

		else:
			#still at an internal genome node.
			#choose which of the two children genome nodes to move into
			if rand>=parentGenomeNode.belowNodes[0].rate:
				rand-=parentGenomeNode.belowNodes[0].rate
				parentGenomeNode.rate=parentGenomeNode.belowNodes[0].rate
				child=parentGenomeNode.belowNodes[1]
				childI=1
			else:
				child=parentGenomeNode.belowNodes[0]
				parentGenomeNode.rate=parentGenomeNode.belowNodes[1].rate
				childI=0
			#if the child we are moving into is not on the same level, but is above, then create a new child at the same level.
			#this is because the rate of the child will be inevitably changed by the mutation event, and we don't want to change the mutation rates for the parent phylogenetic node.
			if child.level<level:
				newChild=phastSim.genomeNode(level=level) #upNode=parentGenomeNode
				parentGenomeNode.belowNodes[childI]=newChild
				newChild.isTerminal=child.isTerminal
				if child.isTerminal:
					newChild.refNode=child.refNode
					newChild.allele=child.allele
				else:
					newChild.belowNodes=list(child.belowNodes)
				
				mutEvent=findPos(rand,newChild,level)
				parentGenomeNode.rate+=newChild.rate
			else:
				#in this case the child is already on the same level, so no need to create another one, just update its mutation rate.
				mutEvent=findPos(rand,child,level)
				parentGenomeNode.rate+=child.rate
			return mutEvent
		
	
	#Function to simulate evolution on one branch,using ETE tree structure and using the genome-wide hierarchy structure.
	# given details of the parent node, it generates details of the child node, and updates the hierarchy accordingly.
	#To simulate evolution on the whole tree, it needs to be called on the root.
	def mutateBranchETEhierarchy(childNode,parentGenomeNode,level):
		#branch length above the current node
		bLen=childNode.dist
		currTime=0.0
		#if newick output is requested, prepare format
		if createNewick:
			childNode.mutAnnotation=[]
		#Initialize child rate and allele numbers with parent ones
		rate=parentGenomeNode.rate
		childNode.mutations=[]
	
		#Sample new mutation event with Gillespie algorithm
		currTime+=np.random.exponential(scale=1.0/rate)
		if verbose:
			print("\n Node "+childNode.name+" BLen: "+str(bLen)+" first sampled time: "+str(currTime)+" ; mutation rate: "+str(rate))
		#for the first mutation event at this node, create a new root genome node of the appropriate level. 
		#otherwise, use the one you already have.
		if currTime<bLen:
			newGenomeNode =phastSim.genomeNode(level=level)
			newGenomeNode.belowNodes=list(parentGenomeNode.belowNodes)
		else:
			newGenomeNode=parentGenomeNode
		while currTime<bLen:
			#Now, sample which type of mutation event it is (from which nucleotide to which nucleotide)
			rand=np.random.random()*rate
			if verbose:
				print("Selecting new mutation event. Rate "+str(rate)+" random value "+str(rand))
			mutEvent=findPos(rand,newGenomeNode,level)
			childNode.mutations.append(mutEvent)
			if createNewick:
				childNode.mutAnnotation.append(allelesList[mutEvent[1]]+str(mutEvent[0]+1)+allelesList[mutEvent[2]])
			rate=newGenomeNode.rate
			if verbose:
				print("New total rate "+str(rate))
			currTime+=np.random.exponential(scale=1.0/rate)
			if verbose:
				print("new time "+str(currTime)+", rate "+str(rate)+" mutation events:")
				print(childNode.mutations)
		
		if verbose:
			print("mutations at the end:")
			print(childNode.mutations)
		#now mutate children of the current node, calling this function recursively on the node children.
		for c in childNode.children:
			mutateBranchETEhierarchy(c,newGenomeNode,level+1)
	



#use simpler approach that collates same rates along the genome - less efficient with more complex models.
else:
	#List of dictionaries that let you know at wich position (starting from 0) is the 1st A of the genome, the second A of the genome, etc (well, actually the 0th, the 1st, etc..).
	#here I could use arrays instead of dictionaries, it might be more efficient.
	positions=[]
	for c in range(sim_run.nCat):
		positions.append([[],[],[],[]])
	#Total mutation rates cumulatively across the genome
	totMutMatrix=np.zeros((sim_run.nCat,4,4),dtype=float)
	#hyper rates of currently affected sites
	extras=[]
	#total mutation rate
	totMut=0.0
	#tot number of A's, C's, G's and T's in the reference genome
	#Number of alleles in each category across the genome.
	totAlleles=np.zeros((sim_run.nCat,4))
	for pos in range(len(ref)):
		a=alleles[ref[pos]]
		cat=sim_run.categories[pos]
		hyp=hyperCategories[pos]
		file.write(str(pos+1)+"\t"+str(cat)+"\t"+str(hyp)+"\t")
		#now sample which alleles are affected by hypermutation and store info in positions and extra vectors
		if hyp>0:
			i=np.random.choice(4)
			j=np.random.choice(3)
			j=(i+j+1)%4
			positions[cat][a].append([pos,hyp,i,j])
			if i==a:
				extras.append([mutMatrix[a][j]*categoryRates[cat]*(hyperMutRates[hyp-1]-1.0),pos,len(positions[cat][a])-1,cat,hyp,i,j])
				totMut+=mutMatrix[a][j]*categoryRates[cat]*(hyperMutRates[hyp-1]-1.0)
			file.write(allelesList[i]+"\t"+allelesList[j]+"\n")
		else:
			positions[cat][a].append([pos,0])
			file.write(".\t"+".\n")
		for j in range(4):
			if j!=a:
				totMutMatrix[cat][a][j]+=mutMatrix[a][j]*categoryRates[cat]
				totMut+=mutMatrix[a][j]*categoryRates[cat]
		totAlleles[sim_run.categories[pos]][a]+=1
	file.close()
	
	print("\n Number of each nucleotide in the genome:")
	print(totAlleles)

	norm=totMut/len(ref)
	print("\n Total cumulative mutation rate per site before normalization: "+str(norm))
	#We rescale by the input normalization factor, this is the same as rescaling all the branch lengths by this rescaling factor
	norm=norm/scale
	print("\n After rescaling: "+str(norm))
		
	#Now normalize mutation rates so that, at the start, the expected number of substitutions per unit branch length is 1?
	for i in range(4):
		for j in range(4):
			mutMatrix[i][j]=mutMatrix[i][j]/norm
			for c in range(sim_run.nCat):
				if j!=i:
					totMutMatrix[c][i][j]=totMutMatrix[c][i][j]/norm
	for e in range(len(extras)):
		extras[e][0]=extras[e][0]/norm
	totMut=totMut/norm
	
	if verbose:
		print("Base-wise total mutation rates:")
		print(totMutMatrix)
		print("Normalized mutation rates:")
		print(mutMatrix)

	#Function to simulate evolution on one branch,using ETE tree structure;
	# given details of the parent node, it generates details of the child node.
	#To simulate whole tree, it needs to be called on root.
	def mutateBranchETE(childNode,parentMuts,parentTotAlleles,parentRate,extrasParent):
		bLen=childNode.dist
		currTime=0.0
		#if newick output is requested, prepare format
		if createNewick:
			childNode.mutAnnotation=[]
		#Initialize child rate and allele numbers with parent ones
		rate=parentRate
		childTotAlleles=[]
		childMutations=[]
		for c in range(sim_run.nCat):
			childTotAlleles.append(list(parentTotAlleles[c]))
			childMutations.append([[],[],[],[]])
			#Initialize child mutation list with parent one
			for i in range(4):
				for k in range(len(parentMuts[c][i])):
					childMutations[c][i].append(list(parentMuts[c][i][k]))
		extrasChild=[]
		for e in range(len(extrasParent)):
			extrasChild.append(list(extrasParent[e]))
	
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
			if verbose:
				print("Selecting new mutation event. Rate "+str(rate)+" random value "+str(rand))
			tot=0.0
			found=False
			hyperExtra=False
			for c in range(sim_run.nCat):
				for i in range(4):
					for j in range(4):
						if j!=i:
							tot+=childTotAlleles[c][i]*mutMatrix[i][j]*categoryRates[c]
							#print("i "+str(i)+" j "+str(j)+" tot "+str(tot))
							if rand<(tot):
								found=True
								break
					if found:
						break
				if found:
					break
			if not found:
				#now use extras vector for hypermutable sites
				for e in range(len(extrasChild)):
					tot+=extrasChild[e][0]
					if rand<(tot):
						found=True
						hyperExtra=True
						extra=extrasChild.pop(e)
						break
			if not found:
				print("Error in selecting mutation type")
				exit()
			
			if hyperExtra:
				#element already removed from extraChild list, now use info in extra to add, remove or modify entry from mutations list, amend total rates, etc.
				a=alleles[ref[extra[1]]]
				extraRate=extra[0]
				rate-=extraRate
				pos=extra[1]
				pos2=extra[2]
				c=extra[3]
				hyp=extra[4]
				i=extra[5]
				j=extra[6]
				if verbose:
					print("Hyermutation genome position "+str(pos+1)+" category "+str(c)+" hypCat "+str(hyp)+" allele "+str(i)+" allele position "+str(pos2)+" to "+str(j))
				if createNewick:
					childNode.mutAnnotation.append(allelesList[i]+str(pos+1)+allelesList[j])
				m=0
				included=False
				while m<len(childMutations[c][a]):
					if childMutations[c][a][m][0]==pos2:
						if verbose:
							print("Position had already mutated")
						if childMutations[c][a][m][1]!=i:
							print("Error, mutation should have been recorded as into hypermutable allele.")
							print(childMutations[c][a][m])
							print(extra)
							print(ref[extra[1]])
							print(c)
							print(a)
							exit()
						childMutations[c][a][m][1]=j
						included=True
						break
					elif childMutations[c][a][m][0]<pos2:
						m+=1
					else:
						if ref[pos]!=allelesList[i]:
							print("Error, hypermutable allele should have been the reference.")
							exit()
						childMutations[c][a].insert(m,[pos2,j])
						if verbose:
							print("New mutation has been introduced")
							print(childMutations[c][a][m])
						included=True
						break
				if not included:
					childMutations[c][a].append([pos2,j])
					
			
			else:
				#Now, sample the specific position of the genome (among those with the mutated alleles) that mutates.
				mutatedBasePos=np.random.randint(childTotAlleles[c][i])
				if verbose:
					print("mutation position "+str(mutatedBasePos)+" category "+str(c)+" allele "+str(i)+" to "+str(j))
		
				#in this case, we are mutating a position that was not mutated before. We need to add one entry to the mutation list without removing any old one.
				if mutatedBasePos<totAlleles[c][i]-len(childMutations[c][i]):
					#print("New mutation")
					newMutPos=mutatedBasePos
					m=0
					while m<len(childMutations[c][i]):
						if childMutations[c][i][m][0]<=newMutPos:
							newMutPos+=1
							m+=1
						else:
							break
					childMutations[c][i].insert(m,[newMutPos,j])
					if createNewick:
						childNode.mutAnnotation.append(allelesList[i]+str(positions[c][i][newMutPos][0]+1)+allelesList[j])
					
					infoExtra=positions[c][i][newMutPos]
					if infoExtra[1]>0:
						if infoExtra[2]==i:
							#remove hypermutation from extras
							if verbose:
								print("removing from extras "+str(infoExtra[0])+" "+str(i)+" "+str(j))
							removed=False
							for e in range(len(extrasChild)):
								if infoExtra[0]==extrasChild[e][1]:
									extra=extrasChild.pop(e)
									rate-=extra[0]
									removed=True
									break
							if not removed:
								print("Error, hyper mutation not removed")
								exit()
							elif verbose:
								print("Hyper mutation removed by normal mutation")
								print(extra)
						elif infoExtra[2]==j:
							extraRate=mutMatrix[i][j]*categoryRates[c]*(hyperMutRates[infoExtra[1]-1]-1.0)
							extrasChild.append([extraRate,infoExtra[0],newMutPos,c,infoExtra[1],infoExtra[2],infoExtra[3]])
							rate+=extraRate
									
				#in this case, we are mutating a position that was already mutated.
				#this means that one item from the mutation list needs to be removed or modified, and no item needs to be added.
				else:
					newMutatedBasePos=mutatedBasePos-(totAlleles[c][i]-len(childMutations[c][i]))
					if verbose:
						print("Modifying pre-existing mutation")
						print(newMutatedBasePos)
					added=False
					for i2 in range(4):
						if i2!=i:
							for m in range(len(childMutations[c][i2])):
								if childMutations[c][i2][m][1]==i:
									if newMutatedBasePos==0:
										if createNewick:
											childNode.mutAnnotation.append(allelesList[i]+str(positions[c][i2][childMutations[c][i2][m][0]][0]+1)+allelesList[j])
										infoExtra=positions[c][i2][childMutations[c][i2][m][0]]
										if infoExtra[1]>0:
											if infoExtra[2]==i:
												#remove hypermutation from extras
												for e in range(len(extrasChild)):
													if infoExtra[0]==extrasChild[e][1]:
														extra=extrasChild.pop(e)
														rate-=extra[0]
														break
											elif infoExtra[2]==j:
												extraRate=mutMatrix[i][j]*categoryRates[c]*(hyperMutRates[infoExtra[1]-1]-1.0)
												extrasChild.append([extraRate,infoExtra[0],childMutations[c][i2][m][0],c,infoExtra[1],infoExtra[2],infoExtra[3]])
												rate+=extraRate
										if j==i2:
											if verbose:
												print("Deleted mutation \n\n\n")
											del childMutations[c][i2][m]
										else:
											childMutations[c][i2][m][1]=j
										added=True
										break
									newMutatedBasePos-=1
						if added:
							break
		
			childTotAlleles[c][i]-=1
			childTotAlleles[c][j]+=1
			rate+=mutMatrix[i][i]*categoryRates[c]
			rate-=mutMatrix[j][j]*categoryRates[c]
				
			currTime+=np.random.exponential(scale=1.0/rate)
			if verbose:
				print("new time "+str(currTime)+", rate "+str(rate)+" mutation events:")
				print(childMutations)
		
		if verbose:
			print("mutations at the end:")
			print(childMutations)
		childNode.mutations=childMutations
		#now mutate children of the current node, calling this function recursively on the node children.
		for c in childNode.children:
			mutateBranchETE(c,childMutations,childTotAlleles,rate,extrasChild)




time2 = time.time() - start
print("Total time after preparing for simulations: "+str(time2))




#Run sequence evolution simulation along tree
if hierarchy:
	mutateBranchETEhierarchy(t,genome_tree.genomeRoot,1)
else:
	muts=[]
	for c in range(sim_run.nCat):
		muts.append([[],[],[],[]])
	mutateBranchETE(t,muts,totAlleles,totMut,extras)
time2 = time.time() - start
print("Total time after simulating sequence evolution along tree with Gillespie approach: "+str(time2))





#Create a succint output
file=open(pathSimu+outputFile+".txt","w")

if hierarchy:
	mutDict={}
	#function to write a succint output iteratively
	def writeGenomeShort(node,file,mutDict):
		#update dictionary
		for m in node.mutations:
				nuc=allelesList[m[2]]
				if nuc!=ref[m[0]]:
					mutDict[m[0]+1]=nuc
				else:
					del mutDict[m[0]+1]
		#print leaf entry to file
		if node.is_leaf():
			file.write(">"+node.name+"\n")
			mutList=list(mutDict.keys())
			mutList.sort()
			for m in mutList:
				file.write(str(m)+" "+mutDict[m]+"\n")
		#pass dictionary to children
		else:
			for c in node.children:
				writeGenomeShort(c,file,mutDict)
		#de-update the dictionary so it can be used by siblings etc. 
		for n in range(len(node.mutations)):
			m=node.mutations[len(node.mutations)-(n+1)]
			nuc=allelesList[m[1]]
			if nuc!=ref[m[0]]:
				mutDict[m[0]+1]=nuc
			else:
				del mutDict[m[0]+1]
	writeGenomeShort(t,file,mutDict)
else:
	#function to write a succint output iteratively
	def writeGenomeShort(node,file):
		if node.is_leaf():
			file.write(">"+node.name+"\n")
			mutDict={}
			for c in range(sim_run.nCat):
				for i in range(4):
					for m in node.mutations[c][i]:
						mutDict[positions[c][i][m[0]][0]+1]=allelesList[m[1]]
			mutList=list(mutDict.keys())
			mutList.sort()
			for m in mutList:
				file.write(str(m)+" "+mutDict[m]+"\n")
		for c in node.children:
			writeGenomeShort(c,file)
			
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
	file=open(pathSimu+outputFile+".tree","w")
	newickTree=writeGenomeNewick(t)+";\n"
	file.write(newickTree)
	file.close()
	
	time3 = time.time() - start
	print("Total time after writing newick file: "+str(time3))






if createFasta or createPhylip:
	if not hierarchy:
		#generate the genome sequence of a sample using node mutations and the reference
		#useful, for example, for generating a fasta file.
		def genomeSeq(mutations):
			for c in range(sim_run.nCat):
				for i in range(4):
					for m in mutations[c][i]:
						refList[positions[c][i][m[0]][0]]=allelesList[m[1]]
			newGenome=''.join(refList)
			for c in range(sim_run.nCat):
				for i in range(4):
					for m in mutations[c][i]:
						pos=positions[c][i][m[0]][0]
						refList[pos]=ref[pos]
			return newGenome


#If requested, create a fasta output
if createFasta:
	file=open(pathSimu+outputFile+".fasta","w")
	
	if hierarchy:
		#function to write a complete sequence output iteratively
		def writeGenome(node,file,nRefList):
			#update list
			for m in node.mutations:
				nRefList[m[0]]=allelesList[m[2]]
			#print leaf entry to file
			if node.is_leaf():
				file.write(">"+node.name+"\n"+(''.join(nRefList))+"\n")
			#pass to children
			else:
				for c in node.children:
					writeGenome(c,file,nRefList)
			#de-update the list so it can be used by siblings etc. 
			for n in range(len(node.mutations)):
				m=node.mutations[len(node.mutations)-(n+1)]
				nRefList[m[0]]=allelesList[m[1]]
		writeGenome(t,file,refList)

	else:
		#function to write a fasta output iteratively
		def writeGenome(node,file):
			if node.is_leaf():
				seq=genomeSeq(node.mutations)
				file.write(">"+node.name+"\n"+seq+"\n")
			for c in node.children:
				writeGenome(c,file)
		writeGenome(t,file)
	file.close()
	
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








