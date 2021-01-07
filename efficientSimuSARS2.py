import sys
import os
import math
import numpy as np
import os.path
from os import path
import argparse
from ete3 import Tree
import time

#NOW RATE VARIATION AND HYPERMUTABILITY

parser = argparse.ArgumentParser(description='Run simulations about the effective alignment size (EASi).')
parser.add_argument('--path',default="", help='Path where to run simulations.')
parser.add_argument('--reference',default="MN908947.3.fasta", help='File containing the reference genome to be used as root genome. To be found in the folder specified with --path.')
parser.add_argument('--treeFile',default="exampleTree.tree", help='Name of file containing the tree used to simulate sequences (assumed within the --path and in newick format).')
parser.add_argument('--scale',default=1.0,type=float, help='Scale the simulation tree by this amount (default 1.0). Branch lengths are assumed in terms of expected substitutions per site (more or less, as frequencies changes through time, total mutation rate might also  change).')
parser.add_argument("--seed", help="Seed for random simulator in genes simulations", type=int, default=1)
parser.add_argument("--mutationRates", help="Mutation rates, by default using the neutral rates estimated from SARS-CoV-2; so far only exactly 12 input values allowed (r_AC, r_AG, r_AT, r_CA, etc...) corresponding to an UNREST nucleotide substitution model.", type=float, nargs='+',  default=[0.039, 0.310, 0.123, 0.140, 0.022, 3.028, 0.747, 0.113, 2.953, 0.056, 0.261, 0.036])
parser.add_argument("--categoryProbs", help="Probabilities of rate categories. They are supposed to sum up to one, but if they don't they are normalized to do so. By default only 1 category is simulated.", type=float, nargs='+',  default=[1.0])
parser.add_argument("--categoryRates", help="Rates of site categories. The overall mutation rate is renormalized so that the expected number of bstitutions per branch length unit per site is 1. By default only 1 category of rate 1.0 is simulated. The number of rates has to be the same as the number of categgory probabilities, otherwise an error is thrown.", type=float, nargs='+',  default=[1.0])
parser.add_argument("--hyperMutProbs", help="Proportions of sites with frequently recurring mutations. They are supposed to be rare, and not sum up to one. By default no recurring mutations are simulated.", type=float, nargs='+',  default=[])
parser.add_argument("--hyperMutRates", help="Rates of recurring mutations. Different classes of recurring mutations with different hypermutabilities can be specified. By defaultno recurrent mutations are simulated. The number of rates has to be the same as the number of recurring mutation site proportions, otherwise an error is thrown.", type=float, nargs='+',  default=[])
parser.add_argument("--verbose", help="Turns on verbose mode.", action="store_true")
parser.add_argument('--outputFile',default="sars-cov-2_simulation_output", help='Output file name containing the simulated genomes in succint format. The file will be created within the folder specified with --path.')
parser.add_argument("--createNewick", help="Create a newick file annotated with the simulated mutation events (default name sars-cov-2_simulation_output.tree).", action="store_true")
#parser.add_argument('--newickO',default="sars-cov-2_simulation_output.tree", help='Output file name containing the tree annotated with simulated mutation events. The file will be created within the folder specified with --path.')
parser.add_argument("--createFasta", help="Create a fasta file with the simulated genomes (default name sars-cov-2_simulation_output.fasta).", action="store_true")
#parser.add_argument('--fastaO',default="sars-cov-2_simulation_output.fasta", help='Output file name containing the simulated genomes. The file will be created within the folder specified with --path.')
parser.add_argument("--createPhylip", help="Create a phylip file with the simulated genomes (default name sars-cov-2_simulation_output.phy).", action="store_true")
#parser.add_argument('--phylipO',default="sars-cov-2_simulation_output.phy", help='Output phylip file name containing the simulated genomes. The file will be created within the folder specified with --path.')
parser.add_argument("--pyvolveSim", help="run simulations using pyvolve (not the default, it's only to allow easy comparison)", action="store_true")
#parser.add_argument('--pyvolveFastaO',default="pyvolve_simulation_output.fasta", help='Output file name containing the genomes simulated by pyvolve, if pyvolve is used. The file will be created within the folder specified with --path.')
parser.add_argument("--treeFormat", help="remove branch support values from tree file to use in pyvolve", action="store_true")
args = parser.parse_args()

pathSimu=args.path
reference=args.reference
np.random.seed(args.seed)
seed=args.seed
scale=args.scale
mutationRates=args.mutationRates
categoryProbs=args.categoryProbs
categoryRates=args.categoryRates
hyperMutProbs=args.hyperMutProbs
hyperMutRates=args.hyperMutRates
verbose=args.verbose
outputFile=args.outputFile
createNewick=args.createNewick
#newickO=args.newickO
createFasta=args.createFasta
#fastaO=args.fastaO
createPhylip=args.createPhylip
#phylipO=args.phylipO
pyvolveSim=args.pyvolveSim
#pyvolveFastaO=args.pyvolveFastaO
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

#print information about a tree node
def printNode(node):
	print("\n Node name " + node.name+" , sampled="+str(node.is_leaf()))
	print(str(len(node.children))+ " children, distance from parent "+str(node.dist))
	print("Mutations:")
	for i in range(4):
		print(node.mutations[i])


nCat=len(categoryProbs)
nHyper=len(hyperMutProbs)

if nCat!=len(categoryRates):
	print("Issue with number of category probs "+str(len(categoryProbs))+" and number of category rates "+str(len(categoryRates)))
	exit()
if nHyper!=len(hyperMutRates):
	print("Issue with number of hypermutation category probs "+str(len(hyperMutProbs))+" and number of hypermutation category rates "+str(len(hyperMutRates)))
	exit()
	
sum=0.0
for i in categoryProbs:
	sum+=i
for i in range(nCat):
	categoryProbs[i]=categoryProbs[i]/sum
if sum>1.000001 or sum<0.999999:
	print("\n Normalizing probabilities of site categories. New probabilities:")
	print(categoryProbs)
for i in hyperMutRates:
	if i<=1.0:
		print("It doesn't make sense to have hypermutability class with mutability <=1.0 . hyperMutRates:")
		print(hyperMutRates)
		exit()
sumHyper=0.0
for i in hyperMutProbs:
	sumHyper+=i
if sumHyper>0.1:
	print("WARNING: hypermutable sites are supposed to be rare, but total proportion is "+str(sum))
if sumHyper>1.0:
	exit()
newHyperMutProbs=[1.0-sumHyper]+hyperMutProbs
	
#sample category for each site of the genome
categories=np.random.choice(nCat,size=len(ref),p=categoryProbs)
#sample hypermutability for each site of the genome.
hyperCategories=np.random.choice(nHyper+1,size=len(ref),p=newHyperMutProbs)


#MAYBE TRY DIFFERENT APPROACH WITH DIVIDE ET IMPERA ALONG THE GENOME,
# WHEN THE RATES ARE UPDATED, UPDATE ONLY RATE OF THE SITE AND OF ALL THE NODES ON TOP OF THE HYRARCHY.
# THAT IS, DEFINE A PIRAMID STRUCTURE, WITH TERMINAL NODES BEING GENOME LOCI AND WITH INTERNAL NODES BING MERGING GROUPS OF LOCI CONTAINING INFORMATION ABOUT THEIR CUMULATIVE RATES.
# THIS WAY UPDATING A MUTATION EVENT REQUIRES COST LOGARITHMIC IN GENOME SIZE EVEN IF EVERY SITE HAS A DIFFERENT RATE. 

#List of dictionaries that let you know at wich position (starting from 0) is the 1st A of the genome, the second A of the genome, etc (well, actually the 0th, the 1st, etc..).
#here I could use arrays instead of dictionaries, it might be more efficient.
positions=[]
for c in range(nCat):
	positions.append([[],[],[],[]])
#Total mutation rates cumulatively across the genome
totMutMatrix=np.zeros((nCat,4,4),dtype=float)
#hyper rates of currently affected sites
extras=[]
#total mutation rate
totMut=0.0
#tot number of A's, C's, G's and T's in the reference genome
#Number of alleles in each category across the genome.
totAlleles=np.zeros((nCat,4))
#save information about the categories of each site on a file
file=open(pathSimu+outputFile+".info","w")
file.write("pos\t"+"cat\t"+"hyperCat\t"+"hyperAlleleFrom\t"+"hyperAlleleTo\n")
for pos in range(len(ref)):
	a=alleles[ref[pos]]
	cat=categories[pos]
	hyp=hyperCategories[pos]
	file.write(str(pos+1)+"\t"+str(cat)+"\t"+str(hyp)+"\t")
	#now sample which alleles are affected by hypermutation and store info in positions and extra vectors
	if hyp>0:
		i=np.random.choice(4)
		j=np.random.choice(3)
		j=(i+j+1)%4
		positions[cat][a].append([pos,hyp,i,j])
		if i==a:
			#here we store the contribution of hypermutability
			#print("Before: "+str(totMut)+" "+str(i)+" "+str(j))
			#print(hyperMutRates[hyp-1])
			#print(mutMatrix[a][j])
			#print(categoryRates[cat])
			extras.append([mutMatrix[a][j]*categoryRates[cat]*(hyperMutRates[hyp-1]-1.0),pos,len(positions[cat][a])-1,cat,hyp,i,j])
			totMut+=mutMatrix[a][j]*categoryRates[cat]*(hyperMutRates[hyp-1]-1.0)
			#print("After: "+str(totMut)+"\n")
		file.write(allelesList[i]+"\t"+allelesList[j]+"\n")
	else:
		positions[cat][a].append([pos,0])
		file.write(".\t"+".\n")
	for j in range(4):
		if j!=a:
			totMutMatrix[cat][a][j]+=mutMatrix[a][j]*categoryRates[cat]
			totMut+=mutMatrix[a][j]*categoryRates[cat]
	totAlleles[categories[pos]][a]+=1
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
		for c in range(nCat):
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
		
			
#generate the genome sequence of a sample using node mutations and the reference
#useful, for example, for generating a fasta file.
def genomeSeq(mutations):
	#seqList=list(refList)
	#print(mutations)
	for c in range(nCat):
		for i in range(4):
			for m in mutations[c][i]:
				#print("\n"+str(i)+" "+str(m[0])+" "+str(m[1]))
				refList[positions[c][i][m[0]][0]]=allelesList[m[1]]
	newGenome=''.join(refList)
	for c in range(nCat):
		for i in range(4):
			for m in mutations[c][i]:
				#print(m)
				#print("\n"+str(i)+" "+str(m[0])+" "+str(m[1]))
				pos=positions[c][i][m[0]][0]
				refList[pos]=ref[pos]
	return newGenome


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
	for c in range(nCat):
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
		#print("Random value: "+str(rand))
		for c in range(nCat):
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
			#positions[cat][a].append([pos,hyp,i,j])
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
					#childMutations[c][i].insert(m,[newMutPos,j])
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
				#for m in range(len(childMutations[i])):
					if childMutations[c][i][m][0]<=newMutPos:
						newMutPos+=1
						m+=1
					else:
						break
				#print(m)
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
				#exit()
				added=False
				for i2 in range(4):
					if i2!=i:
						for m in range(len(childMutations[c][i2])):
							if childMutations[c][i2][m][1]==i:
								if newMutatedBasePos==0:
									if createNewick:
										#print(childMutations[c][i2][m])
										#print(positions[c][i2][childMutations[c][i2][m][0]])
										#print(allelesList[i])
										#print(allelesList[i]+str(positions[c][i2][childMutations[c][i2][m][0]][0]+1)+allelesList[j])
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





start = time.time()

# Loads a tree structure from a newick string in ETE2. The returned variable t is the root node for the tree.
t = Tree(pathSimu+treeFile)
time1 = time.time() - start
print("Time for reading tree with ETE3: "+str(time1))

muts=[]
for c in range(nCat):
	muts.append([[],[],[],[]])
#Run sequence evolution simulation along tree
mutateBranchETE(t,muts,totAlleles,totMut,extras)
time2 = time.time() - start
print("Total time after simulating sequence evolution along tree with Gillespie approach: "+str(time2))






#function to write a succint output iteratively
def writeGenomeShort(node,file):
	if node.is_leaf():
		file.write(">"+node.name+"\n")
		mutDict={}
		for c in range(nCat):
			for i in range(4):
				for m in node.mutations[c][i]:
					mutDict[positions[c][i][m[0]][0]+1]=allelesList[m[1]]
		mutList=list(mutDict.keys())
		mutList.sort()
		for m in mutList:
			file.write(str(m)+" "+mutDict[m]+"\n")
	for c in node.children:
		writeGenomeShort(c,file)

#Create a succint output
file=open(pathSimu+outputFile+".txt","w")
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

#function to write a fasta output iteratively
def writeGenome(node,file):
	if node.is_leaf():
		seq=genomeSeq(node.mutations)
		file.write(">"+node.name+"\n"+seq+"\n")
	for c in node.children:
		writeGenome(c,file)

#If requested, create a fasta output
if createFasta:
	file=open(pathSimu+outputFile+".fasta","w")
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
	file=open(pathSimu+outputFile+".phy","w")
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
	
	my_evolver(seqfile=pathSimu+outputFile+"_pyvolve.txt")
	
	elapsedTime = time.time() - start
	print("Time for simulating with pyvolve: "+str(elapsedTime))
	
	
	

exit()


































#OLD VERSION
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

#generate the genome sequence of a sample using node mutations and the reference
#useful, for example, for generating a fasta file.
def genomeSeq(mutations):
	seqList=list(refList)
	#print(mutations)
	for i in range(4):
			for m in mutations[c][i]:
				#print("\n"+str(i)+" "+str(m[0])+" "+str(m[1]))
				seqList[positions[i][m[0]]]=allelesList[m[1]]
	return ''.join(refList)

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

