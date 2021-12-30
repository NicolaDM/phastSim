#import sys
import os
import subprocess
#import math
import numpy as np
import os.path
#from os import path
import argparse
#from ete3 import Tree
import time

#script to compare the running time of phastSim with pyvolve, seq-gen and indelible.

# python compareSimulators.py --nLeaves 100 --createFasta --pyvolveSim --seqgenSim --indelibleSim --indelibleSim2

parser = argparse.ArgumentParser(description='Compare simulators, in particular phastSim, seq-gen, indelible and pyvolve.')
#parser.add_argument('--path',default="/Users/demaio/Desktop/coronavirus/simulations/", help='Path where to write simulation output.')
#parser.add_argument('--reference',default="/Users/demaio/Documents/GitHub/phastSim/phastSim/example/MN908947.3.fasta", help='File containing the reference genome to be used as root genome.')
parser.add_argument('--path',default="/home/will/Desktop/projects/embl/phastSim/simulation_output_8/", help='Path where to write simulation output.')
parser.add_argument('--indeliblePath', default="/home/will/Downloads/INDELibleV1.03/bin/indelible")
parser.add_argument('--seqgenPath', default="seq-gen")
parser.add_argument('--pyvolvePath', default="/home/will/Desktop/projects/embl/phastSim/scripts/runPyvolve.py")
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
parser.add_argument("--monitorMemory", action="store_true")
args = parser.parse_args()

pathSimu=args.path
indeliblePath=args.indeliblePath
seqgenPath=args.seqgenPath
pyvolvePath=args.pyvolvePath
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

monitorMemory=args.monitorMemory

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
    

def monitor_memory(string_run):
    string_run = string_run.split(">")[0]
    result = subprocess.run(["/usr/bin/time", "-v"] +  string_run.split(), capture_output=True, text=True)
    print("*********************")
    print(string_run)
    print("**********************")
    print(result.stderr.split("\n"))
    memory_usage = int(result.stderr.split("\n")[-15].split(":")[-1])
    return memory_usage

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
memory=[[],[],[],[],[],[],[],[]]

for r in range(nReplicates):
    #times2=[]
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

    memory[0].append("NaN")
    

    

    
    
    stringRun="phastSim --outpath "+pathSimu+"  --seed "+str(seed+r)+" --treeFile "+pathSimu+treeFile+" --reference "+reference+" --outputFile simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+".txt" # --scale "+str(scale)+"
    #if createFasta:
    #    stringRun+=" --createFasta "
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
        
        #run phastSim
        print("Running phastSim")
        start = time.time()
        os.system(stringRun +" >/dev/null")
        #print(stringRun)
        time2 = time.time() - start
        times[1].append(time2)
        #times2.append(time2)
        print("Total time after simulating sequence evolution along tree with phastSim: "+str(time2))

        if monitorMemory:
            result = monitor_memory(stringRun)
            memory[1].append(result)

    else:
        times[1].append("NaN")
        memory[1].append("NaN")
    
    if createFasta:
        start = time.time()
        #print(stringRun)
        os.system(stringRun+" --createFasta"+" >/dev/null")
        time2 = time.time() - start
        #times2.append(time2)
        times[2].append(time2)
        print("Total time after simulating sequence evolution along tree with phastSim and writing fasta file: "+str(time2))

        if monitorMemory:
            result = monitor_memory(stringRun)
            memory[2].append(result)
    else:
        times[2].append("NaN")
        memory[2].append("NaN")
        
    if noHierarchy:
        start = time.time()
        os.system(stringRun+" --noHierarchy "+" >/dev/null")
        time2 = time.time() - start
        #times2.append(time2)
        times[3].append(time2)
        print("Total time after simulating sequence evolution along tree with phastSim and no multilayer genome tree: "+str(time2))

        if monitorMemory:
            result = monitor_memory(stringRun)
            memory[3].append(result)
    else:
        times[3].append("NaN")
        memory[3].append("NaN")
    
    #run pyvolve
    #Run simulations using pyvolve if requested
    if pyvolveSim:
        print("Running pyvolve")
        start = time.time()
                
        stringRun="python "+ pyvolvePath + " --path "+ str(args.path) +"  --seed "+str(seed+r)+" --scale "+str(scale)+" --treeFile "+treeFile2+" --outputFile simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_pyvolve.txt"

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

            if monitorMemory:
                result = monitor_memory(stringRun)
                memory[4].append(result)

    else:
        #times2.append("NaN")
        times[4].append("NaN")
        memory[4].append("NaN")



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
        stringRun=f"{seqgenPath} -l "+str(length)+" -m GTR -f 0.3 0.2 0.2 0.3 -r "+str(mutMatrix[0][1]/mutMatrix[2][3])+" "+str(mutMatrix[0][2]/mutMatrix[2][3])+" "+str(mutMatrix[0][3]/mutMatrix[2][3])+" "+str(mutMatrix[1][2]/mutMatrix[2][3])+" "+str(mutMatrix[1][3]/mutMatrix[2][3])+" 1.0 "+categoryString+" -z "+str(seed+r)+" "+pathSimu+treeFile+" > "+pathSimu+"simulationOutput_repl"+str(r+1)+"_nLeaves"+str(nLeaves)+"_scale"+str(scale)+"_seqgenOutput.txt" # -s "+str(scale)+"
        os.system(stringRun+" ") #&>/dev/null
        time2 = time.time() - start
        #times2.append(time2)
        times[5].append(time2)
        print("Total time after simulating sequence evolution along tree with seqgen: "+str(time2))
        if monitorMemory:
            result = monitor_memory(stringRun)
            memory[5].append(result)
    else:
        #times2.append("NaN")
        times[5].append("NaN")
        memory[5].append("NaN")
    
    
    #run simulations with indelible
    if indelibleSim:

        
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
                file.write(" [insertmodel]    NB    0.5    1 \n")
                file.write(" [deletemodel]    NB    0.5    1 \n")
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
                file.write(" [insertmodel]    NB    0.5    1 \n")
                file.write(" [deletemodel]    NB    0.5    1 \n")
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

        print("Running indelible method 1")
        start = time.time()
        #os.system("cd "+pathSimu+" ; "+"/Users/demaio/Desktop/INDELibleV1.03/bin/indelible ")
        os.system("cd "+pathSimu+" ; "+  indeliblePath+ " ")
        
        time2 = time.time() - start
        #times2.append(time2)
        times[6].append(time2)
        print("Total time after simulating sequence evolution along tree with indelible method 1: "+str(time2))
        if monitorMemory:
            memres = subprocess.run("cd "+pathSimu+" ; "+ " /usr/bin/time -v " + indeliblePath + " ", shell=True, text=True, capture_output=True)
            print("*********************")
            print(memres.stderr.split("\n"))
            print("**********************")
            memory_usage = int(memres.stderr.split("\n")[-15].split(":")[-1])
            memory[6].append(memory_usage)

    else:
        #times2.append("NaN")
        times[6].append("NaN")
        memory[6].append("NaN")
        
    if indelibleSim2:
        #now version with full Gillespie

        
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
                file.write(" [insertmodel]    NB    0.5    1 \n")
                file.write(" [deletemodel]    NB    0.5    1 \n")
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
                file.write(" [insertmodel]    NB    0.5    1 \n")
                file.write(" [deletemodel]    NB    0.5    1 \n")
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

        print("Running indelible method 2")
        start = time.time()

        #os.system("cd "+pathSimu+" ; "+"/Users/demaio/Desktop/INDELibleV1.03/bin/indelible ")
        os.system("cd "+pathSimu+" ; "+  indeliblePath + " ")
        
        time2 = time.time() - start
        #times2.append(time2)
        times[7].append(time2)
        print("Total time after simulating sequence evolution along tree with indelible method 2: "+str(time2))

        if monitorMemory:
            memres = subprocess.run("cd "+pathSimu+" ; "+ " /usr/bin/time -v " + indeliblePath + " ", shell=True, text=True, capture_output=True)
            #print("*********************")
            #print(memres.stderr.split("\n"))
            #print("**********************")
            memory_usage = int(memres.stderr.split("\n")[-15].split(":")[-1])
            memory[7].append(memory_usage)
    else:
        #times2.append("NaN")
        times[7].append("NaN")
        memory[7].append("NaN")
        
    
    #times.append(times2)

print("Total matrix of running times:")
print(times)

if monitorMemory:
    print("Total matrix of max memory usage:")
    print(memory)