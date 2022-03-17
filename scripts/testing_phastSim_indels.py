
import os
import numpy as np
from numpy import mean
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse

#©EMBL-European Bioinformatics Institute, 2022

#script that simulates indels in phastSim and compares their distribution to INDELible and makes plots showing the comparison in a few scenarios.
parser = argparse.ArgumentParser(description='Run simulations of fragments from the human genome to be sequenced pentamer-wise; try to assemble using the pentamers.')
parser.add_argument('--path',default="/Users/demaio/Desktop/coronavirus/phastSim_resumission_March2022/extensive_testing_indels/", help='path where to find files and plot results.')
args = parser.parse_args()
path=args.path



def errplot(times,labels,plotFileName,n_leaves,colors,topPlot,linestyles,title="Comparison of simulated indel distributions",yAxisLabel="Time (seconds)",logY=False,ymin=None,ymax=None,violin=True,xticks=None):
	
	values=[]
	mean_times = []
	errors = []
	for i in range(len(times[0])):
		values.append([])
		for j in range(len(times)):
			if len(times[j][i])>1:
				values[-1].append([])
				for k in range(len(times[j][i])):
					values[-1][-1].append(times[j][i][k])
	for t1 in times:
		mean_times.append([])
		errors.append([])
		for t2 in t1:
			print(t2)
			mean_times[-1].append(np.mean([float(t) for t in t2]))
			if violin:
				errors[-1].append(0.0)
			else:
				errors[-1].append(np.std([float(t) for t in t2]))
			#print([float(t) for t in t2])
			#print(np.mean([float(t) for t in t2]))
			#print(np.std([float(t) for t in t2]))
			#print("\n")

	mean_times = np.array(mean_times).T
	errors = np.array(errors).T    
	
	if violin:
		patches=[]
		for i in range(len(labels)):
			patches.append(mpatches.Patch(color=colors[i]))

	x = range(1,len(n_leaves)+1)
	wids=[]
	for xi in x:
		#wids.append(xi/4)
		wids.append(1.0/4)
	fig = plt.figure(figsize=(15, 9))
	ax1 = fig.add_subplot(111)
	ax1.set_axisbelow(True)
	ax1.set_title(title, fontsize=26)
	ax1.set_xlabel(topPlot, fontsize=22)
	ax1.set_ylabel(yAxisLabel, fontsize=22)
	#ax1.xticks(x, n_leaves, rotation=40)
	ax1.set_xticks(x, minor=False)
	ax1.set_xticklabels(n_leaves, rotation=40, fontdict=None, minor=False)
	#plt.xticks(fontsize=18)
	#plt.yticks(fontsize=18)
	ax1.tick_params(labelsize=16)

	if ymin!=None:
		ax1.set_ylim([ymin, ymax])
	#ax1.set_xscale('log')
	if logY:
		ax1.set_yscale('log')
	ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',   alpha=0.5)
	
	#reordered_numbers = [0, 1, 2, 4, 6, 5, 3]
	for i in range(len(labels)):
		#ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='x', c=colors[i], label=labels[i], capsize=2)
		if violin:
			#print(values[i])
			#print(x[:len(values[i])])
			violin_parts = ax1.violinplot(values[i], positions=x[:len(values[i])], vert=True, widths=wids[:len(values[i])], showmeans=False, showextrema=False, showmedians=False, quantiles=None, points=100, bw_method=None, data=None)
			
			for pc in violin_parts['bodies']:
				pc.set_facecolor(colors[i])
				pc.set_edgecolor(colors[i])
			#ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='.', c=colors[i], label=labels[i], capsize=2, linestyle=linestyles[i])
			if linestyles!=None:
				ax1.plot(x, mean_times[i, :], marker='.', c=colors[i], label=labels[i], linestyle=linestyles[i])
			else:
				ax1.plot(x, mean_times[i, :], marker='.', c=colors[i], label=labels[i])
		else:
			if linestyles!=None:
				ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='x', c=colors[i], label=labels[i], capsize=2, linestyle=linestyles[i])
			else:
				ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='x', c=colors[i], label=labels[i], capsize=2)
	#if violin:
	#	plt.legend(patches, labels, loc='upper left')
	#else:
	plt.legend(loc='upper left', prop={'size': 22})
	fig.savefig(plotFileName,bbox_inches='tight')




#baseLength=0.001
parameters=[(0.0001,1.0,0.3,1.0,0.3), (0.0003,1.0,0.3,1.0,0.3), (0.0001,3.0,0.3,1.0,0.3), (0.0001,1.0,0.3,3.0,0.3), (0.0001,1.0,0.1,1.0,0.3), (0.0001,1.0,0.3,1.0,0.1)]
names=["base","longBranches","highInsertion","highDeletion","longInsertion","longDeletion"]
numReplicates=0
minimumRep=50
runPhastSim=False
runIndelible=False
#rootGenomeLength=29903
rootGenomeLength=1000000
dataPhastSim=[]
dataIndelible=[]
for i in range(len(parameters)):
    print("Parameter "+str(i+1))
    dataPhastSim.append([])
    #newickTree2="(((S1:"+str(parameters[i][0])+",S2:"+str(parameters[i][0])+"):"+str(parameters[i][0])+",(S3:"+str(parameters[i][0])+",S4:"+str(parameters[i][0])+"):"+str(parameters[i][0])+"):"+str(parameters[i][0])+",((S5:"+str(parameters[i][0])+",S6:"+str(parameters[i][0])+"):"+str(parameters[i][0])+",(S7:"+str(parameters[i][0])+",S8:"+str(parameters[i][0])+"):"+str(parameters[i][0])+"):"+str(parameters[i][0])+");"
    #newickTree="(((S1:"+str(parameters[i][0])+",S2:"+str(parameters[i][0])+"):"+str(parameters[i][0])+",(S3:"+str(parameters[i][0])+",S4:"+str(parameters[i][0])+"):"+str(parameters[i][0])+"):"+str(parameters[i][0])+",((S5:"+str(parameters[i][0])+",S6:"+str(parameters[i][0])+"):"+str(parameters[i][0])+",(S7:"+str(parameters[i][0])+",S8:"+str(parameters[i][0])+"):"+str(parameters[i][0])+"):"+str(parameters[i][0])+"):"+str(parameters[i][0])+";"
    #newickTree2="(S1:"+str(parameters[i][0])+",S2:0.0);"
    newickTree="(S1:"+str(parameters[i][0])+",S2:0.0);"
    
    file=open(path+"tree_par"+str(i+1)+".tree","w")
    file.write(newickTree+"\n")
    file.close()
    for r in range(numReplicates+minimumRep):
        dataPhastSim[i].append(([],[]))
        #run phastSim
        if runPhastSim:
            if r>=minimumRep:
                print(r)
                os.chdir('/Users/demaio/Documents/GitHub/phastSim')
                os.system("python3 bin/phastSim --rootGenomeLength "+str(rootGenomeLength)+" --outpath "+path+" --seed "+str(r+1)+" --treeFile "+path+"tree_par"+str(i+1)+".tree --scale 1.0 --outputFile simulationPhastSim_par"+str(i+1)+"_rep"+str(r+1)+"  --indels --insertionRate CONSTANT "+str(parameters[i][1])+" --deletionRate CONSTANT "+str(parameters[i][3])+" --insertionLength GEOMETRIC "+str(parameters[i][2])+" --deletionLength GEOMETRIC "+str(parameters[i][4])+"  \n")
        #investigate phastSim output
        file=open(path+"simulationPhastSim_par"+str(i+1)+"_rep"+str(r+1)+".txt")
        line=file.readline()
        line=file.readline()
        while line!="\n" and line!="" and line[0]!=">":
            linelist=line.split()
            len1=len(linelist[0])
            len2=len(linelist[2])
            if len1==1 and len2>1:
                dataPhastSim[i][r][0].append(len2-1)
            elif len1>1:
                dataPhastSim[i][r][1].append(len1)
            elif len1==1 and linelist[2]=="-":
                dataPhastSim[i][r][1].append(1)
            line=file.readline()
        file.close()

    #run INDELible
    dataIndelible.append([])
    if runIndelible:
        file=open(path+"control.txt","w")
        file.write("[TYPE] NUCLEOTIDE 1\n")
        file.write("[MODEL]    modelname\n")
        file.write("  [submodel]     JC\n")
        file.write("  [insertmodel]   NB  "+str(1.0-parameters[i][2])+" 1\n")
        file.write("  [deletemodel]   NB  "+str(1.0-parameters[i][4])+" 1\n")
        file.write("  [insertrate]   "+str(parameters[i][1])+"\n")
        file.write("  [deleterate]   "+str(parameters[i][3])+"\n")
        file.write("[TREE] treename  "+newickTree+"\n")
        file.write("[PARTITIONS] partitionname\n")
        file.write("  [treename modelname "+str(rootGenomeLength)+"]\n")
        file.write("[EVOLVE] partitionname "+str(numReplicates+minimumRep)+" INDELible_sim_par"+str(i+1)+"\n")
        file.close()
        os.chdir(path)
        os.system("/Applications/INDELibleV1.03/src/indelible \n")
    #now read INDELible output
    file=open(path+"INDELible_sim_par"+str(i+1)+"_TRUE.phy")
    line=file.readline()
    for r in range(numReplicates+minimumRep):
        dataIndelible[i].append(([],[]))
        line=file.readline()
        seq=line.split()[1]
        k=0
        while k<len(seq):
            if seq[k]=="-":
                lengthDel=1
                k+=1
                while k<len(seq) and seq[k]=="-":
                    lengthDel+=1
                    k+=1
                dataIndelible[i][r][1].append(lengthDel)
            else:
                k+=1
        line=file.readline()
        seq=line.split()[1]
        k=0
        while k<len(seq):
            if seq[k]=="-":
                lengthIns=1
                k+=1
                while k<len(seq) and seq[k]=="-":
                    lengthIns+=1
                    k+=1
                dataIndelible[i][r][0].append(lengthIns)
            else:
                k+=1
        line=file.readline()
        line=file.readline()
        line=file.readline()

    file.close()

#Print values to screen
numInsP=[]
numInsI=[]
numDelP=[]
numDelI=[]
lenInsP=[]
lenInsI=[]
lenDelP=[]
lenDelI=[]
for i in range(len(parameters)):
    print("Parameter "+str(i+1))
    numInsP.append([])
    numInsI.append([])
    numDelP.append([])
    numDelI.append([])
    lenInsP.append([])
    lenInsI.append([])
    lenDelP.append([])
    lenDelI.append([])
    #print("phastSim values ")
    for r in range(numReplicates+minimumRep):
        numInsP[i].append(len(dataPhastSim[i][r][0]))
        numInsI[i].append(len(dataIndelible[i][r][0]))
        numDelP[i].append(len(dataPhastSim[i][r][1]))
        numDelI[i].append(len(dataIndelible[i][r][1]))
        lenInsP[i].append(mean(dataPhastSim[i][r][0]))
        lenInsI[i].append(mean(dataIndelible[i][r][0]))
        lenDelP[i].append(mean(dataPhastSim[i][r][1]))
        lenDelI[i].append(mean(dataIndelible[i][r][1]))
        #print("Insertions "+str(len(dataPhastSim[i][r][0]))+" "+str(mean(dataPhastSim[i][r][0])))
        #print("Deletions "+str(len(dataPhastSim[i][r][1]))+" "+str(mean(dataPhastSim[i][r][1])))

    #print("INDELible values ")
    #for r in range(numReplicates):
        #print("Insertions "+str(len(dataIndelible[i][r][0]))+" "+str(mean(dataIndelible[i][r][0])))
        #print("Deletions "+str(len(dataIndelible[i][r][1]))+" "+str(mean(dataIndelible[i][r][1])))
    print("Comparison Insertion Numbers:")
    print(stats.mannwhitneyu(numInsP[i], numInsI[i]))
    print("Mean phastSim: "+str(mean(numInsP[i]))+"   mean INDELible: "+str(mean(numInsI[i]))+"\n")
    print("Comparison Insertion Lengths:")
    print(stats.mannwhitneyu(lenInsP[i], lenInsI[i]))
    print("Mean phastSim: "+str(mean(lenInsP[i]))+"   mean INDELible: "+str(mean(lenInsI[i]))+"\n")
    print("Comparison Deletion Numbers:")
    print(stats.mannwhitneyu(numDelP[i], numDelI[i]))
    print("Mean phastSim: "+str(mean(numDelP[i]))+"   mean INDELible: "+str(mean(numDelI[i]))+"\n")
    print("Comparison Deletion Lengths:")
    print(stats.mannwhitneyu(lenDelP[i], lenDelI[i]))
    print("Mean phastSim: "+str(mean(lenDelP[i]))+"   mean INDELible: "+str(mean(lenDelI[i]))+"\n")


#make plot of file sizes
formats=["PhastSim","INDELible"]
colors=["green","blue"]
nLeaves=["Base","Longer branch","More insertions","More deletions","Longer insertions","Longer deletions"]
leavesRange=range(len(nLeaves))

title="Comparison of number of simulated insertions"
yAxisLabel="Number of simulated events"
values=[]
values.append([])
values.append([])
for n in leavesRange:
    values[0].append([])
    values[1].append([])
    for r in range(numReplicates+minimumRep):
        values[0][n].append(numInsP[n][r])
        values[1][n].append(numInsI[n][r])
plotValues=[]
for n in range(len(nLeaves)):
	plotValues.append([])
	for m in range(len(formats)):
		plotValues[n].append(values[m][n])
errplot(plotValues,formats,path+"plotNumInsetions.pdf",nLeaves,colors,'Simulation Scenario',None,title=title,yAxisLabel=yAxisLabel)

title="Comparison of simulated insertion lengths"
yAxisLabel="Length of insertions"
values=[]
values.append([])
values.append([])
for n in leavesRange:
    values[0].append([])
    values[1].append([])
    for r in range(numReplicates+minimumRep):
        values[0][n].append(lenInsP[n][r])
        values[1][n].append(lenInsI[n][r])
plotValues=[]
for n in range(len(nLeaves)):
	plotValues.append([])
	for m in range(len(formats)):
		plotValues[n].append(values[m][n])
errplot(plotValues,formats,path+"plotLenInsetions.pdf",nLeaves,colors,'Simulation Scenario',None,title=title,yAxisLabel=yAxisLabel)

title="Comparison of number of simulated deletions"
yAxisLabel="Number of simulated events"
values=[]
values.append([])
values.append([])
for n in leavesRange:
    values[0].append([])
    values[1].append([])
    for r in range(numReplicates+minimumRep):
        values[0][n].append(numDelP[n][r])
        values[1][n].append(numDelI[n][r])
plotValues=[]
for n in range(len(nLeaves)):
	plotValues.append([])
	for m in range(len(formats)):
		plotValues[n].append(values[m][n])
errplot(plotValues,formats,path+"plotNumDeletions.pdf",nLeaves,colors,'Simulation Scenario',None,title=title,yAxisLabel=yAxisLabel)

title="Comparison of simulated deletion lengths"
yAxisLabel="Length of deletions"
values=[]
values.append([])
values.append([])
for n in leavesRange:
    values[0].append([])
    values[1].append([])
    for r in range(numReplicates+minimumRep):
        values[0][n].append(lenDelP[n][r])
        values[1][n].append(lenDelI[n][r])
plotValues=[]
for n in range(len(nLeaves)):
	plotValues.append([])
	for m in range(len(formats)):
		plotValues[n].append(values[m][n])
errplot(plotValues,formats,path+"plotLenDeletions.pdf",nLeaves,colors,'Simulation Scenario',None,title=title,yAxisLabel=yAxisLabel)
