#!/usr/bin/env bash

#BSUB -n 1          # number of tasks/CPUs in job
#BSUB -q standard   # queue
#BSUB -e /nfs/research/goldman/will/%J.err   # error file name in which %J is replaced by the job ID
#BSUB -o /nfs/research/goldman/will/%J.out  # output file name in which %J is replaced by the job ID
#BSUB -J will-phastSim-general-expt
#BSUB -M 4096

# Some preamble to run this script on the cluster, comment out these 2 lines if you are running on a laptop
export HOME="/nfs/research/goldman/will"
source $HOME/.bashrc

# YOU NEED TO SET THESE VARIABLES FOR YOUR OWN SETUP!
OUTPATH="/nfs/research/goldman/will/sim_output/"
INDELIBLE_PATH="/hps/software/users/goldman/will/INDELibleV1.03/bin/indelible"
SEQ_GEN_PATH="/hps/software/users/goldman/will/Seq-Gen-1.3.4/source/seq-gen"
REF_PATH="/hps/software/users/goldman/will/phastSim/phastSim/example/MN908947.3.fasta"


# Rest of script - no need to edit below here
TIMES="10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000"

# I just want to get the directory of this script and apparently this stupid code is the best way to do it
CURRENT_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PYVOLVE_SCRIPT_PATH=$CURRENT_FOLDER/runPyvolve.py


for t in $TIMES; do
    echo "Working on experiment "$t
    if [ $t -gt 51 ]
    then
        pyvolveOptions=""
    else
        pyvolveOptions="--pyvolveSim"
    fi 
    echo $pyvolveOptions
    
    if [ $t -gt 10001 ]
    then
        indelibleOptions=""
    else
        indelibleOptions="--indelibleSim --indelibleSim2"
    fi 
    echo $indelibleOptions
    
    if [ $t -gt 20001 ]
    then
        seqgenOptions=""
    else
        seqgenOptions="--seqgenSim"
    fi         
    echo $seqgenOptions
    
    if [ $t -gt 50001 ]
    then
        phastSimOptions="--phastSim"
    else
        phastSimOptions="--phastSim --createFasta"
    fi 
    echo $phastSimOptions
    python $CURRENT_FOLDER/compareSimulators.py --nLeaves $t $pyvolveOptions $seqgenOptions $phastSimOptions $indelibleOptions --seed 123 --replicates 10 \
        --pyvolvePath $PYVOLVE_SCRIPT_PATH --indeliblePath $INDELIBLE_PATH --seqgenPath $SEQ_GEN_PATH \
        --reference $REF_PATH --path $OUTPATH --monitorMemory
done