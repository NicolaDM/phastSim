#!/usr/bin/env bash

#BSUB -n 1          # number of tasks/CPUs in job
#BSUB -q standard   # queue
#BSUB -e /nfs/research/goldman/will/%J.err   # error file name in which %J is replaced by the job ID
#BSUB -o /nfs/research/goldman/will/%J.out  # output file name in which %J is replaced by the job ID
#BSUB -J will-phastSim-scale-expt
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
TIMES="1000 5000 100000"
SCALES="0.1 1 10"
CURRENT_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

for t in $TIMES; do
    for s in $SCALES; do
        if [ $t -gt 5001 ]
        then
            seqgenOptions=""
        else
            seqgenOptions="--seqgenSim"
        fi

        if [ $t -gt 5001 ]
        then
            indelibleOptions=""
        else
            indelibleOptions="--indelibleSim --indelibleSim2"
        fi     

        python $CURRENT_FOLDER/compareSimulators.py --nLeaves $t $seqgenOptions $indelibleOptions \
            --phastSim --seed 126 --replicates 10 --scale $s --noHierarchy \
            --seqgenPath $SEQ_GEN_PATH --indeliblePath $INDELIBLE_PATH \
            --reference $REF_PATH --path $OUTPATH --monitorMemory
    done
done
