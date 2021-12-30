#!/usr/bin/env bash

#BSUB -n 1          # number of tasks/CPUs in job
#BSUB -q standard   # queue
#BSUB -e /nfs/research/goldman/will/%J.err   # error file name in which %J is replaced by the job ID
#BSUB -o /nfs/research/goldman/will/%J.out  # output file name in which %J is replaced by the job ID
#BSUB -J will-phastSim-indel-expt
#BSUB -M 2000

# Some preamble to run this script on the cluster, comment out these 2 lines if you are running on a laptop
export HOME="/nfs/research/goldman/will"
source $HOME/.bashrc

# YOU NEED TO SET THESE VARIABLES FOR YOUR OWN SETUP!
OUTPATH="/nfs/research/goldman/will/sim_output/"
INDELIBLE_PATH="/hps/software/users/goldman/will/INDELIBLEV1.03/bin/indelible"
SEQ_GEN_PATH="/hps/software/users/goldman/will/Seq-Gen-1.3.4/source/seq-gen"
REF_PATH="/hps/software/users/goldman/will/phastSim/phastSim/example/MN908947.3.fasta"


# Rest of script - no need to edit below here
CURRENT_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
TIMES="10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000"


for t in $TIMES; do
    echo "Working on experiment "$t
    if [ $t -gt 5000 ]
    then
        indelibleOptions=""
    else
        indelibleOptions="--indelibleSim --indelibleSim2"
    fi 

    echo $indelibleOptions
    python $CURRENT_FOLDER/compareSimulators.py --nLeaves $t --phastSim --useIndels --createFasta $indelibleOptions --seed 125 --replicates 10 \
        --indeliblePath $INDELIBLE_PATH \
        --reference $REF_PATH --path $OUTPATH --monitorMemory
done