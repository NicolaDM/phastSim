#!/usr/bin/env bash

#BSUB -n 1          # number of tasks/CPUs in job
#BSUB -q standard   # queue
#BSUB -e /nfs/research/goldman/will/%J.err   # error file name in which %J is replaced by the job ID
#BSUB -o /nfs/research/goldman/will/%J.out  # output file name in which %J is replaced by the job ID
#BSUB -J will-phastSim-codon-expt
#BSUB -M 4096

# Some preamble to run this script on the cluster, comment out these 2 lines if you are running on a laptop
export HOME="/nfs/research/goldman/will"
source $HOME/.bashrc

# YOU NEED TO SET THESE VARIABLES FOR YOUR OWN SETUP!
OUTPATH="/nfs/research/goldman/will/sim_output/"
INDELIBLE_PATH="/hps/software/users/goldman/will/INDELIBLEV1.03/bin/indelible"
SEQ_GEN_PATH="/hps/software/users/goldman/will/Seq-Gen-1.3.4/source/seq-gen"
REF_PATH="/hps/software/users/goldman/will/phastSim/phastSim/example/MN908947.3.fasta"


# Rest of script - no need to edit below here
TIMES="10 100 1000 20000"
CURRENT_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

for t in $TIMES; do
    if [ $t -gt 101 ]
    then
        indelibleOptions=""
    else
        indelibleOptions="--indelibleSim --indelibleSim2"
    fi    
    if [ $t -gt 1001 ]
    then
        seqgenOptions=""
    else
        seqgenOptions="--seqgenSim"
    fi

    # nucleotide model
    python $CURRENT_FOLDER/compareSimulators.py --nLeaves $t $seqgenOptions $indelibleOptions --phastSim  --seed 127 --replicates 10 \
        --seqgenPath $SEQ_GEN_PATH --indeliblePath $INDELIBLE_PATH \
        --reference $REF_PATH --path $OUTPATH --monitorMemory

    # nucleotide model + 10 cat
    python $CURRENT_FOLDER/compareSimulators.py --nLeaves $t $seqgenOptions $indelibleOptions --phastSim  --seed 127 --replicates 10 \
        --seqgenPath $SEQ_GEN_PATH --indeliblePath $INDELIBLE_PATH \
        --categoryRates 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 \
        --categoryProbs 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 \
        --reference $REF_PATH --path $OUTPATH --monitorMemory

    # nucleotide model + alpha
    python $CURRENT_FOLDER/compareSimulators.py --nLeaves $t $seqgenOptions $indelibleOptions --phastSim  --seed 127 --replicates 10 \
        --seqgenPath $SEQ_GEN_PATH --indeliblePath $INDELIBLE_PATH \
        --alpha 1.0 \
        --reference $REF_PATH --path $OUTPATH --monitorMemory

    # codon model
    python $CURRENT_FOLDER/compareSimulators.py --nLeaves $t $indelibleOptions --phastSim  --seed 127 --replicates 10 \
        --indeliblePath $INDELIBLE_PATH --codon \
        --reference $REF_PATH --path $OUTPATH --monitorMemory

    # codon model + 10 cat
    python $CURRENT_FOLDER/compareSimulators.py --nLeaves $t $indelibleOptions --phastSim  --seed 127 --replicates 10 \
        --indeliblePath $INDELIBLE_PATH --codon \
        --omegaCategoryRates 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 \
        --omegaCategoryProbs 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 \
        --reference $REF_PATH --path $OUTPATH --monitorMemory

    # codon model + alpha
    python $CURRENT_FOLDER/compareSimulators.py --nLeaves $t --phastSim  --seed 127 --replicates 10 \
        --codon --alpha 1.0 --omegaAlpha 1.0 \
        --reference $REF_PATH --path $OUTPATH --monitorMemory
done
