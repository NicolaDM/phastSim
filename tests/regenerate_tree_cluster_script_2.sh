#!/bin/bash

# I need this preamble to get my conda environment
# to work on the cluster for some reason
export HOME=/nfs/research/goldman/will
source $HOME/.bashrc

module add raxml-8.2.11-gcc-9.3.0-mjwrm3x

python /hps/software/users/goldman/will/phastSim/tests/regenerate_tree_cluster.py \
    --nSims 100 \
    --genomeLength 1000 \
    --randomSeed 199 \
    --nLeaves 128 \
    --mutationsPerBranch 100 \
    --raxmlModelString +-m+GTRCAT+-c+3+--JC69 \
    --rootGenomeFrequencies 0.25+0.25+0.25+0.25 \
    --phastSimOptions +--categoryRates+0.25+0.5+0.25+categoryProbs+0.1+1.0+10.0+--mutationRates+JC69 \
    --outputFolder $HOME/test_nucleotide_sims_categories
