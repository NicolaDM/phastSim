#!/bin/bash

# I need this preamble to get my conda environment
# to work on the cluster for some reason
export HOME=/nfs/research/goldman/will
source $HOME/.bashrc

module add raxml-8.2.11-gcc-9.3.0-mjwrm3x

python /hps/software/users/goldman/will/phastSim/tests/regenerate_tree_cluster.py \
    --nSims 1 \
    --genomeLength 1000 \
    --randomSeed 198 \
    --nLeaves 128 \
    --mutationsPerBranch 100 \
    --raxmlModelString +-m+GTRCATI+-V+--JC69 \
    --rootGenomeFrequencies 0.25+0.25+0.25+0.25 \
    --phastSimOptions +--invariable+0.05+--mutationRates+JC69 \
    --outputFolder $HOME/test_nucleotide_sims_invariable
