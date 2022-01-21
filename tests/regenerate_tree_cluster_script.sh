#!/bin/bash

# I need this preamble to get my conda environment
# to work on the cluster for some reason
export HOME=/nfs/research/goldman/will
source $HOME/.bashrc

python /hps/software/users/goldman/will/phastSim/tests/regenerate_tree_cluster.py \
    --nSims 100 \
    --outputFolder $HOME/test_simple_nucleotide_sims
