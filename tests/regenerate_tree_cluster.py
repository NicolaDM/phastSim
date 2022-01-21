import os
from ete3 import Tree
import numpy as np
import argparse
from collections import Counter
from Bio import SeqIO


N_SIMS = 1
RANDOM_SEED = 5
EXPECTED_N_MUTATIONS_PER_BRANCH = 1000
N_LEAVES = 8
N_BRANCHES = 2 * N_LEAVES - 2
GENOME_LENGTH = 1000000
OUTPUT_FOLDER = "/nfs/research/goldman/will/test_output_1"
ROOT_GENOME_FREQUENCIES_STRING = "0.1 0.2 0.3 0.4"


def setup_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nSims", default=N_SIMS)
    parser.add_argument("--randomSeed", default=RANDOM_SEED)
    parser.add_argument("--mutationsPerBranch", default=EXPECTED_N_MUTATIONS_PER_BRANCH)
    parser.add_argument("--nLeaves", default=N_LEAVES)
    parser.add_argument("--genomeLength", default=GENOME_LENGTH)
    parser.add_argument("--outputFolder", default=OUTPUT_FOLDER)
    parser.add_argument("--rootGenomeFrequencies", default=ROOT_GENOME_FREQUENCIES_STRING)
    return parser

def load_args():

    parser = setup_args()
    args = parser.parse_args()

    global N_SIMS
    N_SIMS = int(args.nSims)

    global RANDOM_SEED
    RANDOM_SEED = int(args.randomSeed)
    np.random.seed(RANDOM_SEED)

    global EXPECTED_N_MUTATIONS_PER_BRANCH
    EXPECTED_N_MUTATIONS_PER_BRANCH = int(args.mutationsPerBranch)
    
    global N_LEAVES
    N_LEAVES = int(args.nLeaves)

    global N_BRANCHES
    N_BRANCHES = 2 * N_LEAVES - 2

    global GENOME_LENGTH
    GENOME_LENGTH = int(args.genomeLength)

    global OUTPUT_FOLDER
    OUTPUT_FOLDER =args.outputFolder

    global ROOT_GENOME_FREQUENCIES_STRING
    ROOT_GENOME_FREQUENCIES_STRING = args.rootGenomeFrequencies

def get_tree_length(t):
    
    return sum([x.dist + get_tree_length(x) for x in t.get_children()])
        
def _rescale_tree(t, scale_factor):

    for x in t.get_children():
        x.dist *= scale_factor
        _rescale_tree(x, scale_factor)

def rescale_tree(t, new_total_length):

    length = get_tree_length(t)
    scale_factor = float(new_total_length) / length
    _rescale_tree(t, scale_factor)

def get_raxml_rates(filepath):

    results = []

    with open(filepath) as f:
        for line in f:
            if len(line.split("Base frequencies:")) > 1:
                results = (line.split("Base frequencies: ")[-1])

            if len(line.split("rates[0]")) > 1:
                results = line.split("ac ag at cg ct gt: ")[-1] + " " + results

    return results


if __name__ == "__main__":

    load_args()

    summary_file = open(f"{OUTPUT_FOLDER}/output.csv", "w")
    summary_file.write("index, input_gtr_rates, output_gtr_rates, input_tree_length, output_tree_length, RF_distance, normalised_gtr_error_percent\n")

    # create N_sim many trees
    for i in range(N_SIMS):
        
        t = Tree()
        t.populate(N_LEAVES, random_branches=True, branch_range=(0.5, 1.5))
        print(get_tree_length(t))
        rescale_tree(t, float(EXPECTED_N_MUTATIONS_PER_BRANCH) * float(N_BRANCHES) / float(GENOME_LENGTH))
        print(get_tree_length(t))
        with open(f"{OUTPUT_FOLDER}/tree_{i}.tree", "w") as tree_out:
            tree_out.write(t.write())


    # create a root genome with given root frequencies, using an empty tree
    # this also does a phastSim simulation, but ignore it
    null_tree = Tree()
    null_tree.populate(0, names_library=["ref"])
    null_tree.dist = 0.0

    with open(f"{OUTPUT_FOLDER}/null_tree.tree", "w") as null_tree_file:
        null_tree_file.write(null_tree.write())

    os.system(
        f"""phastSim \
        --rootGenomeLength {GENOME_LENGTH} \
        --rootGenomeFrequencies {ROOT_GENOME_FREQUENCIES_STRING} \
        --treeFile {OUTPUT_FOLDER}/null_tree.tree \
        --outpath {OUTPUT_FOLDER}/ \
        --outputFile my_ref \
        --createFasta \
        --seed {np.random.randint(1000000000)}
        """)

    reference = SeqIO.read(f"{OUTPUT_FOLDER}/my_ref.fasta", format="fasta")
    observed_frequencies = Counter(reference.seq)
    OBSERVED_ROOT_GENOME_FREQUENCIES_STRING = " ".join([str(float(observed_frequencies[k])/ GENOME_LENGTH) for k in "ACGT"])


    for i in range(N_SIMS):

        # create some random GTR rates
        rates = np.random.uniform(size=6)
        rate_gt = rates[-1]
        gtr_rates = " ".join([str(x) for x in rates]) + " " + OBSERVED_ROOT_GENOME_FREQUENCIES_STRING

        # need this version formatted in the same way as raxml does their rates so they can be compared
        gtr_rates_string_formatted = " ".join([str(x/rate_gt) for x in rates]) + " " + OBSERVED_ROOT_GENOME_FREQUENCIES_STRING

        # do the phastSim simulations
        os.system(
            f"""phastSim \
            --reference {OUTPUT_FOLDER}/my_ref.fasta \
            --treeFile {OUTPUT_FOLDER}/tree_{i}.tree \
            --mutationRates GTR {gtr_rates} \
            --outpath {OUTPUT_FOLDER}/ \
            --outputFile phastSim_{i} \
            --createFasta \
            --seed {np.random.randint(1000000000)}
            """
        )

        # try to regenerate the tree with RAxML
        os.system(f"""raxmlHPC \
            -m GTRCAT -V -F -c 1 \
            -n rax_{i} \
            -s phastSim_{i}.fasta \
            -p {np.random.randint(1000000000)} \
            -w {OUTPUT_FOLDER}/
        """)

        raxml_estimated_rates = get_raxml_rates(f"{OUTPUT_FOLDER}/RAxML_info.rax_{i}")

        # put the stuff that we want into an output file
        # that is:
        # i, input_gtr_rates, estimated_gtr_rates, input_tree_length, output_tree_length, RF_distance
        input_tree = Tree(newick=f"{OUTPUT_FOLDER}/tree_{i}.tree")
        output_tree = Tree(f"{OUTPUT_FOLDER}/RAxML_result.rax_{i}")
        rf_dist = input_tree.robinson_foulds(output_tree, unrooted_trees=True)[0]

        gtr_rates_formatted = np.array([float(x) for x in gtr_rates_string_formatted.split(" ")])
        raxml_rates = np.array([float(x) for x in raxml_estimated_rates.split(" ")])
        
        normalised_gtr_error_pc = 100.0 * np.sqrt(
            np.sum((gtr_rates_formatted - raxml_rates) ** 2) / np.sum(gtr_rates_formatted ** 2)
        )

        summary_file.write(f"{i}, {gtr_rates_string_formatted}, {raxml_estimated_rates}, {get_tree_length(input_tree)}, {get_tree_length(output_tree)}, {rf_dist}, {normalised_gtr_error_pc}\n")
