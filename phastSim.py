import numpy as np
import argparse
from ete3 import Tree
import time


# CONSTANTS
# possible alleles
class Constants:
    def __init__(self):
        self.alleles = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3, "u": 3, "U": 3}
        self.allelesList = ["A", "C", "G", "T"]
        self.nAlleles = 4


def setup_argument_parser():
    parser = argparse.ArgumentParser(
        description='Efficiently simulate sequence evolution along phylogenies with short branches.')
    parser.add_argument('--path', default="", help='Path where to run simulations.')
    parser.add_argument('--reference', default="MN908947.3.fasta",
                        help='File containing the reference genome to be used as root genome. '
                             'To be found in the folder specified with --path.')
    parser.add_argument("--rootGenomeLength",
                        help="A root genome of this size is randomly created. default=0 means that the root genome "
                             "is instead taken as in the reference file.",
                        type=int, default=0)
    parser.add_argument("--rootGenomeFrequencies",
                        help="Frequencies of different states (non-stop codons or nucleotides depending "
                             "on the simulations). If not provided, will use those for SARS-CoV-2 reference genome.",
                        type=float, nargs='+', default=[1.0])
    parser.add_argument('--treeFile', default="exampleTree.tree",
                        help='Name of file containing the tree used to simulate sequences (assumed within the '
                             '--path and in newick format).')
    parser.add_argument('--scale', default=1.0, type=float,
                        help='Scale the simulation tree by this amount (default 1.0). Branch lengths are assumed'
                             ' in terms of expected substitutions per site (more or less, as frequencies changes through'
                             ' time, total mutation rate might also  change).')
    parser.add_argument("--seed", help="Seed for random simulator in genes simulations", type=int, default=1)
    parser.add_argument("--alpha",
                        help="Parameter of the gamma distribution for mutation rate variation; each site will "
                             "then have a separate rate. If specified, continuous rate variation is assumed, otherwise"
                             " homogeneous rates are used unless the --categoryRates option is used.",
                        type=float, default=0.0)
    parser.add_argument("--invariable",
                        help="Proportion of invariable sites, that is, sites that have a mutation rate of 0.0 .",
                        type=float, default=0.0)
    parser.add_argument("--mutationRates",
                        help="Mutation rates, by default using the neutral rates estimated from SARS-CoV-2; so far"
                             " only exactly 12 input values allowed (r_AC, r_AG, r_AT, r_CA, etc...) corresponding to"
                             " an UNREST nucleotide substitution model.",
                        type=float, nargs='+',
                        default=[0.039, 0.310, 0.123, 0.140, 0.022, 3.028, 0.747, 0.113, 2.953, 0.056, 0.261, 0.036])
    parser.add_argument("--codon",
                        help="Run simulations under a codon model, where mutation rates are used to describe nucleotide"
                             " mutation rates, and omegas are used to describe the effect of selection at the amino acid"
                             " level. Default is false (uses a nucleotide substitution model).",
                        action="store_true")
    parser.add_argument("--omegaAlpha",
                        help="Parameter of the gamma distribution for omega variation; each codon will then have a "
                             "separate omega. If specified, continuous omega variation is assumed, otherwise homogeneous"
                             " omegas are used unless the --omegaCategoryRates option is used.",
                        type=float, default=0.0)
    parser.add_argument("--omegaCategoryProbs",
                        help="Probabilities of omega categories. They are supposed to sum up to one, but if they don't"
                             " they are normalized to do so. By default only 1 category is simulated.",
                        type=float, nargs='+', default=[1.0])
    parser.add_argument("--omegaCategoryRates",
                        help="Omegas of different omega categories. The overall evolutionary rate is renormalized so "
                             "that the expected number of substitutions per branch length unit per site is 1."
                             " By default only 1 category of rate 1.0 is simulated. The number of omegas has to be the "
                             "same as the number of omega category probabilities, otherwise an error is thrown.",
                        type=float, nargs='+', default=[1.0])
    parser.add_argument("--categoryProbs",
                        help="Probabilities of rate categories. They are supposed to sum up to one, but if they don't "
                             "they are normalized to do so. By default only 1 category is simulated.",
                        type=float, nargs='+', default=[1.0])
    parser.add_argument("--categoryRates",
                        help="Rates of site categories. The overall mutation rate is renormalized so that the expected "
                             "number of substitutions per branch length unit per site is 1. By default only 1 category"
                             " of rate 1.0 is simulated. The number of rates has to be the same as the number of "
                             "category probabilities, otherwise an error is thrown.",
                        type=float, nargs='+', default=[1.0])
    parser.add_argument("--hyperMutProbs",
                        help="Proportions of sites with frequently recurring mutations. They are supposed to be rare,"
                             " and not sum up to one. By default no recurring mutations are simulated.",
                        type=float, nargs='+', default=[])
    parser.add_argument("--hyperMutRates",
                        help="Rates of recurring mutations. Different classes of recurring mutations with different "
                             "hypermutabilities can be specified. By defaultno recurrent mutations are simulated. "
                             "The number of rates has to be the same as the number of recurring mutation site "
                             "proportions, otherwise an error is thrown.",
                        type=float, nargs='+', default=[])
    parser.add_argument("--noHierarchy",
                        help="Run without hierarchical algorithm; the latter is faster with more complex models.",
                        action="store_true")
    parser.add_argument("--verbose", help="Turns on verbose mode.", action="store_true")
    parser.add_argument('--outputFile', default="sars-cov-2_simulation_output",
                        help='Output file name containing the simulated genomes in succint format. The file will be '
                             'created within the folder specified with --path.')
    parser.add_argument("--createNewick",
                        help="Create a newick file annotated with the simulated mutation events "
                             "(default name sars-cov-2_simulation_output.tree).",
                        action="store_true")
    parser.add_argument("--createFasta",
                        help="Create a fasta file with the simulated genomes "
                             "(default name sars-cov-2_simulation_output.fasta).",
                        action="store_true")
    parser.add_argument("--createPhylip",
                        help="Create a phylip file with the simulated genomes "
                             "(default name sars-cov-2_simulation_output.phy).",
                        action="store_true")
    return parser




class phastSim:
    def __init__(self, args):
        self.args = args
        self.const = Constants()

