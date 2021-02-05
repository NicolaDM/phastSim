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




class phastSim_run:
    def __init__(self, args):
        self.args = args
        self.const = Constants()
        self.hierarchy = not self.args.noHierarchy


    def init_rootGenome(self):
        # initialise a genome either from a file or creating a new one
        if self.args.rootGenomeLength == 0:
            ref, refList = self.load_rootGenome_file()

        else:
            if self.args.codon:
                ref, refList = self.create_rootGenome_codon()
            else:
                ref, refList = self.create_rootGenome_nuc()

        if self.args.codon and (len(ref) % 3) != 0:
            print(
                "Warning: when simulating under a codon model, the ancestral genome length has to be a multiple of 3.")
            print("I will remove the last few bases of the reference and assume the rest is made of coding sequence.")
            ref = ref[:len(ref) - (len(ref) % 3)]

        return ref, refList


    def load_rootGenome_file(self):
        # collect reference
        reference_file = f'{self.args.path}/{self.args.reference}'
        file = open(reference_file)
        line = file.readline()
        ref = ""
        while line != "":
            line = file.readline()
            ref += line.replace("\n", "")
        file.close()
        print(f"\n Finished reading reference genome at {reference_file} with {str(len(ref))} bases.")
        refList = list(ref)
        return ref, refList


    def create_rootGenome_codon(self):
        rootGenomeFrequencies = self.args.rootGenomeFrequencies
        rootGenomeLength = self.args.rootGenomeLength

        if len(rootGenomeFrequencies) < 2:
            print("Using codon frequencies from the SARS-CoV-2 genome to define root genome:")
            rootGenomeFrequencies = [0.03787801, 0.01775209, 0.02012592, 0.03694912, 0.03034369, 0.00701827,
               0.00371555, 0.03292393, 0.01610073, 0.00412839, 0.00485086, 0.01630715,
               0.01620394, 0.00970172, 0.02012592, 0.02652493, 0.02611209, 0.00588296,
               0.01145629, 0.01341728, 0.01620394, 0.00299308, 0.00175457, 0.01971308,
               0.00175457, 0.00350913, 0.0010321, 0.00877284, 0.01042419, 0.0094953,
               0.00464444, 0.02755702, 0.0326143, 0.0188874, 0.01269481, 0.0336464,
               0.01857777, 0.00959851, 0.00258025, 0.03694912, 0.01228197, 0.01063061,
               0.00175457, 0.03478171, 0.01826814, 0.01135308, 0.0115595, 0.03932294, 0.0,
               0.01795851, 0.0, 0.02817628, 0.01878419, 0.0052637, 0.00123852, 0.02229332,
               0.0, 0.00681185, 0.01135308, 0.02353184, 0.02580246, 0.01506863, 0.01692641, 0.03591702]
            print(rootGenomeFrequencies)

        elif len(rootGenomeFrequencies) != 61 and len(rootGenomeFrequencies) != 64:
            print("Error, wrong number of root frequencies given")
            exit()

        if len(rootGenomeFrequencies) == 61:
            rootGenomeFrequencies = rootGenomeFrequencies[:48] + [0.0] + [rootGenomeFrequencies[48]] + [
                0.0] + rootGenomeFrequencies[49:54] + [0.0] + rootGenomeFrequencies[54:]


        sum_freq = np.sum(rootGenomeFrequencies)

        for i in range(64):
            rootGenomeFrequencies[i] = rootGenomeFrequencies[i] / sum_freq
        if sum_freq > 1.000001 or sum_freq < 0.999999:
            print("\n Normalizing root state frequencies. New frequencies:")
            print(rootGenomeFrequencies)


        if (rootGenomeLength % 3) != 0:
            print("Codon model, but root genome length not multiple of 3. Simulating " + str(
                int(rootGenomeLength / 3)) + " codons.")
        cods = np.random.choice(np.arange(64), size=int(rootGenomeLength / 3), replace=True,
                                p=rootGenomeFrequencies)

        allelesList = self.const.allelesList
        codonAllelesList = []
        for i1 in range(4):
            for i2 in range(4):
                for i3 in range(4):
                    codonAllelesList.append(allelesList[i1] + allelesList[i2] + allelesList[i3])
        codList = []
        for i in range(int(rootGenomeLength / 3)):
            codList.append(codonAllelesList[cods[i]])
        ref = "".join(codList)
        refList = list(ref)
        return ref, refList


    def create_rootGenome_nuc(self):
        rootGenomeFrequencies = self.args.rootGenomeFrequencies
        rootGenomeLength = self.args.rootGenomeLength

        if len(rootGenomeFrequencies) < 2:
            print("Using nucleotide frequencies from the SARS-CoV-2 genome to define root genome:")
            rootGenomeFrequencies = [0.29943483931378123, 0.18366050229074005, 0.19606728421897468,
                                     0.32083737417650404]
            print(rootGenomeFrequencies)
        elif len(rootGenomeFrequencies) != 4:
            print("Error, wrong number of root frequencies given")
            exit()

        sum_freq = np.sum(rootGenomeFrequencies)

        for i in range(4):
            rootGenomeFrequencies[i] = rootGenomeFrequencies[i] / sum_freq
        if sum_freq > 1.000001 or sum_freq < 0.999999:
            print("\n Normalizing root state frequencies. New frequencies:")
            print(rootGenomeFrequencies)

        allelesList = self.const.allelesList
        refList = np.random.choice(allelesList, size=int(rootGenomeLength), replace=True, p=rootGenomeFrequencies)
        ref = "".join(refList)
        return ref, refList


    def init_substitution_rates(self):
        # substitution rates
        mutationRates = self.args.mutationRates
        if len(mutationRates) == 12:
            print("\n Assuming UNREST nucleotide mutation model.")
            mutMatrix = np.zeros((4, 4), dtype=float)
            index = 0
            for i in range(4):
                sum = 0.0
                for j in range(4):
                    if j != i:
                        mutMatrix[i][j] = mutationRates[index]
                        sum += mutationRates[index]
                        index += 1
                mutMatrix[i][i] = -sum
        else:
            print("\n Number of mutation rates " + str(len(mutationRates)) + ", model not implemented yet.")
            print(mutationRates)
            exit()
        print("\n Mutation rate matrix:")
        print(mutMatrix)
        return mutMatrix


    def init_gamma_rates(self, ref):
        categoryProbs = self.args.categoryProbs
        categoryRates = self.args.categoryRates

        if self.args.alpha >= 0.000000001:
            print(f"Using a continuous gamma rate distribution with parameter alpha={self.args.alpha}")
            if not self.hierarchy:
                print("Error, continuous rate model only allowed with hierarchical approach")
                exit()
            gammaRates = np.random.gamma(self.args.alpha, 1.0 / self.args.alpha, size=len(ref))

        else:
            nCat = len(categoryProbs)
            nRates = len(categoryRates)
            sum_probs = np.sum(categoryProbs)
            for i in range(nCat):
                categoryProbs[i] = categoryProbs[i] / sum_probs

            if sum_probs > 1.000001 or sum_probs < 0.999999:
                print("\n Normalizing probabilities of site categories. New probabilities:")
                print(categoryProbs)

            if nCat != nRates:
                print(f"Issue with number of category probs {nCat} and number of category rates {nRates}")
                exit()
            print("Using a discrete distribution for variation in rates across the genome.")
            print(categoryProbs)
            print(categoryRates)
            # sample category for each site of the genome
            gammaRates = np.zeros(len(ref))
            categories = np.random.choice(nCat, size=len(ref), p=categoryProbs)
            for i in range(len(ref)):
                gammaRates[i] = categoryRates[categories[i]]

        invariable = self.args.invariable
        if invariable >= 0.000000001:
            print(f"Proportion of invariable {invariable}")
            categoriesInv = np.random.choice(2, size=len(ref), p=[1.0 - invariable, invariable])
            gammaRates[np.nonzero(categoriesInv)[0]] = 0.0

        return gammaRates

