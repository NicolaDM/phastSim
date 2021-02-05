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
        self.stopCodons = ["TAA", "TGA", "TAG"]


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

        self.ref_len = len(ref)
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


    def init_gamma_rates(self):
        categoryProbs = self.args.categoryProbs
        categoryRates = self.args.categoryRates

        if self.args.alpha >= 0.000000001:
            print(f"Using a continuous gamma rate distribution with parameter alpha={self.args.alpha}")
            if not self.hierarchy:
                print("Error, continuous rate model only allowed with hierarchical approach")
                exit()
            gammaRates = np.random.gamma(self.args.alpha, 1.0 / self.args.alpha, size=self.ref_len)

        else:
            nCat = len(categoryProbs)
            self.nCat = nCat
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
            gammaRates = np.zeros(self.ref_len)
            categories = np.random.choice(nCat, size=self.ref_len, p=categoryProbs)
            self.categories = categories
            for i in range(self.ref_len):
                gammaRates[i] = categoryRates[categories[i]]

        invariable = self.args.invariable
        if invariable >= 0.000000001:
            print(f"Proportion of invariable {invariable}")
            categoriesInv = np.random.choice(2, size=self.ref_len, p=[1.0 - invariable, invariable])
            gammaRates[np.nonzero(categoriesInv)[0]] = 0.0

        return gammaRates


    def init_hypermutation_rates(self):
        hyperMutProbs = self.args.hyperMutProbs
        hyperMutRates = self.args.hyperMutRates

        # dealing with hypermutation rates
        nProbs = len(hyperMutProbs)
        nRates = len(hyperMutRates)

        if nProbs != nRates:
            print(f"Issue with number of hypermutation category probs {nProbs} and number of hypermutation"
                  f" category rates {nRates}")
            exit()
        for i in hyperMutRates:
            if i <= 1.0:
                print("It doesn't make sense to have hypermutability class with mutability <=1.0 . hyperMutRates:")
                print(hyperMutRates)
                exit()

        sumHyper = np.sum(hyperMutProbs)
        if sumHyper > 0.1:
            print("WARNING: hypermutable sites are supposed to be rare, but total proportion is " + str(sum))
        if sumHyper > 1.0:
            exit()
        newHyperMutProbs = [1.0 - sumHyper] + hyperMutProbs

        # sample hypermutability for each site of the genome.
        hyperCategories = np.random.choice(nProbs + 1, size=self.ref_len, p=newHyperMutProbs)
        print("Hypermutation class probabilities:")
        print(newHyperMutProbs)
        print("Hypermutation class rates:")
        print(hyperMutRates)
        return hyperCategories


    def init_codon_substitution_model(self):
        # define codon substitution model
        if not self.hierarchy:
            print("Error: codon model only allowed with hierarchical model.")
            exit()

        omegaAlpha = self.args.omegaAlpha
        omegaCategoryProbs = self.args.omegaCategoryProbs
        omegaCategoryRates = self.args.omegaCategoryRates

        nCodons = int(self.ref_len / 3)
        self.nCodons = nCodons
        print("Using a codon model")

        if omegaAlpha >= 0.000000001:
            print("Using a continuous gamma distribution with parameter alpha=" + str(
                omegaAlpha) + " for variation in omega across codons.")
            omegas = np.random.gamma(omegaAlpha, 1.0 / omegaAlpha, size=nCodons)
        else:
            nCatOmega = len(omegaCategoryProbs)
            nRateOmega = len(omegaCategoryRates)

            sum_probs = np.sum(omegaCategoryProbs)
            for i in range(nCatOmega):
                omegaCategoryProbs[i] = omegaCategoryProbs[i] / sum_probs
            if sum_probs > 1.000001 or sum_probs < 0.999999:
                print("\n Normalizing probabilities of omega categories. New probabilities:")
                print(omegaCategoryProbs)

            if nCatOmega != nRateOmega:
                print(f"Issue with number of omega category probs {nCatOmega} and number of "
                      f"omega category rates {nRateOmega}")
                exit()
            print("Using a discrete distribution for variation in omega across codons.")
            print(omegaCategoryProbs)
            print(omegaCategoryRates)

            # sample category for each site of the genome
            omegas = np.zeros(nCodons)
            omegaCategories = np.random.choice(nCatOmega, size=nCodons, p=omegaCategoryProbs)
            for i in range(nCodons):
                omegas[i] = omegaCategoryRates[omegaCategories[i]]

        return omegas



class GenomeTree:
    def __init__(self, nCodons, codon, ref, gammaRates, omegas, mutMatrix, hyperCategories, hyperMutRates, file):
        self.codon = codon
        self.ref = ref
        self.gammaRates = gammaRates
        self.omegas = omegas
        self.mutMatrix = mutMatrix

        self.hyperCategories = hyperCategories
        self.hyperMutRates = hyperMutRates
        self.file = file

        const = Constants()
        self.stopCodons = const.stopCodons
        self.alleles = const.alleles
        self.allelesList = const.allelesList
        self.nAlleles = const.nAlleles

        if not codon:
            self.nTerminalNodes = len(ref)
            # we use numpy to efficiently create many copies of the initial mutation matrix
            self.mutMatrixRepeats = np.repeat(mutMatrix[np.newaxis, ...], len(ref), axis=0)

        else:
            self.nTerminalNodes = nCodons
            # pre-calculating relationships between representations of codons as int and as triplets.
            translationList, codonIndices, codonIndices2, codonAlleles, codonAllelesList = codon_translation_list(
                const.allelesList)

            self.codonIndices = codonIndices
            self.codonIndices2 = codonIndices2
            self.codonAlleles = codonAlleles
            self.codonAllelesList = codonAllelesList

            # Creating quick look-up table that tells you if a mutation is synonymous or not
            isNonsynom, isIntoStop = codon_lookup_table(translationList, codonIndices, codonIndices2)
            self.isNonsynom = isNonsynom
            self.isIntoStop = isIntoStop

        # root of the starting (level 0) genome tree
        self.genomeRoot = genomeNode()
        # position of the genome at which the corresponding node refers to. For the root, it's the whole genome.
        self.genomeRoot.genomePos = [0, self.nTerminalNodes - 1]


    def populateGenomeTree(self, node):
        if node.genomePos[0] == node.genomePos[1]:
            # these are terminal nodes. Only these nodes need information regarding specific alleles and rates.
            pos = node.genomePos[1]
            node.isTerminal = True
            # the corresponding node in the level 0 tree; in this case, it's itself.
            # Only terminal nodes in the level 0 tree need to contain info regarding all mutation rates.
            node.refNode = node
            # the current allele at this terminal node. It starts as the reference allele at the considered position.
            if self.codon:
                if (pos < self.nTerminalNodes - 1) and (self.ref[pos * 3:(pos + 1) * 3] in self.stopCodons):
                    print(f"Warning: stop codon in the middle of the reference {pos * 3 + 1}"
                          f" {self.ref[pos * 3:(pos + 1) * 3]}."
                          f"Continuing simulations but setting omega=0 for this position.")
                    self.omegas[pos] = 0.0
                # print(ref[pos*3:(pos+1)*3])
                node.allele = self.codonAlleles[self.ref[pos * 3:(pos + 1) * 3]]
                # print(node.allele)
                self.file.write(str(pos * 3 + 1) + "-" + str(pos * 3 + 3) + "\t" + str(self.omegas[pos]) + "\t")
                node.rates = {}
                node.hyper = {}
                node.rates[node.allele] = np.zeros(10)
                indeces = self.codonIndices[node.allele]
                node.rate = 0.0
                for i2 in range(3):
                    pos2 = pos * 3 + i2
                    self.file.write(str(self.gammaRates[pos2]) + "\t" + str(self.hyperCategories[pos2]) + "\t")
                    nuc1 = indeces[i2]
                    for i3 in range(3):
                        nuc2 = (nuc1 + i3 + 1) % 4
                        if self.isNonsynom[node.allele, i2, i3]:
                            if self.isIntoStop[node.allele, i2, i3]:
                                node.rates[node.allele][i2 * 3 + i3] = 0.0
                            else:
                                node.rates[node.allele][i2 * 3 + i3] = self.omegas[pos] * self.mutMatrix[nuc1][nuc2] * \
                                                                       self.gammaRates[pos2]
                        else:
                            node.rates[node.allele][i2 * 3 + i3] = self.mutMatrix[nuc1][nuc2] * self.gammaRates[pos2]
                    if self.hyperCategories[pos2] > 0:
                        # site with hypermutation. Sample a random allele i to hypermutate,
                        # and a random destination allele j;
                        # the mutation rate from i to j at the given position is then enhanced.
                        i = np.random.choice(self.nAlleles)
                        j = np.random.choice(self.nAlleles - 1)
                        node.hyper[i2] = (i, j)
                        if i == nuc1:
                            node.rates[node.allele][i2 * 3 + j] *= self.hyperMutRates[self.hyperCategories[pos2] - 1]
                        self.file.write(self.allelesList[i] + "\t" + self.allelesList[(i + j + 1) % 4] + "\t")
                    else:
                        self.file.write(".\t" + ".\t")
                    for i3 in range(3):
                        node.rate += node.rates[node.allele][i2 * 3 + i3]
                    node.rates[node.allele][9] = node.rate
                self.file.write("\n")

            else:
                node.allele = self.alleles[self.ref[pos]]
                # evolutionary categories
                nodeHyper = self.hyperCategories[pos]
                self.file.write(str(pos + 1) + "\t" + str(self.gammaRates[pos]) + "\t" + str(nodeHyper) + "\t")
                # the mutation rates for the considered site
                node.rates = self.mutMatrixRepeats[pos]
                node.rates *= self.gammaRates[pos]
                # node.rates=np.copy(mutMatrix)*categoryRates[node.category]
                # NICOLA: THIS COULD MAYBE BE MADE MORE EFFICIENT.
                # FOR EXAMPLE, ONE COULD INITIALIZE RATES ONLY WHEN SAMPLING AT A LEAF,
                # AND OTHERWISE INITIALIZE ONLY THE TOTAL RATE.
                # ALSO, INSTATED OF COPYING AND STORING MULTIPLE MUTATION MATRICES, ONE COULD LINK TO THE SAME ONE
                # FOR SITES WITH THE SAME RATES.
                # now sample which alleles are affected by hypermutation and store info in positions and extra vectors
                if nodeHyper > 0:
                    # site with hypermutation. Sample a random allele i to hypermutate,
                    # and a random destination allele j;
                    # the mutation rate from i to j at the given position is then enhanced.
                    i = np.random.choice(self.nAlleles)
                    j = np.random.choice(self.nAlleles - 1)
                    j = (i + j + 1) % self.nAlleles
                    node.rates[i][i] -= node.rates[i][j] * (self.hyperMutRates[nodeHyper - 1] - 1.0)
                    node.rates[i][j] *= self.hyperMutRates[nodeHyper - 1]
                    self.file.write(self.allelesList[i] + "\t" + self.allelesList[j] + "\n")
                else:
                    self.file.write(".\t" + ".\n")
                # total mutation rate at the node
                node.rate -= node.rates[node.allele][node.allele]

        else:
            # split the considered part of the genome in two, assign each half to both children,
            # then call the populate function iteratively on the children.
            middle = node.genomePos[0] + int((node.genomePos[1] - node.genomePos[0]) / 2)
            firstHalf = [node.genomePos[0], middle]
            secondHalf = [middle + 1, node.genomePos[1]]
            firstChild = genomeNode()  # upNode=node
            firstChild.genomePos = firstHalf
            secondChild = genomeNode()  # upNode=node
            secondChild.genomePos = secondHalf
            node.belowNodes = [firstChild, secondChild]
            self.populateGenomeTree(node=firstChild)
            self.populateGenomeTree(node=secondChild)
            node.rate = firstChild.rate + secondChild.rate


    def normalize_rates(self, scale):
        # I am assuming the branch lengths are in number of substitutions per nucleotide,
        # even though we might be simulating a codon model.
        norm = self.genomeRoot.rate / len(self.ref)

        print("\n Total cumulative mutation rate per site before normalization: " + str(norm))
        # We rescale by the input normalization factor, this is the same as rescaling all
        # the branch lengths by this rescaling factor
        norm = norm / scale

        # define a function either for codon or nucleotide mode
        if self.codon:
            # function to iteratively normalize all rates
            def normalize(node, norm):
                node.rate /= norm
                if node.isTerminal:
                    node.rates[node.allele] /= norm
                else:
                    normalize(node.belowNodes[0], norm)
                    normalize(node.belowNodes[1], norm)

        else:
            # function to iteratively normalize all rates
            def normalize(node, norm):
                node.rate /= norm
                if node.isTerminal:
                    node.rates /= norm
                else:
                    normalize(node.belowNodes[0], norm)
                    normalize(node.belowNodes[1], norm)

        # normalize the rates
        normalize(node=self.genomeRoot, norm=norm)
        # and return the norm factor
        return norm






# node representing an element of the genome hierarchy, summarizing the total rate of the nodes below it,
# that is, the total rate of part of the genome.
class genomeNode:
    def __init__(self, level=0):  # , upNode=None
        # total mutation rate at this genome node
        self.rate = 0.0
        # the parent node in the genome tree in the same layer.
        # self.upNode = upNode
        # terminal nodes correspond to individual sites of the genome
        # (or individual evolutionary units for whatever they are).
        self.isTerminal = False
        # level of the genome tree hierarchy. higher levels are below, and are not linked by nodes in layers above.
        self.level = level






def check_start_stop_codons(ref, gammaRates, omegas):
    # if first codon is a start codon, or last is a stop codon, don't allow them to evolve
    # if ref[0:3]=="ATG" or ref[0:3]=="atg" or ref[0:3]=="AUG" or ref[0:3]=="aug":
    if ref[0:3] in ["ATG", "atg", "AUG", "aug"]:
        gammaRates[0] = 0.0
        gammaRates[1] = 0.0
        gammaRates[2] = 0.0
    # stopCodons=["TAA","TGA","TAG","tag","tga","taa","UAG","UAA","UGA","uaa","uag","uga"]
    stopCodons = ["TAA", "TGA", "TAG"]
    if ref[-2:] in stopCodons:
        omegas[-1] = 0.0

    return gammaRates, omegas


def codon_translation_list(allelesList):
    # useful for translating codons
    from Bio.Data import CodonTable
    table = CodonTable.ambiguous_dna_by_id[1]
    from Bio.Seq import _translate_str

    # pre-calculating relationships between representations of codons as int and as triplets.
    codonIndices = []
    codonIndices2 = np.zeros((4, 4, 4), dtype=int)
    codonAlleles = {}
    codI = 0
    codonAllelesList = []
    translationList = []
    for i1 in range(4):
        for i2 in range(4):
            for i3 in range(4):
                codonIndices.append((i1, i2, i3))
                codonIndices2[i1, i2, i3] = codI
                codonAlleles[allelesList[i1] + allelesList[i2] + allelesList[i3]] = codI
                codonAllelesList.append(allelesList[i1] + allelesList[i2] + allelesList[i3])
                translationList.append(_translate_str(allelesList[i1] + allelesList[i2] + allelesList[i3], table))
                codI += 1

    return translationList, codonIndices, codonIndices2, codonAlleles, codonAllelesList


def codon_lookup_table(translationList, codonIndices, codonIndices2):
    # Creating quick look-up table that tells you if a mutation is synonymous or not
    isNonsynom = np.zeros((64, 3, 3), dtype=bool)
    isIntoStop = np.zeros((64, 3, 3), dtype=bool)
    for i1 in range(64):
        indeces = codonIndices[i1]
        aa = translationList[i1]
        indeces2 = list(indeces)
        for i2 in range(3):
            for i3 in range(3):
                indeces2[i2] = ((indeces[i2] + i3 + 1) % 4)
                cod2 = codonIndices2[indeces2[0], indeces2[1], indeces2[2]]
                if aa != translationList[cod2]:
                    isNonsynom[i1, i2, i3] = True
                if translationList[cod2] == "*" and aa != "*":
                    isIntoStop[i1, i2, i3] = True
                # print(codonAllelesList[i1]+" "+codonAllelesList[cod2]+" into stop.")
            indeces2[i2] = indeces[i2]

    return isNonsynom, isIntoStop