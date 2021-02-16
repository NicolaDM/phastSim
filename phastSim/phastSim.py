import numpy as np
import argparse
from importlib_resources import files
# from ete3 import Tree



# CONSTANTS
class Constants:
    def __init__(self):
        #self.alleles = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3, "u": 3, "U": 3}
        self.alleles = {"A": 0, "C": 1, "G": 2, "T": 3}
        self.allelesList = ["A", "C", "G", "T"]
        self.nAlleles = 4
        self.stopCodons = ["TAA", "TGA", "TAG"]
        #self.range4=range(4)


def setup_argument_parser():
    parser = argparse.ArgumentParser(
        description='Efficiently simulate sequence evolution along phylogenies with short branches.')
    parser.add_argument('--outpath', default="./", help='Path where to place the results of the simulations.')
    parser.add_argument('--reference', default=None,
                        help='File containing the reference genome to be used as root genome.')
    parser.add_argument("--rootGenomeLength",
                        help="A root genome of this size is randomly created. default=0 means that the root genome "
                             "is instead taken as in the reference file.",
                        type=int, default=0)
    parser.add_argument("--rootGenomeFrequencies",
                        help="Frequencies of different states (non-stop codons or nucleotides depending "
                             "on the simulations). If not provided, will use those for SARS-CoV-2 reference genome.",
                        type=float, nargs='+', default=[1.0])
    parser.add_argument('--treeFile', default=None,
                        help='Name of file containing the tree used to simulate sequences (newick format).')
    parser.add_argument('--scale', default=1.0, type=float,
                        help='Scale the simulation tree by this amount (default 1.0). Branch lengths are assumed'
                             ' in terms of expected substitutions per site '
                             '(more or less, as frequencies changes through time,'
                             ' total mutation rate might also change).')
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
                             "hypermutabilities can be specified. By default no recurrent mutations are simulated. "
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
    parser.add_argument("--createInfo",
                        help="Create a .info file with information regarding the rates at each genome position "
                             "(default name sars-cov-2_simulation_output.info).",
                        action="store_true")
    return parser




class phastSimRun:
    def __init__(self, args):
        self.args = args
        self.const = Constants()
        self.hierarchy = not self.args.noHierarchy
        self.nCodons = None
        self.range4=range(4)

        self.check_default_args()


    def check_default_args(self):
        # if the reference or the tree arguments are not set, then use the default ones shipped with the package
        if not self.args.reference:
            reference = files('phastSim.example').joinpath('MN908947.3.fasta')
            self.args.reference = str(reference)
        if not self.args.treeFile:
            tree = files('phastSim.example').joinpath('exampleTree.tree')
            self.args.treeFile = str(tree)


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
        reference_file = f'{self.args.reference}'
        file = open(reference_file)
        line = file.readline()
        ref = ""
        while line != "":
            line = file.readline()
            ref += line.replace("\n", "")
        file.close()
        print(f"\n Finished reading reference genome at {reference_file} with {str(len(ref))} bases.")
        refList = list(ref)
        #refList = np.array(refList)
        #refList=np.char.split(ref, sep = '')
        return ref, refList


	#Create root genome randomly
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

		#Add stop codon frequencies (set to 0) to list of frequencies
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
                
        range4=self.range4
        allelesList = self.const.allelesList
        codonAllelesList = []
        for i1 in range4:
            for i2 in range4:
                for i3 in range4:
                    codonAllelesList.append(allelesList[i1] + allelesList[i2] + allelesList[i3])
        codList = np.random.choice(codonAllelesList, size=int(rootGenomeLength / 3), replace=True,
                                p=rootGenomeFrequencies)
	
        ref = "".join(codList)
        refList = list(ref)
        #refList = np.array(refList)
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
        refList = list(ref)
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
        categoryRates=np.array(categoryRates)

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
            categories = np.random.choice(nCat, size=self.ref_len, p=categoryProbs)
            self.categories = categories
            gammaRates=categoryRates[categories]

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
        
        omegaCategoryRates=np.array(omegaCategoryRates)

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

            omegaCategories = np.random.choice(nCatOmega, size=nCodons, p=omegaCategoryProbs)
            omegas=omegaCategoryRates[omegaCategories]

        return omegas



class GenomeTree_hierarchical:
    def __init__(self, nCodons, codon, ref, gammaRates, omegas, mutMatrix, hyperCategories, hyperMutRates, file, verbose):
        self.codon = codon
        self.ref = ref
        self.gammaRates = gammaRates
        self.omegas = omegas
        self.mutMatrix = mutMatrix

        self.hyperCategories = hyperCategories
        self.hyperMutRates = hyperMutRates
        self.file = file
        self.verbose = verbose

        const = Constants()
        self.stopCodons = const.stopCodons
        self.alleles = const.alleles
        self.allelesList = const.allelesList
        self.nAlleles = const.nAlleles

        self.alRange = range(const.nAlleles)
        self.range9 = range(9)

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
            
            #create array of all initial rates
            self.allRates=np.zeros((nCodons,10))
            self.range3=range(3)

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
                node.allele = self.codonAlleles[self.ref[pos * 3:(pos + 1) * 3]]
                if node.allele in self.stopCodons:
                    self.omegas[pos] = 0.0
                    if (pos < self.nTerminalNodes - 1):
                    	print(f"Warning: stop codon in the middle of the reference {pos * 3 + 1}"
                          f" {self.ref[pos * 3:(pos + 1) * 3]}."
                          f"Continuing simulations but setting omega=0 for this position.")
                #if (pos < self.nTerminalNodes - 1) and (self.ref[pos * 3:(pos + 1) * 3] in self.stopCodons):
                #    print(f"Warning: stop codon in the middle of the reference {pos * 3 + 1}"
                #          f" {self.ref[pos * 3:(pos + 1) * 3]}."
                #          f"Continuing simulations but setting omega=0 for this position.")
                #    self.omegas[pos] = 0.0
                if self.file!=None:
                	self.file.write(str(pos * 3 + 1) + "-" + str(pos * 3 + 3) + "\t" + str(self.omegas[pos]) + "\t")
                node.rates = {}
                node.hyper = {}
                #Nicola: create initial mega-array of zeros(10) and pass just index to it?
                node.rates[node.allele]= self.allRates[pos]
                #node.rates[node.allele] = np.zeros(10)
                indeces = self.codonIndices[node.allele]
                node.rate = 0.0
                range3=self.range3
                for i2 in range3:
                    pos2 = pos * 3 + i2
                    if self.file!=None:
                    	self.file.write(str(self.gammaRates[pos2]) + "\t" + str(self.hyperCategories[pos2]) + "\t")
                    nuc1 = indeces[i2]
                    for i3 in range3:
                        nuc2 = (nuc1 + i3 + 1) % 4
                        if self.isNonsynom[node.allele, i2, i3]:
                            if self.isIntoStop[node.allele, i2, i3]:
                                node.rates[node.allele][i2 * 3 + i3] = 0.0
                            else:
                                node.rates[node.allele][i2 * 3 + i3] = self.omegas[pos] * self.mutMatrix[nuc1][nuc2] * self.gammaRates[pos2]
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
                        if self.file!=None:
                        	self.file.write(self.allelesList[i] + "\t" + self.allelesList[(i + j + 1) % 4] + "\t")
                    else:
                        if self.file!=None:
                        	self.file.write(".\t" + ".\t")
                    for i3 in range3:
                        node.rate += node.rates[node.allele][i2 * 3 + i3]
                    node.rates[node.allele][9] = node.rate
                if self.file!=None:
                	self.file.write("\n")

            else:
                node.allele = self.alleles[self.ref[pos]]
                # evolutionary categories
                nodeHyper = self.hyperCategories[pos]
                if self.file!=None:
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
                    if self.file!=None:
                    	self.file.write(self.allelesList[i] + "\t" + self.allelesList[j] + "\n")
                else:
                    if self.file!=None:
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
        self.norm = norm

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



    # sample a new allele to mutate to, given the mutation rates at a current node
    # NICOLA: THIS MIGHT BE MADE FASTER, ESPECIALLY WITH LARGE STATE SPACES (E.G. CODON MODELS).
    # INSTEAD OF ITERATING OVER ALL STATES, IN FACT, BELOW ONE COULD USE A DIVIDE AND CONQUER APPROACH
    # OVER THE STATE SPACE, FOR EXAMPLE DEFINING ANOTHER CONSTANT TREE STRUCTURE OVER ALLELE SPACE.
    def sampleMutationCodon(self, rates, rand):
        for j in self.range9:
            if rand < rates[j]:
                return j
            else:
                rand -= rates[j]
        print(
            "You should not have got here - there was a bug somewhere or some unlucky sampling in the machine error area")
        exit()

    def sampleMutation(self, allele, rates, rand):
        for j in self.alRange:
            if j != allele:
                if rand < rates[j]:
                    return j
                else:
                    rand -= rates[j]
        print(
            "You should not have got here - there was a bug somewhere or some unlucky sampling in the machine error area")
        exit()



    def findPos(self, rand, parentGenomeNode, level):
        # find position to mutate along the genome, and update temporary genome tree structure as you go
        if parentGenomeNode.isTerminal:
            # reached a terminal node, now sample the mutation event at the position and update all rates
            node = parentGenomeNode.refNode
            a = parentGenomeNode.allele
            if self.codon:
                j = self.sampleMutationCodon(rates=node.rates[a], rand=rand)
                indeces = self.codonIndices[a]
                i2 = int(j / 3)
                i3 = j % 3
                newIndeces = list(indeces)
                newIndeces[i2] = (newIndeces[i2] + i3 + 1) % 4
                parentGenomeNode.allele = self.codonIndices2[newIndeces[0], newIndeces[1], newIndeces[2]]
                mutEvent = [node.genomePos[0] * 3 + i2, indeces[i2], newIndeces[i2]]
                if self.verbose:
                    print(f"Mutation from {a} {self.codonAllelesList[a]} "
                          f"to {parentGenomeNode.allele} {self.codonAllelesList[parentGenomeNode.allele]}"
                          f" , position {mutEvent[0]} category rate {self.gammaRates[mutEvent[0]]}"
                          f" hyperCat {self.hyperCategories[mutEvent[0]]}"
                          f" omega {self.omegas[node.genomePos[0]]}"
                          f" old rate {node.rates[a][9]} old rates:")
                    print(node.rates[a])
                if not (parentGenomeNode.allele in node.rates):
                    node.rates[parentGenomeNode.allele] = np.zeros(10)
                    indeces = self.codonIndices[parentGenomeNode.allele]
                    parentGenomeNode.rate = 0.0
                    range3=self.range3
                    for i2 in range3:
                        pos2 = node.genomePos[0] * 3 + i2
                        nuc1 = indeces[i2]
                        for i3 in range3:
                            nuc2 = (nuc1 + i3 + 1) % 4
                            if self.isNonsynom[parentGenomeNode.allele, i2, i3]:
                                if self.isIntoStop[parentGenomeNode.allele, i2, i3]:
                                    node.rates[parentGenomeNode.allele][i2 * 3 + i3] = 0.0
                                else:
                                    node.rates[parentGenomeNode.allele][i2 * 3 + i3] = self.omegas[node.genomePos[0]] * \
                                                                                       self.mutMatrix[nuc1][nuc2] * \
                                                                                       self.gammaRates[pos2] / self.norm
                            else:
                                node.rates[parentGenomeNode.allele][i2 * 3 + i3] = self.mutMatrix[nuc1][nuc2] * \
                                                                                   self.gammaRates[pos2] / self.norm
                        if self.hyperCategories[pos2] > 0:
                            if node.hyper[i2][0] == nuc1:
                                node.rates[parentGenomeNode.allele][i2 * 3 + node.hyper[i2][1]] *= self.hyperMutRates[
                                    self.hyperCategories[pos2] - 1]
                        for i3 in range3:
                            node.rates[parentGenomeNode.allele][9] += node.rates[parentGenomeNode.allele][i2 * 3 + i3]
                parentGenomeNode.rate = node.rates[parentGenomeNode.allele][9]
                if self.verbose:
                    print(f" new rate {parentGenomeNode.rate} all rates:")
                    print(node.rates[parentGenomeNode.allele])
            else:
                j = self.sampleMutation(a, node.rates[a], rand)
                mutEvent = [node.genomePos[0], a, j]
                if self.verbose:
                    print(f"Mutation from {a} to {j},"
                          f" position {node.genomePos[0]}"
                          f" category rate {self.gammaRates[node.genomePos[0]]}"
                          f" hyperCat {self.hyperCategories[node.genomePos[0]]}"
                          f" old rate {parentGenomeNode.rate} old rates:")
                    print(node.rates)
                parentGenomeNode.rate = -node.rates[j, j]
                if self.verbose:
                    print(" new rate " + str(parentGenomeNode.rate) + " all rates:")
                    print(node.rates)
                parentGenomeNode.allele = j
            return mutEvent

        else:
            # still at an internal genome node.
            # choose which of the two children genome nodes to move into
            if rand >= parentGenomeNode.belowNodes[0].rate:
                rand -= parentGenomeNode.belowNodes[0].rate
                parentGenomeNode.rate = parentGenomeNode.belowNodes[0].rate
                child = parentGenomeNode.belowNodes[1]
                childI = 1
            else:
                child = parentGenomeNode.belowNodes[0]
                parentGenomeNode.rate = parentGenomeNode.belowNodes[1].rate
                childI = 0

            # if the child we are moving into is not on the same level, but is above,
            # then create a new child at the same level.
            # this is because the rate of the child will be inevitably changed by the mutation event,
            # and we don't want to change the mutation rates for the parent phylogenetic node.
            if child.level < level:
                newChild = genomeNode(level=level)  # upNode=parentGenomeNode
                parentGenomeNode.belowNodes[childI] = newChild
                newChild.isTerminal = child.isTerminal
                if child.isTerminal:
                    newChild.refNode = child.refNode
                    newChild.allele = child.allele
                else:
                    newChild.belowNodes = list(child.belowNodes)

                mutEvent = self.findPos(rand, newChild, level)
                parentGenomeNode.rate += newChild.rate
            else:
                # in this case the child is already on the same level, so no need to create another one, just update its mutation rate.
                mutEvent = self.findPos(rand, child, level)
                parentGenomeNode.rate += child.rate
            return mutEvent


    def mutateBranchETEhierarchy(self, childNode, parentGenomeNode, level, createNewick):
        # Function to simulate evolution on one branch,using ETE tree structure
        # and using the genome-wide hierarchy structure.
        # given details of the parent node, it generates details of the child node, and updates the hierarchy accordingly.
        # To simulate evolution on the whole tree, it needs to be called on the root.

        # branch length above the current node
        bLen = childNode.dist
        currTime = 0.0
        # if newick output is requested, prepare format
        if createNewick:
            childNode.mutAnnotation = []
        # Initialize child rate and allele numbers with parent ones
        rate = parentGenomeNode.rate
        childNode.mutations = []

        # Sample new mutation event with Gillespie algorithm
        currTime += np.random.exponential(scale=1.0 / rate)
        if self.verbose:
            print(f"\n Node {childNode.name} BLen: {bLen} first sampled time: {currTime}; mutation rate: {rate}")
        # for the first mutation event at this node, create a new root genome node of the appropriate level.
        # otherwise, use the one you already have.
        if currTime < bLen:
            newGenomeNode = genomeNode(level=level)
            newGenomeNode.belowNodes = list(parentGenomeNode.belowNodes)
        else:
            newGenomeNode = parentGenomeNode
        while currTime < bLen:
            # Now, sample which type of mutation event it is (from which nucleotide to which nucleotide)
            rand = np.random.random() * rate
            if self.verbose:
                print("Selecting new mutation event. Rate " + str(rate) + " random value " + str(rand))
            mutEvent = self.findPos(rand, newGenomeNode, level)
            childNode.mutations.append(mutEvent)
            if createNewick:
                childNode.mutAnnotation.append(self.allelesList[mutEvent[1]] +
                                               str(mutEvent[0] + 1) +
                                               self.allelesList[mutEvent[2]])
            rate = newGenomeNode.rate
            if self.verbose:
                print(f"New total rate {rate}")
            currTime += np.random.exponential(scale=1.0 / rate)
            if self.verbose:
                print(f"new time {currTime}, rate {rate} mutation events:")
                print(childNode.mutations)

        if self.verbose:
            print("mutations at the end:")
            print(childNode.mutations)
        # now mutate children of the current node, calling this function recursively on the node children.
        for c in childNode.children:
            self.mutateBranchETEhierarchy(c, newGenomeNode, level + 1, createNewick)



    def writeGenomeShort(self, node, file, mutDict):
        # function to write a succint output iteratively
        # update dictionary
        for m in node.mutations:
            nuc = self.allelesList[m[2]]
            if nuc != self.ref[m[0]]:
                mutDict[m[0] + 1] = nuc
            else:
                del mutDict[m[0] + 1]
        # print leaf entry to file
        if node.is_leaf():
            file.write(">" + node.name + "\n")
            mutList = list(mutDict.keys())
            mutList.sort()
            for m in mutList:
                file.write(str(m) + " " + mutDict[m] + "\n")
        # pass dictionary to children
        else:
            for c in node.children:
                self.writeGenomeShort(c, file, mutDict)
        # de-update the dictionary so it can be used by siblings etc.
        for n in range(len(node.mutations)):
            m = node.mutations[len(node.mutations) - (n + 1)]
            nuc = self.allelesList[m[1]]
            if nuc != self.ref[m[0]]:
                mutDict[m[0] + 1] = nuc
            else:
                del mutDict[m[0] + 1]


    def write_genome_short(self, tree, output_path, output_file):
        # open a file a create a container
        genomefile = open(output_path + output_file + ".txt", "w")
        mutDict = {}
        # call the recursive function
        self.writeGenomeShort(node=tree, file=genomefile, mutDict=mutDict)
        genomefile.close()


    def writeGenome(self, node, file, nRefList):
        # function to write a complete sequence output iteratively
        # update list
        for m in node.mutations:
            nRefList[m[0]] = self.allelesList[m[2]]
        # print leaf entry to file
        if node.is_leaf():
            # first write the header then the sequence directly from the array
            #file.write(">" + node.name + "\n")# + (''.join(nRefList)) + "\n")
            #np.savetxt(file, nRefList,fmt='%s')
            #nRefList.tofile(file)
            file.write(">" + node.name + "\n"+''.join(nRefList) + "\n")
            #file.write(">" + node.name + "\n" + np.char.join(nRefList,"") + "\n")
        # pass to children
        else:
            for c in node.children:
                self.writeGenome(c, file, nRefList)
        # de-update the list so it can be used by siblings etc.
        for n in range(len(node.mutations)):
            m = node.mutations[len(node.mutations) - (n + 1)]
            nRefList[m[0]] = self.allelesList[m[1]]


    def write_genome(self, tree, output_path, output_file, refList):
        file = open(output_path + output_file + ".fasta", "w")
        # convert reference list to an array
        #refList = np.array(refList)
        self.writeGenome(tree, file, refList)
        file.close()



    def writeGenomePhylip(self, node, file, nRefList):
        # function to write a complete sequence output iteratively
        # update list
        for m in node.mutations:
            nRefList[m[0]] = self.allelesList[m[2]]
        # print leaf entry to file
        if node.is_leaf():
            file.write(node.name + "\t" + (''.join(nRefList)) + "\n")
        else:
            for c in node.children:
                self.writeGenomePhylip(c, file, nRefList)
        # de-update the list so it can be used by siblings etc.
        for n in range(len(node.mutations)):
            m = node.mutations[len(node.mutations) - (n + 1)]
            nRefList[m[0]] = self.allelesList[m[1]]


    def write_genome_phylip(self, tree, output_path, output_file, refList):
        # open a file for the phylip output
        file = open(output_path + output_file + ".phy", "w")
        file.write("\t" + str(len(tree)) + "\t" + str(len(self.ref)) + "\n")
        # run the recursive function to write the phylip formatted file
        self.writeGenomePhylip(node=tree, file=file, nRefList=refList)
        file.close()


class GenomeTree_vanilla:
    def __init__(self, nCat, ref, mutMatrix, categories, categoryRates, hyperMutRates, hyperCategories, file, verbose):
        self.nCat = nCat
        self.ref = ref
        self.mutMatrix = mutMatrix
        self.categories = categories
        self.categoryRates = categoryRates
        self.hyperMutRates = hyperMutRates
        self.hyperCategories = hyperCategories

        self.file = file
        self.verbose = verbose

        const = Constants()
        self.alleles = const.alleles
        self.allelesList = const.allelesList
        self.range4=range(4)

        self.muts = []
        for c in range(nCat):
            self.muts.append([[], [], [], []])


    def prepare_genome(self):
        # List of dictionaries that let you know at wich position (starting from 0) is the 1st A of the genome,
        # the second A of the genome, etc (well, actually the 0th, the 1st, etc..).
        # here I could use arrays instead of dictionaries, it might be more efficient.
        positions = []
        for c in range(self.nCat):
            positions.append([[], [], [], []])
        # Total mutation rates cumulatively across the genome
        totMutMatrix = np.zeros((self.nCat, 4, 4), dtype=float)
        # hyper rates of currently affected sites
        extras = []
        # total mutation rate
        totMut = 0.0
        # tot number of A's, C's, G's and T's in the reference genome
        # Number of alleles in each category across the genome.
        totAlleles = np.zeros((self.nCat, 4))
        for pos in range(len(self.ref)):
            a = self.alleles[self.ref[pos]]
            cat = self.categories[pos]
            hyp = self.hyperCategories[pos]
            if self.file!=None:
            	self.file.write(str(pos + 1) + "\t" + str(cat) + "\t" + str(hyp) + "\t")
            # now sample which alleles are affected by hypermutation and store info in positions and extra vectors
            if hyp > 0:
                i = np.random.choice(4)
                j = np.random.choice(3)
                j = (i + j + 1) % 4
                positions[cat][a].append([pos, hyp, i, j])
                if i == a:
                    extras.append([self.mutMatrix[a][j] * self.categoryRates[cat] * (self.hyperMutRates[hyp - 1] - 1.0),
                                   pos, len(positions[cat][a]) - 1, cat, hyp, i, j])
                    totMut += self.mutMatrix[a][j] * self.categoryRates[cat] * (self.hyperMutRates[hyp - 1] - 1.0)
                if self.file!=None:
                	self.file.write(self.allelesList[i] + "\t" + self.allelesList[j] + "\n")
            else:
                positions[cat][a].append([pos, 0])
                if self.file!=None:
                	self.file.write(".\t" + ".\n")
            #for j in self.range4:
            #    if j != a:
            #        totMutMatrix[cat][a][j] += self.mutMatrix[a][j] * self.categoryRates[cat]
            #        totMut += self.mutMatrix[a][j] * self.categoryRates[cat]
            totAlleles[self.categories[pos]][a] += 1
        for c in range(self.nCat):
        	for a in self.range4:
        		for j in self.range4:
        			if j != a:
        				totMut += self.mutMatrix[a][j] * self.categoryRates[c]*totAlleles[c][a]
        				totMutMatrix[cat][a][j]+= self.mutMatrix[a][j] * self.categoryRates[c]*totAlleles[c][a]
        if self.file!=None:
        	self.file.close()

        print("\n Number of each nucleotide in the genome:")
        print(totAlleles)

        # set the total mutation rate and cumulative rates as class attributes
        self.totMut = totMut
        self.totMutMatrix = totMutMatrix
        self.extras = extras
        self.totAlleles = totAlleles
        self.positions = positions


    def normalize_rates(self, scale):
        norm = self.totMut / len(self.ref)
        print(f"\n Total cumulative mutation rate per site before normalization: {norm}")

        # We rescale by the input normalization factor, this is the same as rescaling
        # all the branch lengths by this rescaling factor
        norm = norm / scale
        print(f"\n After rescaling: {norm}")

        # Now normalize mutation rates so that, at the start,
        # the expected number of substitutions per unit branch length is 1?
        for i in range(4):
            for j in range(4):
                self.mutMatrix[i][j] = self.mutMatrix[i][j] / norm
                for c in range(self.nCat):
                    if j != i:
                        self.totMutMatrix[c][i][j] = self.totMutMatrix[c][i][j] / norm
        for e in range(len(self.extras)):
            self.extras[e][0] = self.extras[e][0] / norm
        self.totMut = self.totMut / norm

        if self.verbose:
            print("Base-wise total mutation rates:")
            print(self.totMutMatrix)
            print("Normalized mutation rates:")
            print(self.mutMatrix)



    def mutateBranchETE(self, childNode, parentMuts, parentTotAlleles, parentRate, extrasParent, createNewick):
        # Function to simulate evolution on one branch,using ETE tree structure;
        # given details of the parent node, it generates details of the child node.
        # To simulate whole tree, it needs to be called on root.

        bLen = childNode.dist
        currTime = 0.0
        # if newick output is requested, prepare format
        if createNewick:
            childNode.mutAnnotation = []
        # Initialize child rate and allele numbers with parent ones
        rate = parentRate
        childTotAlleles = []
        childMutations = []
        for c in range(self.nCat):
            childTotAlleles.append(list(parentTotAlleles[c]))
            childMutations.append([[], [], [], []])
            # Initialize child mutation list with parent one
            for i in range(4):
                for k in range(len(parentMuts[c][i])):
                    childMutations[c][i].append(list(parentMuts[c][i][k]))
        extrasChild = []
        for e in range(len(extrasParent)):
            extrasChild.append(list(extrasParent[e]))

        # Sample new mutation event with Gillespie algorithm
        currTime += np.random.exponential(scale=1.0 / rate)
        if self.verbose:
            print(f"\n Node {childNode.name} BLen: {bLen}. Mutations at start:")
            print(parentMuts)
            print("num alleles at start:")
            print(parentTotAlleles)
            print(f"tot rate at start: {parentRate}")
            print(f"First sampled time: {currTime}; mutation rate: {rate}")

        while currTime < bLen:
            # Now, sample which type of mutation event it is (from which nucleotide to which nucleotide)
            rand = np.random.random() * rate
            if self.verbose:
                print(f"Selecting new mutation event. Rate {rate} random value {rand}")
            tot = 0.0
            found = False
            hyperExtra = False
            for c in range(self.nCat):
                for i in range(4):
                    for j in range(4):
                        if j != i:
                            tot += childTotAlleles[c][i] * self.mutMatrix[i][j] * self.categoryRates[c]
                            # print("i "+str(i)+" j "+str(j)+" tot "+str(tot))
                            if rand < (tot):
                                found = True
                                break
                    if found:
                        break
                if found:
                    break
            if not found:
                # now use extras vector for hypermutable sites
                for e in range(len(extrasChild)):
                    tot += extrasChild[e][0]
                    if rand < (tot):
                        found = True
                        hyperExtra = True
                        extra = extrasChild.pop(e)
                        break
            if not found:
                print("Error in selecting mutation type")
                exit()

            if hyperExtra:
                # element already removed from extraChild list, now use info in extra to add,
                # remove or modify entry from mutations list, amend total rates, etc.
                a = self.alleles[self.ref[extra[1]]]
                extraRate = extra[0]
                rate -= extraRate
                pos = extra[1]
                pos2 = extra[2]
                c = extra[3]
                hyp = extra[4]
                i = extra[5]
                j = extra[6]
                if self.verbose:
                    print(f"Hypermutation genome position {pos + 1} category {c} hypCat {hyp}"
                          f" allele {i} allele position {pos2} to {j}")
                if createNewick:
                    childNode.mutAnnotation.append(self.allelesList[i] + str(pos + 1) + self.allelesList[j])

                m = 0
                included = False
                while m < len(childMutations[c][a]):
                    if childMutations[c][a][m][0] == pos2:
                        if self.verbose:
                            print("Position had already mutated")
                        if childMutations[c][a][m][1] != i:
                            print("Error, mutation should have been recorded as into hypermutable allele.")
                            print(childMutations[c][a][m])
                            print(extra)
                            print(self.ref[extra[1]])
                            print(c)
                            print(a)
                            exit()
                        childMutations[c][a][m][1] = j
                        included = True
                        break
                    elif childMutations[c][a][m][0] < pos2:
                        m += 1
                    else:
                        if self.ref[pos] != self.allelesList[i]:
                            print("Error, hypermutable allele should have been the reference.")
                            exit()
                        childMutations[c][a].insert(m, [pos2, j])
                        if self.verbose:
                            print("New mutation has been introduced")
                            print(childMutations[c][a][m])
                        included = True
                        break
                if not included:
                    childMutations[c][a].append([pos2, j])


            else:
                # Now, sample the specific position of the genome (among those with the mutated alleles) that mutates.
                mutatedBasePos = np.random.randint(childTotAlleles[c][i])
                if self.verbose:
                    print(f"mutation position {mutatedBasePos} category {c} allele {i} to {j}")

                # in this case, we are mutating a position that was not mutated before.
                # We need to add one entry to the mutation list without removing any old one.
                if mutatedBasePos < self.totAlleles[c][i] - len(childMutations[c][i]):
                    # print("New mutation")
                    newMutPos = mutatedBasePos
                    m = 0
                    while m < len(childMutations[c][i]):
                        if childMutations[c][i][m][0] <= newMutPos:
                            newMutPos += 1
                            m += 1
                        else:
                            break
                    childMutations[c][i].insert(m, [newMutPos, j])
                    if createNewick:
                        childNode.mutAnnotation.append(
                            self.allelesList[i] + str(self.positions[c][i][newMutPos][0] + 1) + self.allelesList[j])

                    infoExtra = self.positions[c][i][newMutPos]
                    if infoExtra[1] > 0:
                        if infoExtra[2] == i:
                            # remove hypermutation from extras
                            if self.verbose:
                                print(f"removing from extras {infoExtra[0]} {i} {j}")
                            removed = False
                            for e in range(len(extrasChild)):
                                if infoExtra[0] == extrasChild[e][1]:
                                    extra = extrasChild.pop(e)
                                    rate -= extra[0]
                                    removed = True
                                    break
                            if not removed:
                                print("Error, hyper mutation not removed")
                                exit()
                            elif self.verbose:
                                print("Hyper mutation removed by normal mutation")
                                print(extra)
                        elif infoExtra[2] == j:
                            extraRate = self.mutMatrix[i][j] * self.categoryRates[c] * \
                                        (self.hyperMutRates[infoExtra[1] - 1] - 1.0)
                            extrasChild.append(
                                [extraRate, infoExtra[0], newMutPos, c, infoExtra[1], infoExtra[2], infoExtra[3]])
                            rate += extraRate

                # in this case, we are mutating a position that was already mutated.
                # this means that one item from the mutation list needs to be removed or modified,
                # and no item needs to be added.
                else:
                    newMutatedBasePos = mutatedBasePos - (self.totAlleles[c][i] - len(childMutations[c][i]))
                    if self.verbose:
                        print("Modifying pre-existing mutation")
                        print(newMutatedBasePos)
                    added = False
                    for i2 in range(4):
                        if i2 != i:
                            for m in range(len(childMutations[c][i2])):
                                if childMutations[c][i2][m][1] == i:
                                    if newMutatedBasePos == 0:
                                        if createNewick:
                                            childNode.mutAnnotation.append(self.allelesList[i] + str(
                                                self.positions[c][i2][childMutations[c][i2][m][0]][0] + 1) +
                                                                           self.allelesList[j])
                                        infoExtra = self.positions[c][i2][childMutations[c][i2][m][0]]
                                        if infoExtra[1] > 0:
                                            if infoExtra[2] == i:
                                                # remove hypermutation from extras
                                                for e in range(len(extrasChild)):
                                                    if infoExtra[0] == extrasChild[e][1]:
                                                        extra = extrasChild.pop(e)
                                                        rate -= extra[0]
                                                        break
                                            elif infoExtra[2] == j:
                                                extraRate = self.mutMatrix[i][j] * self.categoryRates[c] * (
                                                        self.hyperMutRates[infoExtra[1] - 1] - 1.0)
                                                extrasChild.append(
                                                    [extraRate, infoExtra[0], childMutations[c][i2][m][0], c,
                                                     infoExtra[1],
                                                     infoExtra[2], infoExtra[3]])
                                                rate += extraRate
                                        if j == i2:
                                            if self.verbose:
                                                print("Deleted mutation \n\n\n")
                                            del childMutations[c][i2][m]
                                        else:
                                            childMutations[c][i2][m][1] = j
                                        added = True
                                        break
                                    newMutatedBasePos -= 1
                        if added:
                            break

            childTotAlleles[c][i] -= 1
            childTotAlleles[c][j] += 1
            rate += self.mutMatrix[i][i] * self.categoryRates[c]
            rate -= self.mutMatrix[j][j] * self.categoryRates[c]

            currTime += np.random.exponential(scale=1.0 / rate)
            if self.verbose:
                print(f"new time {currTime}, rate {rate} mutation events:")
                print(childMutations)

        if self.verbose:
            print("mutations at the end:")
            print(childMutations)
        childNode.mutations = childMutations
        # now mutate children of the current node, calling this function recursively on the node children.
        for c in childNode.children:
            self.mutateBranchETE(c, childMutations, childTotAlleles, rate, extrasChild, createNewick)



    def writeGenomeShort(self, node, file):
        # function to write a succint output iteratively
        if node.is_leaf():
            file.write(">" + node.name + "\n")
            mutDict = {}
            for c in range(self.nCat):
                for i in range(4):
                    for m in node.mutations[c][i]:
                        mutDict[self.positions[c][i][m[0]][0] + 1] = self.allelesList[m[1]]
            mutList = list(mutDict.keys())
            mutList.sort()
            for m in mutList:
                file.write(str(m) + " " + mutDict[m] + "\n")
        for c in node.children:
            self.writeGenomeShort(c, file)



    def write_genome_short(self, tree, output_path, output_file):
        # open a file
        genomefile = open(output_path + output_file + ".txt", "w")
        # call the recursive function
        self.writeGenomeShort(node=tree, file=genomefile)
        genomefile.close()


    def genomeSeq(self, mutations, refList):
        # generate the genome sequence of a sample using node mutations and the reference
        # useful, for example, for generating a fasta file.
        # used in writing out the full fasta genome
        for c in range(self.nCat):
            for i in range(4):
                for m in mutations[c][i]:
                    refList[self.positions[c][i][m[0]][0]] = self.allelesList[m[1]]
        newGenome = ''.join(refList)
        for c in range(self.nCat):
            for i in range(4):
                for m in mutations[c][i]:
                    pos = self.positions[c][i][m[0]][0]
                    refList[pos] = self.ref[pos]
        return newGenome


    def writeGenome(self, node, file, refList):
        # function to write a fasta output iteratively
        if node.is_leaf():
            seq = self.genomeSeq(mutations=node.mutations, refList=refList)
            file.write(">" + node.name + "\n" + seq + "\n")
        for c in node.children:
            self.writeGenome(c, file, refList)


    def write_genome(self, tree, output_path, output_file, refList):
        # class specific wrapper to write the full fasta genome
        # open a file
        file = open(output_path + output_file + ".fasta", "w")
        # run the recursive function
        self.writeGenome(tree, file, refList)
        file.close()



    def writeGenomePhylip(self, node, file, nRefList):
        # function to write a phylip output iteratively
        if node.is_leaf():
            seq = self.genomeSeq(mutations=node.mutations, refList=nRefList)
            file.write(node.name + "\t" + seq + "\n")
        for c in node.children:
            self.writeGenomePhylip(c, file, nRefList)


    def write_genome_phylip(self, tree, output_path, output_file, refList):
        # open a file for the phylip output
        file = open(output_path + output_file + ".phy", "w")
        file.write("\t" + str(len(tree)) + "\t" + str(len(self.ref)) + "\n")
        # run the recursive function to write the phylip formatted file
        self.writeGenomePhylip(node=tree, file=file, nRefList=refList)
        file.close()




class genomeNode:
    # node representing an element of the genome hierarchy, summarizing the total rate of the nodes below it,
    # that is, the total rate of part of the genome.

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



def writeGenomeNewick(node):
    # function to write a newick output iteratively
    if node.is_leaf():
        outString = node.name + '['
    else:
        outString = '('
        for c in range(len(node.children)):
            outString += writeGenomeNewick(node.children[c])
            if c < len(node.children) - 1:
                outString += ','
        outString += ')['
    stringToAdd = ''
    for i in range(len(node.mutAnnotation)):
        stringToAdd += node.mutAnnotation[i]
        if i < len(node.mutAnnotation) - 1:
            stringToAdd += ','
    stringToAdd += (']:' + str(node.dist))
    return outString + stringToAdd
