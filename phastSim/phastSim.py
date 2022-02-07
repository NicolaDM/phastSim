import numpy as np
import argparse
from importlib_resources import files
from enum import Enum
from collections import Counter
from copy import deepcopy
from itertools import count, chain, product
import phastSim.parsimony_pb2 as protobuf

# CONSTANTS
class Constants:
    def __init__(self):
        self.alleles = {"A": 0, "C": 1, "G": 2, "T": 3}
        self.allelesList = ["A", "C", "G", "T"]
        self.nAlleles = 4
        self.stopCodons = ["TAA", "TGA", "TAG"]


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
                             "on the simulations). If not provided, will use those for SARS-CoV-2 reference genome."
                             "These frequencies are also used when generating indels.",
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
                        help="Mutation rates, by default using the neutral rates estimated from SARS-CoV-2 so far (UNREST model)."
                             "Possible inputs are: "
                             "UNREST r_AC r_AG r_AT r_CA etc... (12 values)"
                             "JC69"
                             "HKY85 r_transition r_transversion pi_A pi_C pi_G pi_T"
                             "GTR r_AC r_AG r_AT r_CG r_CT r_GT pi_A pi_C pi_G pi_T",
                        nargs='+',
                        default=["UNREST", 0.039, 0.310, 0.123, 0.140, 0.022, 3.028, 0.747, 0.113, 2.953, 0.056, 0.261, 0.036])
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
    parser.add_argument("--omegaBeta",
                        help="Second parameter of the gamma distribution for omega variation, in case option omegaAlpha is used."
                             " This one specifies beta. The mean of the distribution will be alpha/beta, and the variance alpha/(beta^2).",
                        type=float, default=1.0)
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
    parser.add_argument("--noNormalization",
                        help="Run without normalizing the mutation rates - they are however still scaled by --scale.",
                        action="store_true")
    parser.add_argument("--verbose", help="Turns on verbose mode.", action="store_true")
    parser.add_argument('--outputFile', default="sars-cov-2_simulation_output",
                        help='Output file name containing the simulated genomes in succint format. The file will be '
                             'created within the folder specified with --path.')
    parser.add_argument('--alternativeOutputFormat', help='Produces a succinct txt output file using a modified VCF format', 
                        action="store_true")
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
    parser.add_argument("--createMAT", 
                        help="Create a mat.pb protocol buffer output file (default name sars-cov-2_simulation_output.mat.pb).",
                        action="store_true")
    parser.add_argument("--indels",
                        help="Allow for insertions and deletions with lengths drawn from a set of configurable parametric distributions",
                        action="store_true")
    parser.add_argument("--insertionRate",
                        help="Distribution for the rates of insertions at each site. Possible values are: "
                             "CONSTANT x - all sites in the genome have a rate of x.\n"
                             "GAMMA ins_alpha ins_beta - rates drawn from a gamma(ins_alpha, ins_beta) distribution,\n"
                             "this distribution has mean ins_alpha/ins_beta",
                        type=str, nargs='+', default=None)
    parser.add_argument("--deletionRate",
                        help="Distribution for the rates of insertions at each site. Possible values are:\n "
                             "CONSTANT x - all sites in the genome have a rate of x.\n "
                             "GAMMA del_alpha del_beta - rates drawn from a gamma(del_alpha, del_beta) distribution,"
                             "this distribution has mean del_alpha/del_beta",
                        type=str, nargs='+', default=None)
    parser.add_argument("--insertionLength",
                        help="Genome-wide distributional parameters for insertion lengths at each site."
                             "Possible values are:\n"
                             "GEOMETRIC p - indels distributed with probability mass function (1-p)^k * p.\n"
                             "NEGBINOMIAL n p - distributed as C(k + n - 1, k) * p^k * (1-p)^-n\n"
                             "ZIPF a - distributed as k^-a / Zeta(a)\n"
                             "LAVALETTE M a\n"
                             "DISCRETE p_1 p_2 p_3 p_4 ... etc. (N values, p_k the probability of an indel of length k)",
                        type=str, nargs='+', default=None)
    parser.add_argument("--deletionLength",
                        help="Genome-wide distributional parameters for deletion lengths at each site."
                             "Possible values are:\n"
                             "GEOMETRIC p - indels distributed with probability mass function (1-p)^k * p.\n"
                             "NEGBINOMIAL n p - distributed as C(k + n - 1, k) * p^k * (1-p)^-n\n"
                             "ZIPF a - distributed as k^-a / Zeta(a)\n"
                             "LAVALETTE M a\n"
                             "DISCRETE p_1 p_2 p_3 p_4 ... etc. (N values, p_k the probability of an indel of length k)",
                        type=str, nargs='+', default=None)
    parser.add_argument('--mutationsTSVinput', 
                        help='Name of (optional) file containing mutation events that are enforced on the tree.'
                        'Sites mutated this way become unmutable through the standard simulated mutation process.', 
                        default=None)
    parser.add_argument('--eteFormat',
                        help='Set an ete3 parsing mode with which to read the input tree newick.',
                        type=int,default=0)
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

        rootGenomeFrequencies = self.init_insertion_frequencies(None)
        rootGenomeLength = self.args.rootGenomeLength


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

        rootGenomeFrequencies = self.init_insertion_frequencies(None)
        rootGenomeLength = self.args.rootGenomeLength

        

        allelesList = self.const.allelesList
        refList = np.random.choice(allelesList, size=int(rootGenomeLength), replace=True, p=rootGenomeFrequencies)
        ref = "".join(refList)
        refList = list(ref)
        return ref, refList


    def init_substitution_rates(self):
        # substitution rates
        mutationInputs = self.args.mutationRates
        model = mutationInputs[0]

        if not model.upper() in ["UNREST", "GTR", "JC69", "HKY85"]:
            print(f"This model, {model}, is not yet implemented, choose one of UNREST, GTR, JC69, HKY85, followed by comma separated parameter values.")
            print("Use phastSim --help for more information.")
            exit()

        mutationInputs = [float(x) for x in mutationInputs[1:]]

        print(f"Assuming {model} nucleotide mutation model.")

        mutationRates = []
        if model.upper() == "UNREST":
            
            if len(mutationInputs) != 12:
                print(f"{model.upper()} requires exactly 12 parameters, r_AC, r_AG, r_AT, r_CA, r_CG, r_CT, r_GA, r_GC, r_GT, r_TA, r_TC, r_TG.")
                print(f"This model had {len(mutationInputs)} parameter values.") 
                exit()
            mutationRates = mutationInputs

        elif model.upper() == "GTR":
            if len(mutationInputs) != 10:
                print(f"{model.upper()} requires exactly 10 parameters, r_AC, r_AG, r_AT, r_CG, r_CT, r_GT, pi_A, pi_C, pi_G, pi_T.")
                print(f"This model had {len(mutationInputs)} parameter values.") 
                exit()

            r_AC, r_AG, r_AT, r_CG, r_CT, r_GT, pi_A, pi_C, pi_G, pi_T = tuple(mutationInputs)
            
            mutationRates = [
                                          r_AC * pi_C, r_AG * pi_G, r_AT * pi_T,
                             r_AC * pi_A,              r_CG * pi_G, r_CT * pi_T, 
                             r_AG * pi_A, r_CG * pi_C,              r_GT * pi_T, 
                             r_AT * pi_A, r_CT * pi_C, r_GT * pi_G]

        elif model.upper() == "JC69":

            mutationRates = [0.25] * 12

        elif model.upper() == "HKY85":
            
            if len(mutationInputs) != 6:
                print(f"{model.upper()} requires exactly 6 parameters, r_transition, r_transversion, pi_A, pi_C, pi_G, pi_T.")
                print(f"This model had {len(mutationInputs)} parameter values.") 
                exit()

            # a = alpha = r_transition, b = beta = r_transversion
            a, b, pi_A, pi_C, pi_G, pi_T = tuple(mutationInputs)
            mutationRates = [
                                       b * pi_C, a * pi_G, b * pi_T,
                             b * pi_A,           b * pi_G, a * pi_T, 
                             a * pi_A, b * pi_C,           b * pi_T, 
                             b * pi_A, a * pi_C, b *pi_G]

        else:
            exit()
        
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
            
            gammaRates = iter(np.random.gamma(self.args.alpha, 1.0 / self.args.alpha, size=self.ref_len))
            if self.args.indels:
                gammaRates = chain(gammaRates, (np.random.gamma(self.args.alpha, 1.0 / self.args.alpha) for _ in count()))

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
            
            categories = iter(np.random.choice(nCat, p=categoryProbs, size=self.ref_len))
            if self.args.indels:
                categories = chain(categories, (np.random.choice(nCat, p=categoryProbs) for _ in count()))

            self.categories = categories
            gammaRates=(categoryRates[c] for c in categories)

        invariable = self.args.invariable
        if invariable >= 0.000000001:
            print(f"Proportion of invariable {invariable}")
            categoriesInv = iter(np.random.choice(2, p=[1.0 - invariable, invariable], size=self.ref_len))
            if self.args.indels:
                categoriesInv = (np.random.choice(2, p=[1.0 - invariable, invariable]) for _ in count())
            gammaRates = (0.0 if c else g for (c,g) in zip(categoriesInv, gammaRates))

        mutationsTSVinput=self.args.mutationsTSVinput
        if (self.args.indels or not self.hierarchy or self.args.codon) and mutationsTSVinput:
            print("TSV mutations is currently not supported alongside indels or a codon model, and must use hierarchy mode.")
            print("Please remove --codon, --indels and --noHierarchy flags.")
            exit()

        if mutationsTSVinput!=None:
            mutation_file = f'{mutationsTSVinput}'
            file = open(mutation_file)
            line=file.readline()
            preMutationsBranches={}
            immutablePositions = []
            while line!="" and line!="\n": 
                linelist=line.split()
                if len(linelist)>1:
                    branch=linelist[0]
                    preMutationsBranches[branch]=[]
                    linelist2=linelist[1].split(",")
                    for i in range(len(linelist2)):
                        pos=int(linelist2[i][1:-1])
                        immutablePositions.append(pos - 1)
                        preMutationsBranches[branch].append((pos,linelist2[i][0],linelist2[i][-1]))
                line=file.readline()
            file.close()
            
            gammaRate = (0.0 if c in immutablePositions else g for (c,g) in zip(count(),gammaRates))
            print(preMutationsBranches)
            return gammaRates, preMutationsBranches

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
            print("WARNING: hypermutable sites are supposed to be rare, but total proportion is " + str(sumHyper))
        if sumHyper > 1.0:
            exit()
        newHyperMutProbs = [1.0 - sumHyper] + hyperMutProbs

        # sample hypermutability for each site of the genome.
        hyperCategories = iter(np.random.choice(nProbs + 1, p=newHyperMutProbs, size=self.ref_len))
        if self.args.indels:
            hyperCategories = chain(hyperCategories, (np.random.choice(nProbs + 1, p=newHyperMutProbs) for _ in count()))
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
        omegaBeta = self.args.omegaBeta
        omegaCategoryProbs = self.args.omegaCategoryProbs
        omegaCategoryRates = self.args.omegaCategoryRates
        
        omegaCategoryRates=np.array(omegaCategoryRates)

        nCodons = int(self.ref_len / 3)
        
        self.nCodons = nCodons
        print("Using a codon model")

        if omegaAlpha >= 0.000000001:
            print("Using a continuous gamma distribution with parameter alpha=" + str(omegaAlpha) + 
                " and beta=" +str(omegaBeta)+ " for variation in omega across codons.")

            omegas = iter(np.random.gamma(omegaAlpha, 1.0 / omegaBeta, size=self.nCodons))

            if self.args.indels:
                omegas = chain(omegas, (np.random.gamma(omegaAlpha, 1.0 / omegaBeta) for _ in count()))
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

            omegas = iter(np.random.choice(omegaCategoryRates, p=omegaCategoryProbs, size=self.nCodons))
            if self.args.indels:
                omegas = chain(omegas, (np.random.choice(omegaCategoryRates, p=omegaCategoryProbs) for _ in count()))

        return omegas

    def init_indel_rates(self):

        if not (self.hierarchy and self.args.indels):
            print("error, indel model only allowed with hierarchical mode, remove --noHierarchy")
            print("exiting")
            exit()

        if (not self.args.insertionRate) or (not self.args.deletionRate):
            print("error, indel models must specify --insertionRate and --deletionRate")
            print("exiting")
            exit()

        def parse_indel_rates(cli_parameters):
            model = cli_parameters[0]
            parameter = [float(x) for x in cli_parameters[1:]]

            if not model.upper() in ["GAMMA", "CONSTANT"]:
                print("error - only indel models allowed are GAMMA or CONSTANT, exiting.")
                exit()

            if model == "CONSTANT":
                generator = iter(np.full(self.ref_len, parameter[0]))
                generator = chain(generator, (parameter[0] for _ in count()))
            
            else:
                if len(parameter) != 2 or parameter[0] <= 0 or parameter[1] <= 0:
                    print("Gamma distribution requires 2 parameters, alpha and beta (shape and rate). Exiting")
                    exit()

                if parameter[0]/parameter[1] > 0.25:
                    print("Warning: indels are meant to be rare but the average per nucleotide indel rate is the ratio of parameters:")
                    print("Ratio: " + str(parameter[0]/parameter[1]))

                generator = iter(np.random.gamma(parameter[0], 1./ parameter[1], size=self.ref_len))
                generator = chain(generator, (np.random.gamma(parameter[0], 1./ parameter[1]) for _ in count()))
            
            return generator

        insertionGenerator = parse_indel_rates(self.args.insertionRate)
        deletionGenerator = parse_indel_rates(self.args.deletionRate)

        return (insertionGenerator, deletionGenerator)


    def init_indel_lengths(self):

        if not (self.hierarchy and self.args.indels):
            print("error, indel model only allowed with hierarchical mode, remove --noHierarchy")
            print("exiting")
            exit()

        # generate a stream of random indel lengths to be passed to the genomeTree_hierarchical
        def parse_indel_args(argument_name):
            indel_data = getattr(self.args, argument_name)

            if not indel_data:
                print("error, indel model requires --insertionLength and --deletionLength")
                print("exiting")
                exit()

            distribution = indel_data[0].upper()
            parameters = indel_data[1:]

            if distribution == "GEOMETRIC":

                if len(parameters) != 1:
                    print(f"{distribution} requires exactly 1 parameter, p.")
                    print(f"This model had {len(parameters)} parameter values.") 
                    exit()

                p = float(parameters[0])

                return (np.random.geometric(p) for _ in count())

            elif distribution == "NEGBINOMIAL":

                if len(parameters) != 2:
                    print(f"{distribution} requires exactly 2 parameters, r p.")
                    print(f"This model had {len(parameters)} parameter values.") 
                    exit()

                r = int(parameters[0])
                p = float(parameters[1])  

                return (1 + np.random.negative_binomial(r, p) for _ in count())

            elif distribution == "ZIPF":

                if len(parameters) != 1:
                    print(f"{distribution} requires exactly 1 parameter, a.")
                    print(f"This model had {len(parameters)} parameter values.") 
                    exit()

                a = float(parameters[0])

                if a <= 1.0:
                    print(f"Parameter a in Zipf distribution must be > 1 whereas, value is: {a}")
                    exit()

                return (np.random.zipf(a) for _ in count())

            elif distribution == "LAVALETTE":

                if len(parameters) != 2:
                    print(f"{distribution} requires exactly 2 parameters, M a.")
                    print(f"This model had {len(parameters)} parameter values.") 
                    exit()

                M = int(parameters[0])
                a = float(parameters[1])

                probs = [(float(i) * float(M)/ float(M - i + 1)) ** (-1 * a) for i in range(1, M+1)]

                total = sum(probs)

                probs = [x/ total for x in probs]

                return (np.random.choice(range(1, M+1), p=probs) for _ in count())

            elif distribution == "DISCRETE":

                parameters = [float(x) for x in parameters]

                if len(parameters) < 1:
                    print("DISCRETE option requires parameters p_1 p_2, ...etc.")
                    exit()

                if abs(sum(parameters) - 1) > 0.00001:
                    print(f"Sum of discrete probabilities adds to {sum(parameters)}, normalising so that the total is 1.")
                    parameters = [x/sum(parameters) for x in parameters]
                    print(f"New values are: {parameters}")

                return (np.random.choice(range(1, len(parameters)+1), p=parameters) for _ in count())

            else:
                print(f"This distribution, {distribution}, is not yet implemented, choose one of ")
                print("GEOMETRIC, NEGBINOMIAL, ZIPF, LAVALETTE, DISCRETE, followed by comma separated parameter values.")
                print("Use phastSim --help for more information.")
                exit()

        insertionLength = parse_indel_args("insertionLength")
        deletionLength = parse_indel_args("deletionLength")

        return insertionLength, deletionLength
        

    def init_insertion_frequencies(self, ref):
        # the model for insertions will be the same as for generating a random root genome, using random choices from
        # some genome frequency distribution

        # in this situation we need to take a look at the reference genome 
        # and count the occurrences of each symbol (nucleotide or codon)
        if self.args.indels and self.args.rootGenomeLength == 0 and self.args.rootGenomeFrequencies == [0.0]:
            if self.args.codon:
                print("Computing the frequencies of each codon assuming the whole genome is in frame 0.")
                codons = ["".join(z) for z in product("ACGT", "ACGT", "ACGT")]
                ref_as_codons = (ref[n:n+3] for n in range(0,len(ref),3))
                n_codons = int(len(ref)/3)
                c = Counter(ref_as_codons)
                self.args.rootGenomeFrequencies = [float(c[x])/float(n_codons) for x in codons]
                print("Codon frequencies: ", self.args.rootGenomeFrequencies)

            else:
                print("Counting occurences of each nondegenerate nucleotide in the genome")
                c = Counter(ref)
                self.args.rootGenomeFrequencies = [float(c[x])/float(self.ref_len) for x in "ACGT"]
                print("ACGT frequencies: ", self.args.rootGenomeFrequencies)

        if not self.args.codon:
            rootGenomeFrequencies = self.args.rootGenomeFrequencies


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

            return rootGenomeFrequencies


        # generate frequencies under codon model
        else:
            rootGenomeFrequencies = self.args.rootGenomeFrequencies

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

                    
            return rootGenomeFrequencies


class GenomeTree_hierarchical:
    def __init__(self, nCodons, codon, ref, gammaRates, omegas, mutMatrix, hyperCategories, hyperMutRates, 
                indels, insertionRate, insertionLength, insertionFrequencies, deletionRate, deletionLength, scale, infoFile, verbose, noNorm):

        self.codon = codon
        self.ref = ref
        self.refList = list(ref)
        self.gammaRates = gammaRates
        self.omegas = omegas
        self.mutMatrix = mutMatrix

        self.hyperCategories = hyperCategories
        self.hyperMutRates = hyperMutRates
        self.infoFile = infoFile
        self.verbose = verbose

        const = Constants()
        self.stopCodons = const.stopCodons
        self.alleles = const.alleles
        self.allelesList = const.allelesList
        self.nAlleles = const.nAlleles

        self.alRange = range(const.nAlleles)
        self.range9 = range(9)
        self.indels = indels
        self.insertionRate = insertionRate
        self.insertionLength = insertionLength
        self.insertionFrequencies = insertionFrequencies
        self.deletionRate = deletionRate
        self.deletionLength = deletionLength
        self.scale = scale
        self.noNorm = noNorm

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

        if self.indels:
            # keep track of the number of insertions so far
            self.n_insertions = 0
            # we will add an extra 'dummy' symbol at the start of the genome. 
            # insertions can happen at this extra 0th site. 
            self.nTerminalNodes +=1

    def __str__(self):
        return str(self.genomeRoot)

    def populate_genome_tree(self):
        # public function - just calls populateGenomeTree with the correct parameters. 
        pos_left = (-1 if self.indels else 0)
        pos_right = (self.nTerminalNodes - 1) if self.indels else self.nTerminalNodes
        self.populateGenomeTree(self.genomeRoot, self.ref, pos_left, pos_right - 1)

    def populateGenomeTree(self, node, ref, pos_left, pos_right):
        if pos_left == pos_right:
            # these are terminal nodes. Only these nodes need information regarding specific alleles and rates.
            node.genomePos = pos_left
            pos = pos_left
            node.isTerminal = True
            # the corresponding node in the level 0 tree; in this case, it's itself.
            # Only terminal nodes in the level 0 tree need to contain info regarding all mutation rates.
            node.refNode = node

            if self.indels:
                # set the insertion rate of the node, and the insertion index (i.e. which number insert - 0 meaning the node was on the reference).
                node.insertionPos = self.n_insertions
                node.insertionRate = next(self.insertionRate)

                # the first node is a 'dummy' node with allele -1, representing a blank character. 
                # An insertion can happen at this node, but nothing else. 
                if (pos == -1):
                    node.deletionRate = 0.0
                    node.rate = 0.0 
                    node.allele = -1
                    node.rates = np.zeros(1)
                    return
            
                # set the deletion rates of the (non-lead-character) node
                node.deletionRate = next(self.deletionRate)


            # the current allele at this terminal node. It starts as the reference allele at the considered position.
            if self.codon:
                node.allele = self.codonAlleles[ref[pos * 3:(pos + 1) * 3]]
                node.omega = next(self.omegas)
                if node.allele in self.stopCodons:
                    node.omega = 0.0
                    if (pos < self.nTerminalNodes - 1):
                        print(f"Warning: stop codon in the middle of the reference {pos * 3 + 1}"
                          f" {ref[pos * 3:(pos + 1) * 3]}."
                          f"Continuing simulations but setting omega=0 for this position.")


                if self.infoFile!=None:
                    if self.indels:
                        self.infoFile.write(str(node.insertionPos) + "\t")
                    self.infoFile.write(str(pos * 3 + 1) + "-" + str(pos * 3 + 3) + "\t" + str(node.omega) + "\t")
                node.rates = {}
                node.hyper = {}
                node.gammaRates = [0.0,0.0,0.0]
                node.hyperCategories = [0,0,0]
                #Nicola: create initial mega-array of zeros(10) and pass just index to it?
                if self.indels and node.insertionPos:
                    node.rates[node.allele] = np.zeros(10)
                else:
                    node.rates[node.allele]= self.allRates[pos]

                indices = self.codonIndices[node.allele]
                node.rate = 0.0
                range3=self.range3
                
                for i2 in range3:

                    node.hyperCategories[i2] = next(self.hyperCategories)
                    node.gammaRates[i2] = next(self.gammaRates)
                    if self.infoFile!=None:
                        self.infoFile.write(str(node.gammaRates[i2]) + "\t" + str(node.hyperCategories[i2]) + "\t")
                    nuc1 = indices[i2]
                    for i3 in range3:
                        nuc2 = (nuc1 + i3 + 1) % 4
                        if self.isNonsynom[node.allele, i2, i3]:
                            if self.isIntoStop[node.allele, i2, i3]:
                                node.rates[node.allele][i2 * 3 + i3] = 0.0
                            else:
                                node.rates[node.allele][i2 * 3 + i3] = node.omega * self.mutMatrix[nuc1][nuc2] * node.gammaRates[i2]
                        else:
                            node.rates[node.allele][i2 * 3 + i3] = self.mutMatrix[nuc1][nuc2] * node.gammaRates[i2]
                    if node.hyperCategories[i2] > 0:
                        # site with hypermutation. Sample a random allele i to hypermutate,
                        # and a random destination allele j;
                        # the mutation rate from i to j at the given position is then enhanced.
                        i = np.random.choice(self.nAlleles)
                        j = np.random.choice(self.nAlleles - 1)
                        node.hyper[i2] = (i, j)
                        if i == nuc1:
                            node.rates[node.allele][i2 * 3 + j] *= self.hyperMutRates[node.hyperCategories[i2] - 1]
                        if self.infoFile!=None:
                            self.infoFile.write(self.allelesList[i] + "\t" + self.allelesList[(i + j + 1) % 4] + "\t")
                    else:
                        if self.infoFile!=None:
                            self.infoFile.write(".\t" + ".\t")
                    for i3 in range3:
                        node.rate += node.rates[node.allele][i2 * 3 + i3]
                    node.rates[node.allele][9] = node.rate
                if self.infoFile!=None:
                    if self.indels:
                        self.infoFile.write(str(node.insertionRate) + "\t" + str(node.deletionRate))
                    self.infoFile.write("\n")

            else:
                node.allele = self.alleles[ref[pos]]
                # evolutionary categories
                node.gammaRate = next(self.gammaRates)
                nodeHyper = next(self.hyperCategories)
                node.hyperCategories = nodeHyper
                if self.infoFile!=None:
                    if self.indels:
                        self.infoFile.write(str(node.insertionPos) + "\t")
                    self.infoFile.write(str(pos + 1) + "\t" + str(node.gammaRate) + "\t" + str(nodeHyper) + "\t")
                # the mutation rates for the considered site
                if self.indels and node.insertionPos:
                    node.rates = self.mutMatrix.copy()
                else:
                    node.rates = self.mutMatrixRepeats[pos]
                node.rates *= node.gammaRate
                # node.rates=np.copy(mutMatrix)*categoryRates[node.category]
                # NICOLA: THIS COULD MAYBE BE MADE MORE EFFICIENT.
                # FOR EXAMPLE, ONE COULD INITIALIZE RATES ONLY WHEN SAMPLING AT A LEAF,
                # AND OTHERWISE INITIALIZE ONLY THE TOTAL RATE.
                # ALSO, INSTATED OF COPYING AND STORING MULTIPLE MUTATION MATRICES, ONE COULD LINK TO THE SAME ONE
                # FOR SITES WITH THE SAME RATES.
                # now sample which alleles are affected by hypermutation and store info inscale positions and extra vectors
                if nodeHyper > 0:
                    # site with hypermutation. Sample a random allele i to hypermutate,
                    # and a random destination allele j;
                    # the mutation rate from i to j at the given position is then enhanced.
                    i = np.random.choice(self.nAlleles)
                    j = np.random.choice(self.nAlleles - 1)
                    j = (i + j + 1) % self.nAlleles
                    node.rates[i][i] -= node.rates[i][j] * (self.hyperMutRates[nodeHyper - 1] - 1.0)
                    node.rates[i][j] *= self.hyperMutRates[nodeHyper - 1]
                    if self.infoFile!=None:
                        self.infoFile.write(self.allelesList[i] + "\t" + self.allelesList[j])
                        if self.indels:
                            self.infoFile.write("\t" + str(node.insertionRate) + "\t" + str(node.deletionRate))
                        self.infoFile.write("\n")
                else:
                    if self.infoFile!=None:
                        self.infoFile.write(".\t" + ".")
                        if self.indels:
                            self.infoFile.write("\t" + str(node.insertionRate) + "\t" + str(node.deletionRate))
                        self.infoFile.write("\n")
                # total mutation rate at the node
                node.rate -= node.rates[node.allele][node.allele]


        else:
            # split the considered part of the genome in two, assign each half to both children,
            # then call the populate function iteratively on the children.
            middle = pos_left + int((pos_right - pos_left) / 2)
            firstChild = genomeNode(level=node.level)
            secondChild = genomeNode(level=node.level)
            node.belowNodes = [firstChild, secondChild]
            self.populateGenomeTree(node=firstChild, ref=ref, pos_left=pos_left, pos_right=middle)
            self.populateGenomeTree(node=secondChild, ref=ref, pos_left=middle+1, pos_right=pos_right)
            node.rate = firstChild.rate + secondChild.rate

    def check_start_stop_codons(self):

        if self.codon:
            # if first codon is a start codon, or last is a stop codon, don't allow them to evolve
            # if ref[0:3]=="ATG" or ref[0:3]=="atg" or ref[0:3]=="AUG" or ref[0:3]=="aug":
            if self.ref[0:3] in ["ATG", "AUG"]:
                # set gammaRates for the first node to 0
                parent = self.genomeRoot
                child = parent.belowNodes[0]
                while (not child.isTerminal):
                    parent = child
                    child = parent.belowNodes[0]

                if self.indels:
                    # the 1st node is a dummy node when using indels so the 2nd node refers to the 1st codon
                    parent.belowNodes[1].gammaRates = [0.0, 0.0, 0.0]
                else:
                    parent.belowNodes[0].gammaRates = [0.0, 0.0, 0.0]


            stopCodons = ["TAA", "TGA", "TAG"]
            if self.ref[-2:] in stopCodons:
                node = self.genomeRoot
                while (not node.isTerminal):
                    node = node.belowNodes[1]
                
                node.omega = 0.0


    def normalize_rates(self):

        # When normalizing, branch lengths are in number of substitutions per nucleotide,
        # even though we might be simulating a codon model.
        if self.noNorm:
            self.norm = 1.0/self.scale
            print("\n Not normalizing mutation rates, as required by user. ")

        else:
            # We rescale by the input normalization factor, this is the same as rescaling all
            # the branch lengths by this rescaling factor
            norm = self.genomeRoot.rate / (len(self.ref) * self.scale)
            self.norm = norm
            print("\n Total cumulative substitution rate per site before normalization: " + str(norm))
        
        self.normalizeRates(self.genomeRoot)


    def normalizeRates(self, rootNode):
        # This is an internal function implementation that can be reused elsewhere (e.g. when creating indel inserts).

        norm = self.norm
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
        normalize(node=rootNode, norm=norm)

        if self.indels:
            # if there are indels then some further calculation is require to convert the rates from 
            # substitution rates to overall rates.
            

            def include_indels(node, scale):
                
                if node.isTerminal:
                    node.insertionRate *= scale
                    node.deletionRate *= scale
                    node.rate += (node.insertionRate + node.deletionRate)
                
                else:
                    node.rate = include_indels(node.belowNodes[0], scale) + include_indels(node.belowNodes[1], scale)

                return node.rate

            include_indels(node=rootNode, scale=self.scale)

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


    def deleteNodes(self, rand, node, remaining_deletions, level):
        # This function should be called once findPos has sampled a deletion at a particular node. 
        # It will recursively traverse the genome tree starting from the root, 'deleting' nodes by 
        # setting their rates to 0 and their allele as -1, creating the node at a new level if necessary. 
        # The function then returns the total number of deleted nodes. 
        # The function avoids traversing the whole genome tree by exiting early if the total number of deletions
        # has been carried out, or if it is on a branch to the left of the site being deleted. 

        if node.isTerminal:
            if node.allele == -1:
                # nodes with this allele are 'empty' symbols and cannot/ should not count as being deleted
                return 0, ""

            deleted_string = (self.codonAllelesList[node.allele] if self.codon else self.allelesList[node.allele])
            node.allele = -1
            node.rate = 0.0
            node.insertionRate = 0.0
            node.deletionRate = 0.0
            return 1, deleted_string
        
        else:
            # nodes need to be created if on a different level, otherwise can directly change nodes
            child0, child1 = node.belowNodes[0], node.belowNodes[1]
            rate0, rate1 = child0.rate, child1.rate
            newRate0, newRate1 = rate0, rate1
            newChild0, newChild1 = child0, child1

            # deal with the left child
            deleted_leaves_left = 0
            deleted_string_left = ""
            if rand > rate0 or remaining_deletions == 0:
                # under these conditions we can short circuit and ignore this whole branch
                pass

            else:
                if child0.level < level or not (child0.isTerminal and child0.allele == -1):
                    newChild0 = genomeNode(level=level)
                    newChild0.isTerminal = child0.isTerminal
                    newChild0.rate = child0.rate

                    if child0.isTerminal:
                        newChild0.refNode = child0.refNode
                        newChild0.allele = child0.allele
                    else:
                        newChild0.belowNodes = list(child0.belowNodes)

                
                deleted_leaves_left, deleted_string_left = self.deleteNodes(rand, newChild0, remaining_deletions, level)
                newRate0 = newChild0.rate

            # deal with right child but with slightly different exit conditions
            deleted_leaves_right = 0
            deleted_string_right = ""
            if rand > rate0 + rate1 or remaining_deletions - deleted_leaves_left == 0:
                # can short circuit this child under slightly different conditions
                pass

            else:
                if child1.level < level or not (child1.isTerminal and child1.allele == -1):
                    newChild1 = genomeNode(level=level)
                    newChild1.isTerminal = child1.isTerminal
                    newChild1.rate = child1.rate

                    if child1.isTerminal:
                        newChild1.refNode = child1.refNode
                        newChild1.allele = child1.allele

                    else:
                        newChild1.belowNodes = list(child1.belowNodes)
       

                deleted_leaves_right, deleted_string_right = self.deleteNodes(rand - rate0, newChild1, remaining_deletions - deleted_leaves_left, level)
                newRate1 = newChild1.rate

            node.rate = newRate0 + newRate1
            node.belowNodes = [newChild0, newChild1]

            return deleted_leaves_left + deleted_leaves_right, deleted_string_left + deleted_string_right


    def sampleInsertion(self, node, level):
        

        insertion_length = next(self.insertionLength)

        # increment the total insertion counter
        self.n_insertions += 1
        
        insertion = None
        genomeSeq = ""
        if self.codon:
            insertion = np.random.choice(self.codonAllelesList, insertion_length, self.insertionFrequencies)
        else:
            insertion = np.random.choice(self.allelesList, insertion_length, self.insertionFrequencies)
            # turn the nucleotide sequence into a genomeTree
        
        genomeSeq = "".join(insertion)
        # generate the genome tree of this insert from (node, genomeSeq)
        new_subtree_root = genomeNode(level=level)
        self.populateGenomeTree(new_subtree_root, genomeSeq, 0, insertion_length-1)
        self.normalizeRates(new_subtree_root)

        # build a bigger tree of the form node = (node_copied, subTree). 
        node_copied = deepcopy(node)

        # delete all of the old properties of the node
        for prop in node_copied.__dict__:
            delattr(node, prop)

        # set the new properties of the node
        node.belowNodes = [node_copied, new_subtree_root]
        node.isTerminal = False       
        node.level = node_copied.level
        node.rate = node_copied.rate + new_subtree_root.rate
        
        pos = node_copied.refNode.genomePos
        if self.codon:
            pos = 3*pos + 2

        mutEvent = mutation(
                    mType=mType.INS, 
                    genomePos=pos, 
                    insertionPos=node_copied.refNode.insertionPos, 
                    source="",
                    target=genomeSeq, 
                    index=self.n_insertions)

        return mutEvent

    def findPos(self, rand, parentGenomeNode, level):
        # find position to mutate along the genome, and update temporary genome tree structure as you go
        if parentGenomeNode.isTerminal:
            # reached a terminal node, now sample the mutation event at the position and update all rates
            node = parentGenomeNode.refNode
            a = parentGenomeNode.allele
            
            # we'll choose an indel if rand is less than a threshold determined by the insertion and deletion rates for this site
            if self.indels:

                # deletions
                if rand < node.deletionRate:

                    attempted_deletion_length = next(self.deletionLength)
                    pos = node.genomePos
                    deleted_string = ((self.codonAllelesList[node.allele] if self.codon else self.allelesList[node.allele])
                                         + "?" * (attempted_deletion_length - 1) * (3 if self.codon else 1))

                    mutEvent = mutation(
                                mType=mType.DEL, 
                                genomePos=(3 if self.codon else 1) * pos, 
                                insertionPos=node.insertionPos, 
                                source=deleted_string, 
                                target="", 
                                length=attempted_deletion_length)

                    # In the tidy up function (self.deleteNodes), the length (and also the deleted string) of this mutEvent may be changed 
                    # if we are in an edge case where we delete at a site very close to the end of the genome. 
                    # The rates for all nodes are also updated there. 
                    return mutEvent
                
                # else insertion:
                # rescale rand so that it is uniformly between 0 and insertionRate
                rand -= node.deletionRate
                if rand < node.insertionRate:
                
                    # sampleInsertion will mutate the parentGenomeNode
                    mutEvent = self.sampleInsertion(parentGenomeNode, level)

                    return mutEvent

                # rescale rand so that it is uniformly between 0 and the substitution rate := node.rate - (node.insertionRate + node.deletionRate)
                rand -= node.insertionRate
            
            # otherwise we are doing a substitution and can proceed as normal
            if self.codon:
                j = self.sampleMutationCodon(rates=node.rates[a], rand=rand)
                indices = self.codonIndices[a]
                i2 = int(j / 3)
                i3 = j % 3
                newindices = list(indices)
                newindices[i2] = (newindices[i2] + i3 + 1) % 4
                parentGenomeNode.allele = self.codonIndices2[newindices[0], newindices[1], newindices[2]]
                mutEvent = mutation(
                            mType=mType.SUB, 
                            genomePos=node.genomePos * 3 + i2, 
                            source=self.allelesList[indices[i2]], 
                            target=self.allelesList[newindices[i2]])

                if self.indels:
                    mutEvent.insertionPos = node.insertionPos

                if self.verbose:
                    print(f"Mutation from {a} {self.codonAllelesList[a]} "
                          f"to {parentGenomeNode.allele} {self.codonAllelesList[parentGenomeNode.allele]}"
                          f" , position {mutEvent.genomePos} category rate {node.gammaRates[i2]}"
                          f" hyperCat {node.hyperCategories[i2]}"
                          f" omega {node.omega}"
                          f" old rate {node.rates[a][9]} old rates:")
                    print(node.rates[a])
                if not (parentGenomeNode.allele in node.rates):
                    node.rates[parentGenomeNode.allele] = np.zeros(10)
                    indices = self.codonIndices[parentGenomeNode.allele]
                    parentGenomeNode.rate = 0.0
                    range3=self.range3
                    for i2 in range3:
                        nuc1 = indices[i2]
                        for i3 in range3:
                            nuc2 = (nuc1 + i3 + 1) % 4
                            if self.isNonsynom[parentGenomeNode.allele, i2, i3]:
                                if self.isIntoStop[parentGenomeNode.allele, i2, i3]:
                                    node.rates[parentGenomeNode.allele][i2 * 3 + i3] = 0.0
                                else:
                                    node.rates[parentGenomeNode.allele][i2 * 3 + i3] = node.omega * \
                                                                                       self.mutMatrix[nuc1][nuc2] * \
                                                                                       node.gammaRates[i2] / self.norm
                            else:
                                node.rates[parentGenomeNode.allele][i2 * 3 + i3] = self.mutMatrix[nuc1][nuc2] * \
                                                                                   node.gammaRates[i2] / self.norm
                        if node.hyperCategories[i2] > 0:
                            if node.hyper[i2][0] == nuc1:
                                node.rates[parentGenomeNode.allele][i2 * 3 + node.hyper[i2][1]] *= self.hyperMutRates[
                                    node.hyperCategories[i2] - 1]
                        for i3 in range3:
                            node.rates[parentGenomeNode.allele][9] += node.rates[parentGenomeNode.allele][i2 * 3 + i3]
                parentGenomeNode.rate = node.rates[parentGenomeNode.allele][9]
                if self.indels:
                    parentGenomeNode.rate += (node.insertionRate + node.deletionRate)
                if self.verbose:
                    print(f" new rate {parentGenomeNode.rate} all rates:")
                    print(node.rates[parentGenomeNode.allele])
            else:
                j = self.sampleMutation(a, node.rates[a], rand)

                mutEvent = mutation(mType=mType.SUB, genomePos=node.genomePos, source=self.allelesList[a], target=self.allelesList[j])
                if self.verbose:
                    print(f"Mutation from {a} to {j},"
                          f" position {node.genomePos}"
                          f" category rate {node.gammaRate}"
                          f" hyperCat {node.hyperCategories}"
                          f" old rate {parentGenomeNode.rate} old rates:")
                    print(node.rates)
                parentGenomeNode.rate = -node.rates[j, j]
                if self.indels:
                    parentGenomeNode.rate += (node.insertionRate + node.deletionRate)
                    mutEvent.insertionPos = node.insertionPos
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

            # if the child we are moving into is not on the same level, but is below,
            # then create a new child at the same level.
            # this is because the rate of the child will be inevitably changed by the mutation event,
            # and we don't want to change the mutation rates for the parent phylogenetic node.
            if child.level < level:
                newChild = genomeNode(level=level)  # upNode=parentGenomeNode
                parentGenomeNode.belowNodes[childI] = newChild
                newChild.isTerminal = child.isTerminal
                newChild.rate = child.rate
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


    def applyMutation(self, parentGenomeNode, level, mutEvent):
        # find position to mutate along the genome, and update temporary genome tree structure as you go.
        # this one is only used for mutations that are forced by the user
        node = parentGenomeNode.refNode
        a = parentGenomeNode.allele

        if parentGenomeNode.isTerminal:
            # reached a terminal node, now force the mutation event at the position and update all rates
            j=self.alleles[mutEvent.target]

            if self.verbose:
                print(f"Forced Mutation from {a} to {j},"
                        f" position {node.genomePos}"
                        f" category rate {node.gammaRate}"
                        f" hyperCat {node.hyperCategories}"
                        f" old rate {parentGenomeNode.rate} old rates:")
                print(node.rates)
                print(" new rate " + str(parentGenomeNode.rate) + " all rates:")
                print(node.rates)
            parentGenomeNode.allele = j

        else:
            # still at an internal genome node.
            # choose which of the two children genome nodes to move into
            pos=mutEvent.genomePos
            child1Pos=parentGenomeNode.belowNodes[0].refNode.genomePos

            if pos>=child1Pos[0] and pos<=child1Pos[1]:
                parentGenomeNode.rate = parentGenomeNode.belowNodes[1].rate
                child = parentGenomeNode.belowNodes[0]
                childI = 0
            else:
                parentGenomeNode.rate = parentGenomeNode.belowNodes[0].rate
                child = parentGenomeNode.belowNodes[1]
                childI = 1

            # if the child we are moving into is not on the same level, but is above,
            # then create a new child at the same level.
            # this is because the rate of the child will be inevitably changed by the mutation event,
            # and we don't want to change the mutation rates for the parent phylogenetic node.
            if child.level < level:
                newChild = genomeNode(level=level)  
                parentGenomeNode.belowNodes[childI] = newChild
                newChild.isTerminal = child.isTerminal
                newChild.refNode=child.refNode
                if child.isTerminal:
                    newChild.refNode = child.refNode
                    newChild.allele = child.allele
                else:
                    newChild.belowNodes = list(child.belowNodes)

                self.applyMutation(newChild, level, mutEvent)

                parentGenomeNode.rate += newChild.rate
            else:
                # in this case the child is already on the same level, so no need to create another one, just update its mutation rate.
                self.applyMutation(child, level, mutEvent)
                parentGenomeNode.rate += child.rate


    def mutateBranchETEhierarchy(self, childNode, parentGenomeNode, level, createNewick, preMutationsBranches):
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

        #Go through the input mutation events at the current branch and apply them to the genome tree.
        if childNode.name in preMutationsBranches:
            newGenomeNode = genomeNode(level=level)
            newGenomeNode.belowNodes = list(parentGenomeNode.belowNodes)
            for m in preMutationsBranches[childNode.name]:
                pos=m[0]
                fromA=m[1]
                toA=m[2]
                print("applying mutation")
                mutEvent=mutation(mType.SUB, pos-1,source=fromA,target=toA)
                print(mutEvent)
                self.applyMutation(newGenomeNode, level, mutEvent)
                childNode.mutations.append(mutEvent)
                if createNewick:
                    childNode.mutAnnotation.append(str(mutEvent))
            rate = newGenomeNode.rate

        # Sample new mutation event with Gillespie algorithm
        currTime += np.random.exponential(scale=1.0 / rate)
        if self.verbose:
            print(f"\n Node {childNode.name} BLen: {bLen} first sampled time: {currTime}; mutation rate: {rate}")
        # for the first mutation event at this node, create a new root genome node of the appropriate level.
        # otherwise, use the one you already have.
        
        if currTime < bLen and (not (childNode.name in preMutationsBranches)):
        #if currTime < bLen:
            newGenomeNode = genomeNode(level=level)
            newGenomeNode.belowNodes = list(parentGenomeNode.belowNodes)
        
        elif (not (childNode.name in preMutationsBranches)):
            newGenomeNode = parentGenomeNode
        
        while currTime < bLen:
            # Now, sample which type of mutation event it is (from which nucleotide to which nucleotide)
            rand = np.random.random() * rate
            if self.verbose:
                print("Selecting new mutation event. Rate " + str(rate) + " random value " + str(rand))
            mutEvent = self.findPos(rand, newGenomeNode, level)
            
            # if the mutation event was a deletion then we need to go through the genome tree and delete the right number of nodes
            if mutEvent.mType == mType.DEL:
                # The mutEvent produced by findPos may not actually have managed to produce a deletion of the correct length
                # (this may happen if a long deletion is sampled near the end of the genome and runs out of symbols).
                # To deal with this we update the mutEvent once deleteNodes has finished its tidy up - returning the true number
                # of deleted nodes. 
                length_deleted, string_deleted = self.deleteNodes(rand, newGenomeNode, mutEvent.length, level)
                mutEvent.length = length_deleted
                mutEvent.source = string_deleted

                if self.verbose:
                    print(f"Deletion occurred: {mutEvent.source} at coordinates ({mutEvent.insertionPos},{mutEvent.genomePos})")

            childNode.mutations.append(mutEvent)
            if createNewick:
                childNode.mutAnnotation.append(str(mutEvent))
            rate = newGenomeNode.rate
            if self.verbose:
                print(f"New total rate {rate}")
            currTime += np.random.exponential(scale=1.0 / rate)
            if self.verbose:
                print(f"new time {currTime}, rate {rate} mutation events:")
                print([str(m) for m in childNode.mutations])

        if self.verbose:
            print("mutations at the end:")
            print([str(m) for m in childNode.mutations])
        # now mutate children of the current node, calling this function recursively on the node children.
        for c in childNode.children:
            self.mutateBranchETEhierarchy(c, newGenomeNode, level + 1, createNewick, preMutationsBranches)


    def getAbsoluteCoords(self, mut, mutDict, insertionDict):
        """
        refPos - the position on the reference that this mut appears on. 
        offset - if the mutation is nested inside another mutation, or not at the start of a series of mutations, 
        then it may have a non-zero offset, referring to how far along the mutation is inside another already existing mutation. 
        """
        if mut.insertionPos:
            
            refPos = insertionDict[mut.insertionPos]
            for offset, coords in enumerate(mutDict[refPos][2]):
                if (mut.insertionPos, mut.genomePos) == coords:
                    return refPos, offset

        else:
            return mut.genomePos, 0


    def wrapWriterWithMutationDict(writingFunction):
        """
        This function wraps a writing function, by keeping track of a mutDict and insertionDict; the writing function is there
        just to write the output to the file in the correct format. 

        I've done this because the same wrapper with a mutDict and insertionDict is used in multiple places. 

        This whole function is overly complicated and potentially could be quite slow. 
        The aim is to 'collapse' all of the (possibly overlapping / nested) indels and substitutions
        into a concise format. 

        Any advice on how to avoid all this complicated work (either by simplifying the code or deciding on
        a different but still unambiguous output format) is appreciated!

        mutDict = {
            position: (original nuc, mutated nucs, list of the coordinates of the new nuc),

            1: (A, -AATTTG, [(0,1), (2,0), (2,1),(5, 0), (5,1), (5,2), (2,2)])

            123: (ATT, "", [])

            456: (A, G, [(0, 456)])
        }

        insertionDict = {
            2: 1
            5: 1
        }

        (2, 1) -> refPos = 1, offset = 2

        (0, 123) -> refPos = 0, offset = 123
        
        mutation(source = "", target = T, insertionPos = 0, genomePos = 123, [index = 6])

        mutation(source = AAAAA, target = "", insertionPos = 0, genomePos = 234)

        """

        def wrappedFunction(self, node, file, mutDict, insertionDict, **kwargs):

            # flatten the mutDict + node.mutations into a more concise format, and save indels in insertionDict
            for m in node.mutations:

                # find the 'absolute' coordinate of the mutation i.e. where it lives relative to the 
                # already existing mutations in the mutDict. 
                refPos, offset = self.getAbsoluteCoords(m, mutDict, insertionDict)

                if m.mType == mType.SUB:

                    if refPos in mutDict:
                        target = mutDict[refPos][1]
                        mutDict[refPos][1] = target[:offset] + m.target + target[offset+1:]

                        # if the substitution has exactly reversed a previous substitution then we can delete it
                        if mutDict[refPos][1] == self.ref[refPos]:
                            del mutDict[refPos]

                    else:
                        mutDict[refPos] = [self.ref[refPos], m.target, [(m.insertionPos, m.genomePos)]]

                elif m.mType == mType.INS:

                    # add the insertion to the insertion dictionary
                    insertionDict[m.index] = refPos

                    if refPos in mutDict:
                        # record the indices i.e. the relative genome coordinates of the existing data
                        indices = mutDict[refPos][2]
                        target = mutDict[refPos][1]
                        

                        # the new data is comprised of: original, new insertion, new indices of genome
                        mutDict[refPos] = [mutDict[refPos][0], 
                                        target[:offset+1] + m.target + target[offset+1:],
                                        indices[:offset+1] + [(m.index, i) for i in range(len(m.target))] + indices[offset+1:]]

                    else:
                        # deal with an edge case where we are at the -1 genome position i.e. right at the start of the genome
                        first_symbol = (self.ref[refPos] if refPos != -1 else "")
                        mutDict[refPos] = [first_symbol,
                                        first_symbol + m.target,
                                        [(m.insertionPos, m.genomePos)] + [(m.index, i) for i in range(len(m.target))]]

                # deletion - this is the most difficult case
                else:

                    # need to deal with the deletion one symbol at a time
                    deletedChars = 0
                    m.deletedPositions = []
                    
                    startedInMiddleOfDictionaryItem = False
                    if offset > 0:
                        startedInMiddleOfDictionaryItem = True

                    # are we appending "-" characters to the end of an existing dictionary item or not?
                    appending = False

                    # dictPos - the dictionary key in mutDict that we are reading from / writing to
                    dictPos = refPos

                    while (deletedChars < len(m.source)):

                        # we may be appending characters to the current deletion one at a time
                        if appending:
                            
                            # if the next character is already in the dictionary, we must stop appending 
                            if refPos in mutDict:
                            
                                appending = False
                                
                                # delete a character if it is non-blank
                                if mutDict[dictPos][1][offset] != "-":
                                    mutDict[dictPos][1] = mutDict[dictPos][1][:offset] + "-" + mutDict[dictPos][1][offset+1:]
                                    deletedChars += 1
                                    m.deletedPositions.append(mutDict[dictPos][2][offset])
                            
                            # otherwise we can continue appending to the deletion 
                            else:
                                
                                mutDict[dictPos][0] += self.ref[refPos]
                                mutDict[dictPos][1] += "-"
                                mutDict[dictPos][2] += [(0, refPos)]
                                deletedChars += 1
                                m.deletedPositions.append((0, refPos))
                        
                        # in this case we are working on an existing dictionary item, or have just finished doing so
                        else:
                            
                            # a position not ever seen before, we must add it to the dictionary and 
                            # can append to it
                            if not dictPos in mutDict:
                                
                                mutDict[dictPos] = [self.ref[dictPos], "-", [(0, dictPos)]]
                                deletedChars += 1
                                m.deletedPositions.append((0, dictPos))
                                appending = True
                                
                            # delete a character if it is non-blank
                            else: 
                                if mutDict[dictPos][1][offset] != "-":
                                    mutDict[dictPos][1] = mutDict[dictPos][1][:offset] + "-" + mutDict[dictPos][1][offset+1:]
                                    deletedChars += 1
                                    m.deletedPositions.append(mutDict[dictPos][2][offset])
                        
                        
                        # now move to the next character
                        # increment the reference position if we pass it
                        if (mutDict[dictPos][2][offset] == (0, refPos)):
                            refPos += 1  
                            
                        # we will change to a new dictionary item if we are not in append mode or are not going to be
                        offset += 1
                        if offset == len(mutDict[dictPos][2]):
                            if (not appending) or (refPos in mutDict):
                                
                                if startedInMiddleOfDictionaryItem:
                                    refPos += 1
                                    startedInMiddleOfDictionaryItem = False
                                appending = False
                                dictPos = refPos
                                offset = 0

            # call the writing function - which will write to the file or pass data to the child nodes as appropriate
            writingFunction(self, node, file, mutDict, insertionDict, **kwargs)

            # de-update the list so it can be used by siblings etc.
            # we need to do everything in the first part but backwards - terrific
            for m in reversed(node.mutations):

                refPos, offset = self.getAbsoluteCoords(m, mutDict, insertionDict)

                if m.mType == mType.SUB:

                    if refPos in mutDict:
                        target = mutDict[refPos][1]
                        mutDict[refPos][1] = target[:offset] + m.source + target[offset+1:]

                        # if the substitution has exactly reversed a previous substitution then we can delete it
                        if mutDict[refPos][1] == self.ref[refPos]:
                            del mutDict[refPos]

                    else:
                        mutDict[refPos] = [self.ref[refPos], m.source, [(m.insertionPos, m.genomePos)]]

                elif m.mType == mType.INS:

                    # de-update an insertion by removing it from the insertion dictionary 
                    # and also removing any reference to it in the mutDict
                    del insertionDict[m.index]
                    insertionLength = len(m.target)

                    # remove the insertion
                    mutDict[refPos][1] = mutDict[refPos][1][:offset + 1] + mutDict[refPos][1][offset + insertionLength + 1:]
                    mutDict[refPos][2] = mutDict[refPos][2][:offset + 1] + mutDict[refPos][2][offset + insertionLength + 1:]

                    # if the de-updated mutation is now identical to the reference genome then delete the dictionary item
                    if mutDict[refPos][0] == mutDict[refPos][1]:
                        del mutDict[refPos]

                # deletion - this is the more difficult case - we need to put back everything that had been deleted. 
                # this may mean putting back an arbitrary string of insertions and substitutions
                else:


                    # need to deal with the deletion one symbol at a time
                    deletedChars = 0
                    counter = 0

                    while (deletedChars < len(m.source)):
                        # a deletion may have skipped characters (other nested deletions)
                        # we need to make sure we only 'un-delete' the correct characters
                        if mutDict[refPos][2][offset] == m.deletedPositions[deletedChars]:
                            
                            # insert back symbols in the target mutation
                            mutDict[refPos][1] = mutDict[refPos][1][:offset] + m.source[deletedChars] + mutDict[refPos][1][offset+1:]
                            deletedChars += 1

                        # move to next symbol
                        if mutDict[refPos][2][offset] == (0, refPos + counter):
                            counter += 1
                        offset += 1
                        if offset == len(mutDict[refPos][1]):
                            offset = 0
                            
                            # remove 'null' deletions
                            if mutDict[refPos][0] == mutDict[refPos][1]:
                                del mutDict[refPos]

                            refPos += counter
                            counter = 0

        return wrappedFunction

    @wrapWriterWithMutationDict
    def writeGenomeShortIndels(self, node, file, mutDict, insertionDict, alternative_output_format):

        # print leaf entry to file
        if node.is_leaf():
            file.write(">" + node.name + "\n")

            if self.verbose:
                print(">" + node.name)
                            
            mutList = list(mutDict.keys())
            mutList.sort()

            for pos in mutList:
                m = mutDict[pos]
                if not alternative_output_format:
                    file.write(f"{m[0]}\t{pos+1}\t{m[1]}\n")

                else:
                    if '-' == m[1][0]: # deletion
                        file.write(f"-\t{pos+1}\t{len(m[1])}\n")
                    
                    else: # insertion and/or substitution
                        if m[0][0] != m[1][0]: # substitution
                            file.write(f"{m[0]}\t{pos+1}\t{m[1]}\n")

                        if len(m[1]) > 1: # insertion
                            file.write(f"i\t{pos+1}\t{m[1][1:]}\n")

                if self.verbose:
                    print(f"{m[0]}\t{pos+1}\t{m[1]}: debug {m}")

        # pass data to children
        else:
            for c in node.children:
                self.writeGenomeShortIndels(c, file, mutDict, insertionDict, alternative_output_format=alternative_output_format)

                  

    def writeGenomeShortNoIndels(self, node, file, mutDict):
        
        for m in node.mutations:
            nuc = m.target
            if nuc != self.ref[m.genomePos]:
                mutDict[m.genomePos] = nuc
            else:
                del mutDict[m.genomePos]

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
                self.writeGenomeShortNoIndels(c, file, mutDict)

        # de-update the dictionary so it can be used by siblings etc.
        for n in range(len(node.mutations)):
            m = node.mutations[len(node.mutations) - (n + 1)]
            nuc = m.source

            if nuc != self.ref[m.genomePos]:
                mutDict[m.genomePos] = nuc
            else:
                del mutDict[m.genomePos]
        

    def write_genome_short(self, tree, output_path, output_file, alternative_output_format):
        # open a file a create a container
        genomefile = open(output_path + output_file + ".txt", "w")
        mutDict = {}
        
        # call the recursive function - there are two depending on whether or not indels are present
        if self.indels:
            insertionDict = {}
            self.writeGenomeShortIndels(
                tree, 
                genomefile, 
                mutDict, 
                insertionDict, 
                alternative_output_format=alternative_output_format)

        else:
            self.writeGenomeShortNoIndels(node=tree, file=genomefile, mutDict=mutDict)
        genomefile.close()


    @wrapWriterWithMutationDict
    def writeGenomeIndels(self, node, file, mutDict, insertionDict):
        
        # print leaf entry to file
        if node.is_leaf():
            # first write the header then the sequence directly from the array

            for k, m in mutDict.items():

                # this deals with insertions and substitutions, and deletions of length 1
                if len(m[0]) == 1:
                    self.refList[k] = m[1].replace("-", "")
                
                else:
                    # deletions of length > 1
                    for i, character in enumerate(m[0]):
                        self.refList[k + i] = ""

            file.write(">" + node.name + "\n"+''.join(self.refList) + "\n")

            # de-update self.refList
            for k, m in mutDict.items():

                # this deals with insertions and substitutions, and deletions of length 1
                if len(m[0]) == 1:
                    self.refList[k] = m[0]
                
                else:
                    # deletions of length > 1
                    for i, character in enumerate(m[0]):
                        self.refList[k + i] = m[0][i]            

        # pass to children
        else:
            for c in node.children:
                self.writeGenomeIndels(c, file, mutDict, insertionDict)

    def writeGenomeNoIndels(self, node, file, nRefList):
        # function to write a complete sequence output iteratively
        # update list
        for m in node.mutations:
            nRefList[m.genomePos] = m.target
        # print leaf entry to file
        if node.is_leaf():
            # first write the header then the sequence directly from the array
            file.write(">" + node.name + "\n"+''.join(nRefList) + "\n")

        # pass to children
        else:
            for c in node.children:
                self.writeGenomeNoIndels(c, file, nRefList)
        # de-update the list so it can be used by siblings etc.
        for n in range(len(node.mutations)):
            m = node.mutations[len(node.mutations) - (n + 1)]
            nRefList[m.genomePos] = m.source

    def write_genome(self, tree, output_path, output_file, refList):
        file = open(output_path + output_file + ".fasta", "w")
        # convert reference list to an array
        #refList = np.array(refList)
        if self.indels:
            mutDict, insertionDict = {}, {}
            self.writeGenomeIndels(node=tree, file=file, mutDict=mutDict, insertionDict=insertionDict)

        else:
            self.writeGenomeNoIndels(tree, file, refList)
        file.close()



    def writeGenomePhylip(self, node, file, nRefList):
        # function to write a complete sequence output iteratively
        # update list
        for m in node.mutations:
            nRefList[m.genomePos] = m.target
        # print leaf entry to file
        if node.is_leaf():
            file.write(node.name + "\t" + (''.join(nRefList)) + "\n")
        else:
            for c in node.children:
                self.writeGenomePhylip(c, file, nRefList)
        # de-update the list so it can be used by siblings etc.
        for n in range(len(node.mutations)):
            m = node.mutations[len(node.mutations) - (n + 1)]
            nRefList[m.genomePos] = m.source


    def write_genome_phylip(self, tree, output_path, output_file, refList):

        if self.indels:
            print("Phylip output format is not supported yet with indels.")
        
        else:
            # open a file for the phylip output
            file = open(output_path + output_file + ".phy", "w")
            file.write("\t" + str(len(tree)) + "\t" + str(len(self.ref)) + "\n")
            # run the recursive function to write the phylip formatted file
            self.writeGenomePhylip(node=tree, file=file, nRefList=refList)
            file.close()


    def write_genome_mat(self, tree, output_path, output_file):

        mat = protobuf.data()
        mat.newick = tree.write(format=1)

        self.writeGenomeMAT(tree, mat)

        f = open(output_path + output_file + ".mat.pb", "wb")
        f.write(mat.SerializeToString())
        f.close()

    def writeGenomeMAT(self, node, mat):

        _ = mat.metadata.add()
        mutations_protobuf = mat.node_mutations.add()
        for m in node.mutations:
            m_pb = mutations_protobuf.mutation.add()
            m_pb.position = m.genomePos + 1

            if getattr(m, "insertionPos", 0):
                m_pb.insertion_position = m.insertionPos
                m_pb.ref_nuc = -1
            else:
                m_pb.ref_nuc = self.alleles[self.ref[m.genomePos]]

            if m.mType == mType.INS:
                m_pb.insertion_index = m.index

            ch = m.source[0]
            if ch == "-":
                m_pb.par_nuc = -1
            else:
                m_pb.par_nuc = self.alleles[ch]

            for ch in m.target:
                if ch == "-":
                    m_pb.mut_nuc.append(-1)
                else:
                    m_pb.mut_nuc.append(self.alleles[ch])

        for c in node.children:
            self.writeGenomeMAT(c, mat)




class GenomeTree_vanilla:
    def __init__(self, nCat, ref, mutMatrix, categories, categoryRates, hyperMutRates, hyperCategories, infoFile, verbose):
        self.nCat = nCat
        self.ref = ref
        self.mutMatrix = mutMatrix
        self.categories = categories
        self.categoryRates = categoryRates
        self.hyperMutRates = hyperMutRates
        self.hyperCategories = hyperCategories

        self.infoFile = infoFile
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
            if self.infoFile!=None:
                self.infoFile.write(str(pos + 1) + "\t" + str(cat) + "\t" + str(hyp) + "\t")
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
                if self.infoFile!=None:
                    self.infoFile.write(self.allelesList[i] + "\t" + self.allelesList[j] + "\n")
            else:
                positions[cat][a].append([pos, 0])
                if self.infoFile!=None:
                    self.infoFile.write(".\t" + ".\n")
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
        if self.infoFile!=None:
            self.infoFile.close()

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
        for i in self.range4:
            for j in self.range4:
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
            for i in self.range4:
                for k in range(len(parentMuts[c][i])):
                    childMutations[c][i].append(list(parentMuts[c][i][k]))
        extrasChild = []
        for e in range(len(extrasParent)):
            extrasChild.append(list(extrasParent[e]))

        # Sample new mutation event with Gillespie algorithm
        currTime += np.random.exponential(scale=1.0 / rate)
        if self.verbose:
            print(f"\n Node {childNode.name} BLen: {bLen}. Mutations at start:")
            print([str(m) for m in parentMuts])
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
                for i in self.range4:
                    for j in self.range4:
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
                print([str(m) for m in childMutations])

        if self.verbose:
            print("mutations at the end:")
            print([str(m) for m in childMutations])
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



    def write_genome_short(self, tree, output_path, output_file, **kwargs):
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

    def __str__(self):
        if self.isTerminal:
            props = self.__dict__
            short_props = {k:v for k, v in props.items() if k in ["rate", "level", "isTerminal", "allele"]}
            return "<" +str(short_props) +">"
        else:
            f = self.belowNodes[0]
            s = self.belowNodes[1]
            return f"({f},rate:{self.rate},level:{self.level},{s})"


class mutation:
    
    def __init__(self, mType, genomePos, **kwargs):
        # mutations should have the following properties (depending on their mType):
        # genomePos - 0 indexed position
        # insertionPos (if indels) - an index specifying which insert the mutation belongs to (0 for the reference)
        # source (i.e. what the sequence was originally)
        # target (i.e. what the sequence was after the mutation)
        # length (if DEL)
        # N.B. I am using 'source' and 'target' rather than other words because 'to' and 'from' are python keywords,
        # so are 'input' and 'output'.
        # index  (if INS) insertions have an index so that other mutations can refer to them
        self.mType = mType
        self.genomePos = genomePos
        for k, v in kwargs.items():
            setattr(self, k, v)



    def __str__(self):

        string = f"{self.source}{self.genomePos + 1}{self.target}" 

        if getattr(self, "insertionPos", 0):
            string += ("@" + str(self.insertionPos))
        
        if self.mType == mType.INS:
            string = str(self.index) + ":" + string

        return string


class mType(Enum):
    SUB = 1
    INS = 2
    DEL = 3



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
        indices = codonIndices[i1]
        aa = translationList[i1]
        indices2 = list(indices)
        for i2 in range(3):
            for i3 in range(3):
                indices2[i2] = ((indices[i2] + i3 + 1) % 4)
                cod2 = codonIndices2[indices2[0], indices2[1], indices2[2]]
                if aa != translationList[cod2]:
                    isNonsynom[i1, i2, i3] = True
                if translationList[cod2] == "*" and aa != "*":
                    isIntoStop[i1, i2, i3] = True
                # print(codonAllelesList[i1]+" "+codonAllelesList[cod2]+" into stop.")
            indices2[i2] = indices[i2]

    return isNonsynom, isIntoStop



def writeGenomeNewick(node):
    # function to write a newick output iteratively
    if node.is_leaf():
        outString = node.name + '[&mutations={'
    else:
        outString = '('
        for c in range(len(node.children)):
            outString += writeGenomeNewick(node.children[c])
            if c < len(node.children) - 1:
                outString += ','
        outString += ')[&mutations={'
    stringToAdd = ''
    for i in range(len(node.mutAnnotation)):
        stringToAdd += node.mutAnnotation[i]
        if i < len(node.mutAnnotation) - 1:
            stringToAdd += ','
    stringToAdd += ('}]:' + str(node.dist))
    return outString + stringToAdd
