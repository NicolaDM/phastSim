import pytest
from pytest import approx
from types import SimpleNamespace
import numpy as np
from ete3 import Tree, TreeNode
from phastSim import phastSim




def test_phastSim_imported_correctly():
    
    # if this test failed something has gone very wrong and most likely it is due to an import failing. 
    c = phastSim.Constants()
    print(c.alleles)
    assert c.alleles == {'A':0, 'C':1, 'G':2, 'T':3}

def create_standard_mock_arguments():
    return SimpleNamespace(
        outpath="_",
        reference="_",
        rootGenomeLength=10,
        rootGenomeFrequencies=[],
        treeFile="_",
        scale=1,
        seed=1,
        alpha=0,
        invariable=0.0,
        mutationRates=[0.039, 0.310, 0.123, 0.140, 0.022, 3.028, 0.747, 0.113, 2.953, 0.056, 0.261, 0.036],
        codon=False,
        omegaAlpha=0.0,
        omegaCategoryProbs=[1.0],
        omegaCategoryRates=[1.0],
        categoryProbs=[1.0],
        categoryRates=[1.0],
        hyperMutProbs=[],
        hyperMutRates=[],
        indels=True,
        insertionRate=0.1,
        deletionRate=0.2,
        insertionLength=0.1, 
        deletionLength=0.9, 
        noHierarchy=False,
        verbose=True,
        outputFile="_",
        createNewick="_",
        createFasta="_",
        createPhylip="_",
        createInfo="_"
    )

def test_phastSimRun_init_rootGenome():
    
    mock_args = create_standard_mock_arguments()

    sim_run = phastSim.phastSimRun(mock_args)

    assert len(sim_run.create_rootGenome_nuc()[1]) == 10
    
    sim_run.args.codon = True
    assert len(sim_run.create_rootGenome_codon()[1]) == 9

    assert sim_run.init_substitution_rates()[0][0] == -0.472

def test_phastSimRun_init_gamma_rates():
    np.random.seed(10)
    mock_args = create_standard_mock_arguments()
    mock_args.alpha = 0.1
    
    sim_run = phastSim.phastSimRun(mock_args)

    ref, _  = sim_run.init_rootGenome()

    assert len(ref) == 10

    gammaRates = sim_run.init_gamma_rates()

    assert next(gammaRates) == approx(0.2286, 0.01)



def setup_genome_tree():
    np.random.seed(10)
    mock_args = create_standard_mock_arguments()
    mock_args.nCodons=3
    sim_run = phastSim.phastSimRun(mock_args)
    ref, _  = sim_run.init_rootGenome()
    mutMatrix = sim_run.init_substitution_rates()
    gammaRates = sim_run.init_gamma_rates()

    # set up hypermutation rates
    hyperCategories = sim_run.init_hypermutation_rates()

    # set up codon substitution model
    omegas = sim_run.init_codon_substitution_model()
    ins, ins_len, dels, dels_len = sim_run.init_indel_rates()
    insertion_freqs = sim_run.init_insertion_frequencies()

    t = Tree("(A:1,(B:1,(E:1,D:1):0.5):0.5);")
    print(t)

    genome_tree = phastSim.GenomeTree_hierarchical(
        nCodons=mock_args.nCodons,
        codon=mock_args.codon,
        ref=ref,
        gammaRates=gammaRates,
        omegas=omegas,
        mutMatrix=mutMatrix,
        hyperCategories=hyperCategories,
        hyperMutRates=sim_run.args.hyperMutRates,
        indels=True,
        insertionRate=ins,
        insertionLength=ins_len,
        insertionFrequencies=insertion_freqs,
        deletionRate=dels,
        deletionLength=dels_len,
        scale=sim_run.args.scale,
        file=None,
        verbose=True
    )
    return genome_tree

def test_GenomeTree_hierarchical_initialises_without_indels():

    np.random.seed(10)
    mock_args = create_standard_mock_arguments()
    mock_args.nCodons=4
    mock_args.rootGenomeLength=12
    mock_args.codon = True

    sim_run = phastSim.phastSimRun(mock_args)
    ref, _  = sim_run.init_rootGenome()
    mutMatrix = sim_run.init_substitution_rates()
    gammaRates = sim_run.init_gamma_rates()

    # set up hypermutation rates
    hyperCategories = sim_run.init_hypermutation_rates()

    # set up codon substitution model
    omegas = sim_run.init_codon_substitution_model()


    genome_tree = phastSim.GenomeTree_hierarchical(
        nCodons=mock_args.nCodons,
        codon=mock_args.codon,
        ref=ref,
        gammaRates=gammaRates,
        omegas=omegas,
        mutMatrix=mutMatrix,
        hyperCategories=hyperCategories,
        hyperMutRates=sim_run.args.hyperMutRates,
        indels=False,
        insertionRate=0,
        insertionLength=0,
        insertionFrequencies=[],
        deletionRate=0,
        deletionLength=0,
        scale=sim_run.args.scale,
        file=None,
        verbose=True
    )    

    genome_tree.populate_genome_tree()
    genome_tree.check_start_stop_codons()
    genome_tree.normalize_rates()
    print(genome_tree)

def test_GenomeTree_hierarchical_initialises():

    genome_tree = setup_genome_tree()

    genome_tree.populate_genome_tree()
    print(genome_tree)
    genome_tree.check_start_stop_codons()
    print(genome_tree)
    genome_tree.normalize_rates()
    print(genome_tree)

def test_genomeTree_hierarchical_findPos_without_indels():
    pass


def test_genomeTree_hierarchical_findPos_single_deletion():

    genome_tree = setup_genome_tree()

    genome_tree.populate_genome_tree()
    genome_tree.normalize_rates()

    print(genome_tree)
    # these initial conditions should produce a deletion
    newGenomeNode = phastSim.genomeNode(level=1)
    newGenomeNode.belowNodes = list(genome_tree.genomeRoot.belowNodes)
    mutation = genome_tree.findPos(0.2566, newGenomeNode, level=1)
    print(newGenomeNode)
    print(mutation.length, mutation.insertionPos, mutation.genomePos, mutation.mType)
    genome_tree.deleteNodes(0.2566, newGenomeNode, mutation.length, level=1)
    print(newGenomeNode)

def test_genomeTree_hierarchical_findPos_single_insertion():

    genome_tree = setup_genome_tree()

    genome_tree.populate_genome_tree()
    genome_tree.check_start_stop_codons()
    genome_tree.normalize_rates()

    print(genome_tree)
    print(genome_tree.alleles)
    n = genome_tree.genomeRoot.belowNodes[0].belowNodes[1].belowNodes[0].belowNodes[0] 
    # assert that this is a terminal node otherwise the test doesn't really make sense
    assert n.isTerminal

    mutation = genome_tree.sampleInsertion(rate=0.8, node=n, rand=0.3, level=1)

    print(mutation)
    
    print(n)

def test_genomeTree_hierarchical_findPos_single_deletion_codon_model():
    pass

def test_genomeTree_hierarchical_findPos_single_insertion_codon_model():
    pass


def test_printing_functions():
    
    genome_tree = setup_genome_tree()

    class MockFile:
        """
        A class which has 1 method called write, 
        so that I can use it as a mock file object in this test.
        """
        def __init__(self):
            self.written_data = ""

        def write(self, string):
            self.written_data += string


    def parameterised_test(mutDict, insertionDict, mutations, expected_output):
        f = MockFile()
        node = TreeNode(name="test_node")
        node.mutations = mutations 

        original_mD = mutDict.copy()
        original_iD = insertionDict.copy()
        genome_tree.writeGenomeShortIndels(node=node, file=f, mutDict=mutDict, insertionDict=insertionDict)
        # the whole point of this function is that the genome tree updates and then de-updates
        # any mutations. So we need the mutDict and insertionDict to remain the same before and after printing. 
        assert mutDict == original_mD
        assert insertionDict == original_iD
        assert f.written_data == expected_output #

    parameterised_test(
        mutDict={}, 
        insertionDict={}, 
        mutations=[phastSim.mutation(mType=phastSim.mType.DEL, genomePos=4, insertionPos=0, source="GAAT", target="")], 
        expected_output=">test_node\nGAAT5----\n"
    )