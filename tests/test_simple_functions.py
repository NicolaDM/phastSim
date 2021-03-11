import pytest
from types import SimpleNamespace
import numpy as np
from ete3 import Tree
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
    assert len(sim_run.create_rootGenome_codon()[1]) == 9

    assert sim_run.init_substitution_rates()[0][0] == -0.472

def test_phastSimRun_init_gamma_rates():

    mock_args = create_standard_mock_arguments()
    mock_args.alpha = 0.1
    
    sim_run = phastSim.phastSimRun(mock_args)

    ref, _  = sim_run.init_rootGenome()

    assert len(ref) == 10

    gammaRates = sim_run.init_gamma_rates()

    print(gammaRates)


def test_GenomeTree_hierarchical_initialises_without_indels():

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
    gammaRates, omegas = phastSim.check_start_stop_codons(ref=ref, gammaRates=gammaRates, omegas=omegas)


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
        deletionRate=0,
        file=None,
        verbose=True
    )

    genome_tree.populateGenomeTree(node=genome_tree.genomeRoot)
    genome_tree.normalize_rates(scale=1)

def test_GenomeTree_hierarchical_initialises():

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
    gammaRates, omegas = phastSim.check_start_stop_codons(ref=ref, gammaRates=gammaRates, omegas=omegas)
    ins, dels = sim_run.init_indel_rates()

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
        deletionRate=dels,
        file=None,
        verbose=True
    )

    genome_tree.populateGenomeTree(node=genome_tree.genomeRoot)
    #print(genome_tree)
    genome_tree.normalize_rates(scale=1)
    #print(genome_tree)

def test_genomeTree_hierarchical_findPos_without_indels():
    pass


def test_genomeTree_hierarchical_findPos():
    pass