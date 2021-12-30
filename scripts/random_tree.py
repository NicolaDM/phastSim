# encoding: utf-8

"""
Random Phylogenetic Tree Generator.
This has been built on top of 
https://raw.githubusercontent.com/tresoldi/ngesh/master/ngesh/random_tree.py
but has been modified to make it more efficient for large tree.

This script provides function to generate random phylogenetic trees in
a Yule (birth only) or Birth-Death model, setting different generation
parameters and limiting the tree in terms of number of leaves and/or
evolution time.
"""

# Import Python standard libraries
import math
import random

# Import 3rd party libraries
import numpy as np
from ete3 import Tree

# Import other modules from this library
#import utils.py
import hashlib

# Define the maximum number of tries for generation
__MAX_ATTEMPTS = 3000


#def __extant(tree):
#     """
#     Internal function returning a list of non-extinct leaves in a tree.
# 
#     Parameters
#     ----------
# 
#     tree: ete3 tree object
#         The tree whose nodes will be checked.
# 
#     Returns
#     -------
#     leaves: list
#         List of extant leaves.
#     """

    # Return a filtered list compiled with a list comprehension; the
    # 'extinct' field is not part of ETE3 defaults, but we use here in
    # order to easily differentiate between alive and extinct leaves in
    # Birth-Death models.
    #return [leaf for leaf in tree.get_leaves() if leaf.extinct is False]

def set_seeds(seed):
    """
    Set seeds globally from the user provided one.

    The function takes care of reproducibility and allows to use strings and
    floats as seed for `numpy` as well.
    """

    random.seed(seed)

    # allows using strings as np seeds, which only takes uint32 or arrays of
    # NOTE: this won't set the seed if it is None: if you want to seed none
    #       as seed, manually call np.random.seed()
    if isinstance(seed, (str, float)):
        np_seed = np.frombuffer(
            hashlib.sha256(str(seed).encode("utf-8")).digest(), dtype=np.uint32
        )
    else:
        np_seed = seed

    # Set the np set
    np.random.seed(np_seed)

def label_tree(leaves, model, seed=None):
    """
    Labels the nodes of a tree according to a model.

    Linguistic labels are unique names generated in a way intended to be
    readable.

    Please note that the `tree` object is changed in place (no return).

    Parameters
    ----------

    tree: ete3 tree object
        The tree whose nodes will be labeled.
    model : str
        A string indicating which model for label generation should be
        used. Possible values are "enum" (for enumerated labels), "human"
        (for random single names), and "bio" (for random biological names).
    seed : value
        An optional seed for the random number generator, only used in case
        of linguistic and biological labels. Defaults to `None`.
    """

    # Cache the leaves, so we can also obtain their number
    #leaves = tree.get_leaves()
# 
#     if model == "bio":
#         # As we are using a simple model with replacements, even if
#         # extremely unlikely, we might have repeated items in the labels.
#         # The execution would not fail as we are using `zip()`, only items
#         # would be unnamed, but we are manually adding missing labels as
#         # enumerations to make sure there are no anynomous nodes.
#         # TODO: decide on better approach or make case explicit in docs
#         species = sorted(set(utils.random_species(len(leaves), seed)))
#         species += ["L%i" % i for i in range(len(leaves) - len(species))]
# 
#         for leaf_node, name in zip(leaves, species):
#             leaf_node.name = name
# 
#     elif model == "human":
#         for leaf, name in zip(leaves, utils.random_labels(len(leaves), seed)):
#             leaf.name = name
# 
#     else:
        # Build the pattern for the label, including the computation of the
        # number of padding zeros needed
    pattern = "L%%0%ii" % (1 + math.floor(math.log10(len(leaves))))

    # Label all leaves first
    for leaf_idx, leaf_node in enumerate(leaves):
        leaf_node.name = pattern % (leaf_idx + 1)
        
def delete_single_child_internal(t):
    """Utility function that removes internal nodes
    with a single child from tree"""

    for node in t.traverse("postorder"):
        if(not node.is_leaf() and len(node.get_children()) < 2):
            node.delete()

    while len(t.get_children()) == 1:
        t = t.children[0]
        t.up = None


def __gen_tree(**kwargs):
    """
    Internal function for tree generation.

    This is an internal function for the tree generation, whose main
    difference to `gen_tree()`, the one exposed to the user, is that it
    does not guarantee that a tree will be generated, as the parameters and
    the random sampling might lead to dead-ends where all the leaves in
    a tree are extinct before any or all the stopping criteria are met.

    As an internal function, it does not set default values to the arguments
    and does not perform any checking on the values. Information on the
    arguments, which have the same variable names and properties, are given
    in the documentation for `gen_tree()`.
    """

    # Initialize the RNG
    set_seeds(kwargs["seed"])

    # Compute the overall event rate (birth plus death), from which the
    # random expovariate will be drawn. `birth` is here normalized in range
    # [0..1] so that we can directly compare with the results of
    # `.random()` and decide if the event is a birth or a death.
    # `death` does not need to be normalized, as it is not used anymore (the
    # only check, below, is `.random() <= birth`).
    event_rate = kwargs["birth"] + kwargs["death"]
    birth = kwargs["birth"] / event_rate

    # Create the tree root as a node. Given that the root is at first set as
    # non-extinct and with a branch length of 0.0, it will be immediately
    # subject to either a speciation or extinction event.
    tree = Tree()
    leaves=[tree]
    tree.dist = 0.0
    tree.extinct = False
    
    #list of currently existing leaves.
    #this saves us recurring along the phylogeny each time a sample is collected.
    leaves=[tree]

    # Iterate until an acceptable tree is generated (breaking the loop with
    # a tree) or all leaves go extinct (breaking the loop with `tree` as None).
    # `total_time`, of which we keep track in case `max_time` is provided,
    # is the total evolution time (sum of branch lengths) from the root to the
    # extant nodes.
    total_time = 0.0
    while True:
        # Get the list of extant leaves
        #leaf_nodes = __extant(tree)

        # Compute the event time before the next birth/death event from a
        # random exporaviate reflecting the number of extant leaves and the
        # combined event probability.
        event_time = random.expovariate(len(leaves) * event_rate)

        # Update the total evolution time. If a maximum alloted time
        # `max_time` is provided and we overshoot it, break the loop
        # without implementing the event (as, by the random event time, it
        # would take place *after* our maximum time, in the future).
        total_time += event_time
        if kwargs["max_time"] and total_time > kwargs["max_time"]:
            for node in leaves:
                node.dist=kwargs["max_time"]-node.dist
            break

        # Select a random node among the extant ones and set it as extinct
        # before simulating either a birth or death event; the type of
        # event is decided based on the comparison of the result of a
        # `random.random()` call with `birth` (here already normalized in
        # relation to `event_rate`)
        #node = np.random.choice(leaf_nodes)
        leafN=np.random.random_integers(len(leaves))-1
        node=leaves[leafN]
        node.extinct = True
        node.dist=total_time - node.dist
        if np.random.random() <= birth:
            # The event will be a birth (i.e., speciation one), with at least
            # two children (the number is increased by a random sample from a
            # Poisson distribution using the `lam` parameter, so that
            # hard politomies are possible). The distance
            # of the children is here initially set to the time of its creation, and will be
            # set to its actual length at the time it's removed or simulation ends.
            child_node = Tree()
            child_node.dist = total_time
            child_node.extinct = False
            node.add_child(child_node)
            leaves[leafN]=child_node
            poiSample=np.random.poisson(kwargs["lam"])
            for _ in range(1 + poiSample):
                child_node = Tree()
                child_node.dist = total_time
                child_node.extinct = False
                node.add_child(child_node)
                leaves.append(child_node)
        else:
            del leaves[leafN]
            

        # (Re)Extract the list of extant nodes, now that we might have new
        # children and that the randomly selected node went extinct
        # (easier than directly manipulating the Python list). From the
        # updated list, we will extend the branch length of all extant leaves
        # (thus including any new children) by the `event_time` computed
        # above.
        #Nicola: I made this more efficient now, plus, there was a mistake here
        # as the event time was added also to the newly created leaves, which instead should have kept
        #dist=0.
        #leaf_nodes = __extant(tree)
        #for leaf in leaf_nodes:
        #    new_leaf_dist = leaf.dist + event_time
        #    leaf.dist = min(
        #        new_leaf_dist, (kwargs["max_time"] or new_leaf_dist)
        #    )

        # If the event above was a death event, we might be in the undesirable
        # situation where all lineages went extinct before we
        # could finish the random generation according to the
        # user-requested parameters, so that one or both stopping criteria
        # cannot be satisfied. A solution could
        # be to recursively call this function, with the same
        # parameters, until a valid tree is found, but this is not
        # optimal (nor elegant) and might get us stuck in a
        # loop if we don't keep track of the number of iterations
        # (especially if we got to this point by using a
        # user-provided random seed and/or set of unfortunate parameters).
        # In face of that, it is preferable to be explicit about the problem by
        # returning a `None` value, with the user (or a wrapper
        # function) being in charge of asserting that the desired
        # number of random trees is collected (even if it is a single one).
        if len(leaves)==0:
            tree = None
            break

        # Check whether the number of leaves stopping criterium was reached
        if kwargs["min_leaves"] and len(leaves) >= kwargs["min_leaves"]:
            for node in leaves:
                node.dist=total_time-node.dist
            break

        #if kwargs["max_time"] and total_time >= kwargs["max_time"]:
        #    break

    # In some cases we might end up with technically valid trees composed
    # only of the root. We make sure at least one speciation event took
    # place, returning `None` as failure in other cases.
    #Nicola: I am removing this bit as it doesn't do what it's supposed to do, and as it's not
    #clear why one would do what this is supposed to do.
    #if tree and len(__extant(tree)) <= 2:
    #    tree = None

    # Prune the tree, removing extinct leaves, if requested and if a
    # tree was found. Remember that the ete3 `prune()` method takes a list
    # of the nodes that will be kept, removing the other ones.
    if kwargs["prune"] and tree:
        tree.prune(leaves)

    # Label the tree before returning it, if it was provided
    if kwargs["labels"] and tree:
        label_tree(leaves, kwargs["labels"], seed=kwargs["seed"])

    return tree


def gen_tree(birth, death, **kwargs):
    """
    Return a random phylogenetic tree.

    At least one stopping criterion must be informed, with the tree being
    returned when such criterion, or either criteria, is/are met.

    This function wraps the internal `__gen_tree()` function which cannot
    guarantee that a valid tree will be generated given the user
    parameters and the random sampling. It will try as many times as
    necessary to provide a valid (and reproducible, given a `seed`) tree,
    within the limits of an internal parameter for maximum number of
    attempts.

    Parameters
    ----------

    birth : float
        The birth rate (lambda) for the generated tree.
    death : float
        The death rate (mu) for the generated tree. Must be explicitly set
        to zero for Yule model (i.e., birth only).
    min_leaves : int
        A stopping criterion with the minimum number of extant leaves.
        The generated tree will have at least the number of requested
        extant leaves (possibly more, as the last speciation event might
        produce more leaves than the minimum specified.
        Defaults to None.
    max_time : float
        A stopping criterion with the maximum allowed time for evolution.
        Defaults to None.
    labels : str or None
        The model to be used for generating random labels, either
        "enum" (for enumerated labels), "human" (for random single names),
        "bio" (for random biological names" or None. Defaults to "enum".
    seed : value
        An optional seed for the random number generator. Defaults to None.
    lam : float
        The expectation of interval for sampling a Poisson distribution
        during speciation, with a minimum of two descendants. Should be used
        if more than two descendants are to be allowed. Defaults to zero,
        meaning that all speciation events will have two and only two
        descendents.
    prune : bool
        A flag indicating whether any non-extant leaves should be pruned from
        the tree before it is returned.

    Returns
    -------

    tree : ete3 tree
        The tree randomly generated according to the parameters.
    """

    # Parse kwargs and set defaults
    min_leaves = kwargs.get("min_leaves", None)
    max_time = kwargs.get("max_time", None)
    labels = kwargs.get("labels", "enum")
    lam = kwargs.get("lam", 0.0)
    prune = kwargs.get("prune", False)
    seed = kwargs.get("seed", None)

    # Confirm that at least one stopping condition was provided
    if not (min_leaves or max_time):
        raise ValueError("At least one stopping criterion is required.")

    # Confirm that a valid `labels` was passed
    if labels not in ["enum", "human", "bio", None]:
        raise ValueError("Invalid label model provided ('%s')" % labels)

    # Generate the random tree
    cur_attempt = 0
    while True:
        # Update the current attempt and drop out if maximum was reached.
        # In the extremely unlikely event we were not able to generate a tree,
        # after all these iterations (probably due to the user parameters),
        # we fail explicitly.
        cur_attempt += 1
        if cur_attempt == __MAX_ATTEMPTS:
            raise RuntimeError("Unable to generate a valid tree.")

        # Seed the RNG with the current seed (which, in the first run, will
        # either be the one provided by the user or `None`) and extract a
        # new seed for future iterations if needed (if the provided seed fails,
        # in case of no tree generation, there is no point in trying the
        # feed again).
        set_seeds(seed)
        seed = np.random.random()

        # Ask for a new tree
        tree = __gen_tree(
            birth=birth,
            death=death,
            min_leaves=min_leaves,
            max_time=max_time,
            labels=labels,
            lam=lam,
            prune=prune,
            seed=seed,
        )

        # Break out of the loop if a valid tree was found, as in most of the
        # cases; if no tree could be generated, `__gen_tree()` will return
        # `None`.
        if tree:
            break
    delete_single_child_internal(tree)
    return tree

