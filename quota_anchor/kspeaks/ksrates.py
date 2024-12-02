# this module slightly modified from ksrates(https://github.com/VIB-PSB/ksrates)

import logging
import sys
from ete3 import Tree

logger = logging.getLogger('main.ksrates')

def get_newick_tree(tree_string):
    """
    Gets the config file field of the Newick tree.
    Checks and exits if the species' names in the Newick tree contain illegal characters (underscore or spaces).

    :return tree_string: the tree object by ete3
    """

    if not (tree_string.endswith(';')):
        tree_string += ";"
    if tree_string == "();" or tree_string == ";":
        logger.error('Parameter "newick_tree" is empty, please fill in')
        sys.exit(1)
    try:
        tree = Tree(tree_string)
    except Exception:
        logger.error('Unrecognized format for parameter "newick_tree" (for example, parentheses do not match)')
        sys.exit(1)

    # Check if species' informal names contain illegal characters (underscore or spaces)
    species_illegal_char=[]
    for informal_name in tree.get_leaf_names():
        if "_" in informal_name or " " in informal_name:
            species_illegal_char.append(informal_name)
    if len(species_illegal_char) != 0:
        logger.error(f"Informal species' names must not contain any spaces or underscores. Please change the following names in commandline or configuration file:")
        for informal_name in species_illegal_char:
            logger.error(f"- {informal_name}")
        sys.exit(1)
    return tree

def check_integrity_newick_tree(tree):
    """
    :param tree: the original tree object

    Checks if there are syntax errors in the newick_tree input. Exists if there are errors.

    - Case 1: The presence of extra unnecessary pairs of parenthesis generates internal nodes with only one child node,
      instead of two children nodes; this will rise problems during the parsing of the tree to obtain the species trios.
      Therefore, the code exists and prompts the user to remove such unnecessary parentheses.

      Example: the input Newick tree contains a subtree whose outermost pair of parenthesis has to be removed.
      String visualization of the subtree:   (((elaeis,oryza),asparagus))
      ASCII visualization of the subtree - note the extra node at the base of this subtree:
                                     /-elaeis
                                  /-|
                            -- /-|   \-oryza
                                 |
                                  \-asparagus

    - Case 2: In presence of unresolved phylogeny (i.e. three or more children nodes branching off from an internal node)
      there will be problems in downstream analysis due to ambiguous outgroup relationships.
      Therefore, the code exists and prompts the user to rearrange the node(s).
      
      Example: the input Newick tree contains a subtree where the basal node has three children nodes.
      String visualization of the subtree: (elaeis,oryza,maize)
      ASCII visualization of the subtree:   
                               /-elaeis
                              |
                            --|--oryza
                              |
                               \-maize
    """
    # For each internal node, check integrity (must have exactly two children)
    # Checking structural integrity of input Newick tree...
    trigger_exit = False
    internal_nodes_with_one_child, internal_nodes_with_three_children = [], []
    for node in tree.traverse():
        if not node.is_leaf():
            number_of_children_nodes = len(node.get_children())
            if number_of_children_nodes == 1:
                internal_nodes_with_one_child.append(node)
            elif number_of_children_nodes > 2:
                internal_nodes_with_three_children.append(node)

    if len(internal_nodes_with_one_child) != 0:
        logger.error(f'The tree structure has one ore more incomplete internal nodes:')
        logger.error(f"likely there are unnecessary pairs of parentheses that generate internal nodes with only one child node instead of two children nodes")
        logger.error(f"Please adjust the input tree as suggested below and rerun the analysis")
        logger.error(f"Such syntax error can be solved by removing the unnecessary outermost pair of parentheses in the following subtree(s):\n")
        for node in internal_nodes_with_one_child:
            logger.error(f'Subtree {internal_nodes_with_one_child.index(node)+1}: {node.write(format=9).rstrip(";")}{node}\n')
        trigger_exit = True

    if len(internal_nodes_with_three_children) != 0:
        logger.error(f'The tree structure contains unresolved phylogenetic relationships')
        logger.error(f"Please adjust the tree so that each internal node has exactly two children nodes")
        logger.error(f"Such structural issue has been encountered at the base of the following subtree(s):\n")
        for node in internal_nodes_with_three_children:
            logger.error(f'Subtree {internal_nodes_with_three_children.index(node)+1}: {node.write(format=9).rstrip(";")}{node}\n')
        trigger_exit = True

    if trigger_exit:
        sys.exit(1)

def get_species_node(species_name, tree):
    """
    :param species_name: the name of the focal species as it appears in the Newick tree
    :param tree: the Newick tree 
    :return species_node: the tree node object associated to the focal species
    """
    species_node = tree.search_nodes(name=species_name)
    if len(species_node) == 0:
        logger.error(f"The focal species {species_name} was not found in the provided tree.")
        sys.exit(1)
    return species_node

def reorder_tree_leaves(tree, focal_species):
    """
    :param tree: the original tree object
    :param focal_species: the name of the focal species
    :return: a new equivalent tree object with the focal species as top leaf
    """
    species_node = get_species_node(focal_species, tree)[0]
    sister_node = species_node.get_sisters()[0]
    sister_node = sister_node.write(format=9).rstrip(";")
    string = f"({species_node.name},{sister_node})"

    node = species_node.up
    while node != tree.get_tree_root():
        sister_node = node.get_sisters()[0]
        sister_node = sister_node.write(format=9).rstrip(";")
        string = "(" + string + "," + sister_node + ")"
        node = node.up
    
    string = string + ";"
    ordered_tree = Tree(string)
    return ordered_tree

def get_species_history(species_node):
    """
    :param species_node: the tree node object associated to the focal species
    :return species_history: list of the ancestor nodes of the focal species, from the focal species
                             included until the root included
    """
    species_history = []
    ancestors = species_node[0].get_ancestors()   # ancestors are in order from leaf to root
    species_history.extend(species_node)
    species_history.extend(ancestors)
    return species_history

def labeling_internal_nodes(species_node):
    """
    Labels the ancestor nodes of the focal species with numbers, starting from 1, until the root.
    Also labels the focal species and its ancestor nodes with a specific type feature.
    
    :param species_node: the tree node object associated to the focal species
    """
    species_node[0].add_feature("type", "species_of_interest")
    node_label = 1
    for ancestor in species_node[0].iter_ancestors():
        ancestor.name = node_label  # the name label to be shown in the ASCII tree will start from 1
        ancestor.add_feature("type", "ancestor")
        node_label = node_label + 1

def get_sister_species_of_a_node(currentnode):
    """
    :param currentnode: the current node 
    :return: the list containing the names of the sister species of the current node.
    """
    # First get the sister node of currentnode
    sis_node = currentnode.get_sisters()     # list with sister NODE(s) (not with the sister LEAVES...) -> there is only one sister node if the tree is binary
    if len(sis_node) != 1:
        sys.exit(print("One or more phylogenetic relationships in the tree are unresolved (>2 branches). Exiting."))
    # Then get the leaves of the sister node (multiple sister leaves are instead allowed of course!)
    sisters = []    # list of all sister species
    sisters_leaves = sis_node[0].get_leaves()
    for sisters_leaf in sisters_leaves:
        sisters.append(sisters_leaf.name)
    return sisters


def get_outspecies_of_a_node(currentnode, max_num_outspecies):
    """
    :param currentnode: the current node
    :param max_num_outspecies: the maximum number (N) of outspecies allowed for the adjustment (only the N closest will be considered)
    :return: a list containing the names of the N outgroup(s) of the current node.
    """

    outspecies = []     # list of all outspecies of the current parent node
    for ancestor in currentnode.iter_ancestors():
        outgroup_node = ancestor.get_sisters()    # get sister NODE(s) of the current parent node
        for branch in outgroup_node:
            outspecies_leaves = branch.get_leaves() # get the outspecies leaves
            for outspecies_leaf in outspecies_leaves:
                if len(outspecies) < max_num_outspecies:
                    outspecies.append(outspecies_leaf.name)
    return outspecies
