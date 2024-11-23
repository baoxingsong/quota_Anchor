import sys
import logging
import warnings
warnings.filterwarnings('ignore', category=SyntaxWarning)
from . import ksrates
import pandas as pd

logger =logging.getLogger('main.trios')
MAX_NUM_OUT_SPECIES = 4


class Trios:
    def __init__(self, config_pra, parameter):
        self.focal_species = ""
        self.nwk_tree = ""
        self.outfile_trios_path = ""
        self.outfile_species_pair_file = ""

        for i in config_pra.sections():
            if i == 'trios':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])

        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)

    def get_trios(self, outfile_trios_path):
        focal_species = self.focal_species
        #TODO: color parameter may be need to add in ks_fitting module

        # divergence_colors = ["#f5493d", "#0022e0", "#b89f00", "#1b5f1c", "#f466e4", "#3593ab"]
        # num_required_colors = len(divergence_colors)
        tree = ksrates.get_newick_tree(self.nwk_tree)
        #TODO: this is print raw nwk tree
        logger.info("The following is raw nwk tree")
        print(tree)
        ksrates.check_integrity_newick_tree(tree)
        ordered_tree = ksrates.reorder_tree_leaves(tree, focal_species)
        logger.info("The following is ordered tree")
        print(ordered_tree)
        species_of_interest_node = ksrates.get_species_node(focal_species, ordered_tree)
        # get the list of ancestors (as tree node objects) in the lineage that lead to the focal species
        sp_history = ksrates.get_species_history(species_of_interest_node)

        if len(sp_history)-2 == 0:
            logger.error("")
            logger.error(f"Species [{focal_species}] has no outgroup in the provided Newick tree "
                        f"and the rate-adjustment can't be performed.")
            logger.error(f"Please add at least one outgroup species or change the focal species.")
            sys.exit(1)

        ksrates.labeling_internal_nodes(species_of_interest_node)
        # TODO: color parameter may be need to add in ks_fitting module

        # if len(divergence_colors) < num_required_colors:
        #     print("")
        #     print(f'Configuration file field "divergence_colors" is missing {num_required_colors - len(divergence_colors)} color(s) ' +
        #                 f"out of {num_required_colors} required for the analysis on focal species [{focal_species}]")
        #     print("Please add the missing color(s) and rerun the analysis")
        #     print("Exiting.")
        #     sys.exit(1)
        trios_array = []  # list of trios
        # TODO: txt file
        outfile_drawing_path = f"tree_{focal_species}.txt"
        with open(outfile_drawing_path, "w+") as outfile_drawing:
            outfile_drawing.write(f"Focal species: {focal_species}\n\n")
            
            node = 0
            while node < len(sp_history)-2:
                # the name label to be shown in the ASCII tree will start from 1 and not from 0
                outfile_drawing.write(f"Node {node+1}:\n")
                currentnode = sp_history[node]

                # GETTING SISTER SPECIES
                sisters = ksrates.get_sister_species_of_a_node(currentnode)
                outfile_drawing.write(f"Sister species:      {', '.join(sisters)}\n")

                # GETTING OUTSPECIES
                outspecies = ksrates.get_outspecies_of_a_node(currentnode, MAX_NUM_OUT_SPECIES)
                outfile_drawing.write(f"Outgroup species:    {', '.join(outspecies)}\n\n")

                # APPENDING TRIOS (a trio is composed of focal species, sister species and outgroup species)
                for s in sisters:
                    for o in outspecies:
                        trios_array.append([node+1, focal_species, s, o])
                node += 1

            print_tree = tree.get_ascii(attributes=["name"], show_internal=True)
            outfile_drawing.write(f"{print_tree}\n\n")

        logger.info(f"- Total number of trios: {len(trios_array)}")
        logger.info("")
        header = ["Node", "Focal_Species", "Sister_Species", "Out_Species"]
        df = pd.DataFrame(trios_array, columns=header)
        df.to_csv(outfile_trios_path, index=False, header=True)
        return df, trios_array

    @staticmethod
    def get_necessary_pairs(trios_array, species_pair_file):
        # THIS CODE JUST TAKES THE PAIRS THAT ARE NECESSARY FOR THE TRIOS
        species_pairs = []
        new_species_pairs = []
        for trio in trios_array:
            combos = [[trio[1], trio[2]], [trio[1], trio[3]], [trio[2], trio[3]]]
            for pair in combos:
                if pair not in species_pairs:
                    species_pairs.append(pair)
                    species_pairs.append([pair[1], pair[0]])
                    new_pair=[pair[0], pair[1], str(1), str(1)]
                    new_species_pairs.append(new_pair)
        
        # species_pairs = [sorted(x, key=str.casefold) for x in species_pairs]
        header = ["Species_1", "Species_2", "r_value", "q_value"]
        df = pd.DataFrame(new_species_pairs, columns=header)
        df.to_csv(species_pair_file, index=False, header=True)
        
        # TODO: ksrates can use unknown species pairs to construct tree branch length
        # species_pairs_unknown = []
        return species_pairs

    def run(self):
        logger.info("Init trios and the following parameters are config information")
        print()
        for key, value in vars(self).items():
            if key != "conf":
                print(key, "=", value)
        print()
        df, trios_array = self.get_trios(self.outfile_trios_path)
        self.get_necessary_pairs(trios_array, self.outfile_species_pair_file)
        logger.info(f"{self.outfile_trios_path} and {self.outfile_species_pair_file} generated done!")
