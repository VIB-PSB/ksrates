import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import sys
from ete3 import Tree, TreeStyle, NodeStyle, Face, RectFace, TextFace, StaticItemFace
from ete3.treeview.qt4_face_render import _TextFaceItem
from PyQt5.QtCore import QRectF
from PyQt5.QtWidgets import QGraphicsRectItem
from PyQt5.QtGui import QColor, QBrush, QFontMetrics
import logging
from statistics import mean
import ksrates.fc_rrt_correction as fcCorrect

# Filenames
_TREE_BRANCH_DISTANCES = "tree_{}_distances.pdf"
_TREE = "tree_{}.pdf"

class _LabelToggleItem(QGraphicsRectItem):
    # Empty item whose sole purpose is to toggle the display of the branch length label depending on whether it fits
    # in the width given by the branch line. Needs to be as a Face with position="branch-top" so that it has access
    # to the node item through the parent item hierarchy.
    # It determines the plotting width of the branch line of a node via its parent item and compares it to the plotting
    # width of the branch length label. If the width of the branch length label is larger than the width of branch line
    # it clears the text of the branch length label.

    def __init__(self, clearable_text_face, *arg, **karg):
        QGraphicsRectItem.__init__(self, *arg, **karg)
        self.clearable_text_face = clearable_text_face

    def boundingRect(self):
        # we do not plot anything
        return QRectF(0, 0, 0, 0)

    def paint(self, painter, option, widget):
        if self.clearable_text_face is not None:
            logging.debug(f'LabelToggleItem paint  \t{self.clearable_text_face.node.name}')
            logging.debug(f'  dist:                {self.clearable_text_face.node.dist}')
            logging.debug(f'  x:                   {self.x()}')
            logging.debug(f'  scene x:             {self.scenePos().x()}')
            logging.debug(f'  parent x:            {self.parentItem().x()}')
            # logging.debug(f'  parent scene x:      {self.parentItem().scenePos().x()}')
            # logging.debug(f'  grandparent x:       {self.parentItem().parentItem().x()}')
            # logging.debug(f'  grandparent scene x: {self.parentItem().parentItem().scenePos().x()}')

            # the parent item of our (empty) item is the node item and the
            # plotting width of the branch line is the x position of that
            # node item (the default filled circle) in the coordinates of
            # the parent item.
            distance_branch_width = round(self.parentItem().x())
            logging.debug(f'  dist branch width:   {distance_branch_width}')

            # painter.setPen(QColor("blue"))
            # painter.drawLine(0, 0, 0, 2)
            # painter.setPen(QColor("green"))
            # painter.drawLine(distance_branch_width, 0, distance_branch_width, 2)

            distance_label_width = self.clearable_text_face.get_text_width()
            logging.debug(f'  dist label width:    {distance_label_width}')
            if distance_label_width >= distance_branch_width:
                logging.debug(f'  dist label width >= dist branch width, clearing distance label!')
                self.clearable_text_face.clear_text()


class _ClearableTextFace(TextFace):
    """
    Text Face for which the text can be cleared.

    Adapted from and extends ete3.treeview.faces.TextFace

    :param text:     Text to be drawn
    :param ftype:    Font type, e.g. Arial, Verdana, Courier (default=Arial)
    :param fsize:    Font size, e.g. 10,12,6, (default=5)
    :param fgcolor:  Foreground font color. RGB code or color name in :data:`SVG_COLORS` (default="black")
    :param penwidth: Pen width used to draw the text (default=0)
    :param fstyle:   "normal" or "italic" (default="normal")
    """

    def __repr__(self):
        return "Clearable Text Face [%s] (%s)" % (self._text, hex(self.__hash__()))

    def __init__(self, text, ftype="Arial", fsize=10, fgcolor="black", penwidth=0, fstyle="normal"):

        Face.__init__(self)
        TextFace.__init__(self, text, ftype, fsize, fgcolor, penwidth, fstyle, False)

        self.type = "item"
        self.item = None

    def update_items(self):
        # uses protected member ete3.treeview.qt4_face_render._TextFaceItem, probably we should write our own
        self.item = _TextFaceItem(self, self.node, self.get_text())
        self.item.setFont(self._get_font())
        self.item.setBrush(QBrush(QColor(self.fgcolor)))

    def get_text_width(self):
        fm = QFontMetrics(self._get_font())
        text_width = fm.width(self.get_text())
        logging.debug(f'ClearableTextFace get_text_width  \t{self.node.name}')
        logging.debug(f'  dist:       {self.node.dist}')
        logging.debug(f'  text:       {self.get_text()}')
        logging.debug(f'  text width: {text_width}')
        return text_width

    def clear_text(self):
        if self.item is not None:
            logging.debug(f'ClearableTextFace clear_text  \t{self.node.name}')
            self.item.setText("")



def counts_expected_line_number_in_correction_table(species, tree, latin_names):
    """
    Returns the expected number of lines of the complete rate-adjustment table:
    there is one line per species in the tree, excluding the species of
    interest itself and the species of the deepest outgroup branching from the root
    and that cannot be adjusted.

    :param species: node of the focal species
    :param tree: input tree object
    :param latin_names: a dictionary-like data structure that associates each informal species name to its latin name
    :return expected_species_in_correction_table: list of expected species in the adjustment table
    """
    species_node = get_species_node(species, tree)
    species_history = get_species_history(species_node)
    expected_line_number_in_correction_table = len(species_history[-2]) - 1  

    leaves = species_history[-2].get_leaves()
    expected_species_in_correction_table = []
    for leaf in leaves:
        leaf = leaf.name
        if leaf != species:
            expected_species_in_correction_table.append(latin_names[leaf])

    return expected_species_in_correction_table


def find_missing_pairs_for_tree_rates(tree, species, species_history, latin_names):
    """
    Finds all species pairs in the tree that are not strictly needed for the adjustment of
    divergence lines, but that are instead needed to be able to compute the branch length 
    for all the branches present in the tree (with the exception of the base-most species).
    
    :param tree: input tree object
    :param species: node of the focal species
    :param species_history: the list of ancestor nodes of the focal species; includes the focal species and goes up to the root included
    :param latin_names: a dictionary-like data structure that associates each informal species name to its latin name
    :return missing_pairs_with_latin_names: list of lists, contains the missing necessary species name (both latin and informal names)
    :return missing_pairs: list of lists, contains the missing informal species names
    """
    missing_pairs_with_latin_names = []
    missing_pairs = []
    for divergence_node in species_history[:-1]:
        sister_node = divergence_node.get_sisters()[0]
        for node in sister_node.traverse():
            if not node.is_leaf():
                children = node.get_children()
                leaf1 = children[0].get_leaves()[0]
                leaf2 = children[1].get_leaves()[0]
                if node not in missing_pairs_with_latin_names:
                    latin1, latin2 = latin_names[leaf1.name], latin_names[leaf2.name]
                    sorted_latin_tag = "_".join(sorted([latin1, latin2], key=str.casefold))
                    missing_pairs_with_latin_names.append([sorted_latin_tag, sorted([leaf1.name, leaf2.name], key=str.casefold)])
                    missing_pairs.append(sorted([leaf1.name, leaf2.name], key=str.casefold))
    return missing_pairs_with_latin_names, missing_pairs
            

def reorder_tree_leaves(tree, species):
    """
    :param tree: the original tree object
    :param species: the name of the focal species
    :return: a new equivalent tree object with the focal species as top leaf
    """
    species_node = get_species_node(species, tree)[0]
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


def label_leaves_with_latin_names(tree, latin_names):
    """
    Changes the name attribute of each leaf node object from the informal name to the latin name of the species and
    adds a TextFace to each leaf node object that displays it.

    :param tree: the Newick tree
    :param latin_names: a dictionary-like data structure that associates each informal species name to its latin name
    """
    for leaf_node in tree.get_leaves():
        leaf_node.name = latin_names[leaf_node.name]
        # The whitespace after the name in the string ("{leaf_node.name} ") is a padding from the branch line
        leaf_node.add_face(TextFace(f"{leaf_node.name} ", ftype="Arial", fsize=18, fstyle="italic"), column=0,
                           position="branch-right")


def node_and_branch_style(tree):
    """
    Sets the style of tree branches (solid lines) and their default size.
    In the final picture, it will be the style of the branch lines with known lengths.

    :param tree: the Newick tree
    """
    nstyle_branch_circle = NodeStyle()
    nstyle_branch_circle["hz_line_width"] = 2
    nstyle_branch_circle["vt_line_width"] = 2
    nstyle_branch_circle["size"] = 0
    for node in tree.traverse():
        node.set_style(nstyle_branch_circle)
        node.dist = 0.25
        if node.is_root():
            # make root branch shorter
            node.dist = 0.125
        # add a type feature to each node, so we can later label specific nodes with
        # specific types, see e.g. labeling_internal_nodes()
        # (if ete's TreeNode class would have a function has_feature() we wouldn't need this)
        node.add_feature("type", "node")


def unknown_branch_len_style(unknown_node, root=False):
    """
    Updates the style of the branch lines with unknown length (dashed lines).

    :param unknown_node: the node whose branch length is unknown (the branch that comes out of it and
                         goes up to its parent node)
    :param root:         whether the unknown_node is the root node (default=False)
    """
    unknown_branch_style = NodeStyle()
    unknown_branch_style["hz_line_type"] = 1
    if root:
        unknown_branch_style["vt_line_type"] = 1
    unknown_branch_style["hz_line_width"] = 2
    unknown_branch_style["vt_line_width"] = 2
    unknown_branch_style["size"] = 0
    unknown_node.set_style(unknown_branch_style)


def get_species_node(species_name, tree):
    """
    :param species_name: the name of the focal species as it appears in the Newick tree
    :param tree: the Newick tree 
    :return species_node: the tree node object associated to the focal species
    """
    species_node = tree.search_nodes(name=species_name)
    if len(species_node) == 0:
        logging.error(f"The focal species ({species_name}) was not found in the provided tree. "
                               f"Exiting.")
        sys.exit(1)
    return species_node


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
        sys.exit(logging.error("One or more phylogenetic relationships in the tree are unresolved (>2 branches). Exiting."))
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
            # Limiting the amount of outgroups to the closest ones, if required in configuration file
            for outspecies_leaf in outspecies_leaves:
                if isinstance(max_num_outspecies, int):
                    if len(outspecies) < max_num_outspecies:
                        outspecies.append(outspecies_leaf.name)
                else:
                    outspecies.append(outspecies_leaf.name)
    return outspecies


def get_branch_length_and_errorbox(species, ancestor_node, correction_table, consensus_strategy_for_multi_outgroups, latin_names, rate_species_dict, rate_sister_dict):
    """
    NOTE: the average_peak_of_divergence_event, margin_error_box and error_text are not used later in the code, for now; they will be removed or re-integrated.
    The main purpose of this function is to update the two dictionaries collecting the branch-specific Ks contributions present in the correction_table (rate_species_dict and
    rate_sister_dict).

    Given an ancestor node belonging to the species history, it takes into account all the sister species that diverged at that node.
    For each sister, takes the adjusted divergence Ks value; then it makes an average adjusted Ks value to have a single representative Ks value for the node.
    For each sister, takes the error margins for the divergence (SD); then it considers the lowest and the highest error margins as the margins for the divergence. 
    The functions takes into account the user choice on how to deal with multiple outgroup.
    Returns the (mean) Ks value for the current divergence node in the tree, and the left and right margin values for delimiting an error box around it.
    
    :param species: the current focal species
    :param ancestor_node: one of the internal nodes belonging to the path that goes from the tree root to the focal species
    :param correction_table: adjustment results in DataFrame format (contains both possible types of consensus strategy for how to deal with multiple outgroups)
    :param consensus_strategy_for_multi_outgroups: user choice about which consensus strategy to use when dealing with multiple outgroups
    :param latin_names: a dictionary-like data structure that associates each informal species name to its latin name
    :param rate_species_dict: empty dictionary that will associate the focal species with its branch-specific Ks contribution at each divergence
    :param rate_sister_dict: empty dictionary that will associate each sister species with its own branch-specific Ks contribution
    :return: average_peak_of_divergence_event, adjusted Ks value of the current divergence (it is a mean value in case of multiple species diverged at that node with the focal species)  
    :return: margin_error_box, dictionary containing the smallest left error margin and the highest right error margin for the divergence when considering all the species diverging from that node
    :return: error_text, same as margin_error_box but in string format (left and right margins within brackets)
    """ 
    ancestor_node_leaves = get_sister_species_of_a_node(ancestor_node) # getting the sisters involved in the current divergence node

    equivalent_peak_list = [] # divergence peaks from all the species that divergence in the current internal node 
    errors_dict = {}

    # Depending on the user choice for dealing with multiple outgroups, we plot the results obtained by using
    # either the "mean among outgroups" strategy or the "best outgroup" strategy.
    # The parameter comes from the config file and it is given as argument to this function ("consensus_strategy_for_multi_outgroups")
    if consensus_strategy_for_multi_outgroups == "mean among outgroups":
        # Then consider the columns in the correction_table that are generated using the average method
        column_header_peak, column_header_SD, column_header_rate_species, column_header_rate_sister = "Adjusted_Mode_Mean", "Adjusted_Mode_Mean_SD", "Ks_Focal_Mean", "Ks_Sister_Mean"
    elif consensus_strategy_for_multi_outgroups == "best outgroup":
        # Then consider the columns in the correction_table that are generated using the best outgroup method
        column_header_peak, column_header_SD, column_header_rate_species, column_header_rate_sister = "Adjusted_Mode_Best", "Adjusted_Mode_Best_SD", "Ks_Focal_Best", "Ks_Sister_Best"

    # All the sister species from the current node provide equivalent data, because they all measure the same divergence (node)
    for sister in ancestor_node_leaves:
        latinSister = latin_names[sister]
        # Getting the adjusted peak (with SD) for the divergence with the current sister species
        corrected_peak = correction_table.loc[correction_table['Sister_Species'] == latinSister, [column_header_peak]]
        corrected_peak_float = corrected_peak.iat[0,0] # to convert from DataFrame type to Float type
        corrected_peak_sd = correction_table.loc[correction_table['Sister_Species'] == latinSister, [column_header_SD]]
        corrected_peak_sd_float = corrected_peak_sd.iat[0,0]
        equivalent_peak_list.append(corrected_peak_float)
        
        # Computing the left and right margins of the divergence for the current sister species 
        left_margin_error = corrected_peak_float - corrected_peak_sd_float
        right_margin_error = corrected_peak_float + corrected_peak_sd_float
        errors_dict[sister] = [left_margin_error, right_margin_error]

        # Get the branch length for the focal species at each divergence node
        rate_species = correction_table.loc[correction_table['Sister_Species'] == latinSister, [column_header_rate_species]]
        rate_species_float = rate_species.iat[0,0] # to convert from DataFrame type to Float type
        rate_species_dict[species] = rate_species_float
        # Get the branch length of the sister species (it's the branch-specific Ks contribution of the species) and adding it to a dictionary
        rate_sister = correction_table.loc[correction_table['Sister_Species'] == latinSister, [column_header_rate_sister]]
        rate_sister_float = rate_sister.iat[0,0] # to convert from DataFrame type to Float type
        rate_sister_dict[sister] = rate_sister_float
    
    # After collecting data of each sister species for the current divergence node, we need a single value to represent this node in the tree.
    # The average among the equivalent peaks of all the sisters is taken.
    # Note: this is not used for now.
    average_peak_of_divergence_event = mean(equivalent_peak_list)

    # The single mean value just computed is just a representative number of a broader range given by the error margins.
    # Let's take into account the margins of all the sister species, and consider the overall information from the lowest margin to the highest margin.
    # Getting the left and right margins to make an "error box" around the divergence node 
    margin_error_box = {}
    smallest_lef_margin = float("inf")
    highest_right_margin = -float("inf")
    for sister in errors_dict.keys():
        smallest_lef_margin = min([smallest_lef_margin, errors_dict[sister][0]])
        highest_right_margin = max([highest_right_margin, errors_dict[sister][1]])
    margin_error_box["left"] = round(smallest_lef_margin, 2)
    margin_error_box["right"] = round(highest_right_margin, 2)
    error_text = [margin_error_box["left"], margin_error_box["right"]]

    return average_peak_of_divergence_event, margin_error_box, error_text


def get_rates_from_current_analysis(rate_dict, correction_table, species, species_history, latin_names):
    """
    Gets the branch-specific Ks contributions obtained from the current analysis, namely the ones computed between the focal species\\
    and each of the other species. Updates rate_dict with such branch-specific Ks contributions.

    :param rate_dict: empty dictionary that will collect the known branch-specific Ks contributions per tree node
    :param correction_table: adjustment data from which the branch-specific Ks contributions coming from the current analysis will be extracted
    :param species: the focal species of the current analysis
    :param species_history: the list of ancestor nodes of the focal species; includes the focal species and goes up to the root included
    :param latin_names: dictionary associating the informal species names to their latin names
    """
    already_used_leaves = []
    for node in species_history[1:-1]:
        rate_dict[node] = {}
        leaves = node.get_leaves()
        for leaf in leaves:
            if leaf.name != species and leaf.name not in already_used_leaves:
                rate_sister = correction_table.loc[correction_table['Sister_Species'] == latin_names[leaf.name], ["Ks_Sister_Mean"]]
                rate_sister_float = rate_sister.iat[0,0] # to convert from DataFrame type to Float type
                rate_dict[node][leaf.name] = rate_sister_float
                already_used_leaves.append(leaf.name)


def get_rates_from_ortholog_peak_db(rate_dict, sister_node, latin_names, ortholog_db, peak_stats, missing_ortholog_data_from_database):
    """
    It's possible that some branches in the tree can't be assigned a length equal to branch-specific Ks contributions based only on the data\\
    coming from the current rate-adjustment analysis. In this case, the code tries to compute the missing branch-specific Ks contributions on the spot\\
    with the relative rate test formulas by looking for the missing ortholog data in the ortholog peak database: \\
    in fact, other rate-adjustments based on other focal species may have already provided such missing ortholog peaks needed for the RRT formulas,\\
    or perhaps the user can decide to run separately from this analysis the wgd ortholog pipeline needed to get the missing ortholog peaks\\
    and then they can try again to obtain the tree figure with complete branch lengths.

    :param rate_dict: dictionary that collects the known branch-specific Ks contributions per tree node
    :param sister_node: the sister node (it's only one!) of the current ancestor node of the focal species (which belongs to species_history)
    :param latin_names: dictionary associating the informal species names to their latin names
    :param ortholog_db: filename/path to the ortholog peak database
    :param peak_stats: flag to specify whether the ortholog distribution peak is the mode or the median
    :param missing_ortholog_data_from_database: flag to state if one or more branches are lacking Ks contributions due to missing data in the ortholog peak database 
    :return: the updated flag (set to True if there are missing ortholog data and some branch lengths are unknown)
    """
    list_of_nodes = []
    for n in sister_node[0].traverse():
        list_of_nodes.append(n)

    for node in list_of_nodes:
        if not node.is_leaf():
            rate_dict[node] = {}
            outspecies_list = get_sister_species_of_a_node(node)
            outspecies_list.extend(get_outspecies_of_a_node(node, None))
            children = node.get_children()
            child1_leaves, child2_leaves = children[0].get_leaves(), children[1].get_leaves()
            
            node_rate_resolved = False
            for child1_leaf in child1_leaves:
                for child2_leaf in child2_leaves:
                    reduced_latin_names = {child1_leaf.name : latin_names[child1_leaf.name],  child2_leaf.name : latin_names[child2_leaf.name]}
                    sorted_names = sorted(reduced_latin_names.items(), key=lambda x: x[1].lower()) # alphabetically sorted according to value and case-insensitive
                    sister1, sister2 = sorted_names[0][0], sorted_names[1][0]

                    latinSister1, latinSister2 = sorted_names[0][1], sorted_names[1][1]
                    latinSister1_latinSister2 = "_".join([latinSister1, latinSister2])

                    # Check if the ortholog distribution peak of the two species is present in the database
                    # If not, warn the user that it has to be computed separately to have all branch length
                    # equal to branch-specific Ks contributions
                    try:
                        __ = ortholog_db.at[latinSister1_latinSister2, 'Mode']
                    except Exception:
                        missing_ortholog_data_from_database = True

                    # Check if there are ortholog data in database to use a species as an outgroup for the two leaves
                    # There should always be at least the focal species of the current analysis, except if deleted from DB for some reasons
                    list_of_successful_outspecies = []
                    for outspecies in outspecies_list:
                        latinSister1_latinOut = "_".join(sorted([latinSister1, latin_names[outspecies]]))
                        latinSister2_latinOut = "_".join(sorted([latinSister2, latin_names[outspecies]]))

                        try:
                            __ = ortholog_db.at[latinSister1_latinOut, 'Mode']
                            __ = ortholog_db.at[latinSister2_latinOut, 'Mode']
                            list_of_successful_outspecies.append(outspecies)
                        except Exception:
                            pass

                        try: # Trying to compute the branch-specific Ks contributions (relative rate test formulas)
                            rate_species, __, rate_sister, __ = fcCorrect.decompose_ortholog_ks(ortholog_db, latinSister1_latinSister2, latinSister1_latinOut, latinSister2_latinOut, peak_stats)
                            if sister1 not in rate_dict[node]:
                                rate_dict[node][sister1] = rate_species
                            if sister2 not in rate_dict[node]:
                                rate_dict[node][sister2] = rate_sister
                            node_rate_resolved = True
                            break
                        except Exception:
                            pass
                    
                    if list_of_successful_outspecies == []:
                        missing_ortholog_data_from_database = True
                    
                    if node_rate_resolved: # no need to investigate other leaves of child2 node;
                        break              # the rate (length) of these branches can be computed with at least one species pair 
                if node_rate_resolved:     # same thing here, no need to investigate other leaves of child1 node
                    break 

        for leaf in node.get_leaves():
            parent = node.up
            try:
                if not node.is_leaf():
                    branch_length = rate_dict[parent][leaf.name] - rate_dict[node][leaf.name]
                    node.dist = branch_length
                    draw_branch_length_label(node, known_distance=True)
                    break
                else:
                    branch_length = rate_dict[parent][leaf.name]
                    node.dist = branch_length
                    draw_branch_length_label(node, known_distance=True)
                    break
            except Exception: # Not enough data to compute the rate / set a length
                node.dist = 10 # impossible number to flag an unknown length
                draw_branch_length_label(node, known_distance=False)
                unknown_branch_len_style(node)

            if not node.is_leaf(): # just to avoid that there are multiple AttrFace labels on the internal node branches
                break
    return missing_ortholog_data_from_database


def draw_branch_length_label(node, known_distance):
    """
    Places a label with the branch length on top of a node's branch line.

    :param node: the node to place the label on
    :param known_distance: flag to state if the branch length of this node is known or unknown
    """
    # Actually places the real label on top of the branch line and an empty ghost label below it,
    # both using "float" positions. This is a trick, because otherwise the label would overlap with line.
    # I didn't use the "branch-top" position because that mode automatically adapts and extends the
    # branch line to the width of the label with a dotted line, visually distorting the actual branch lengths.

    # Generate the branch length label/face based on whether the distance is known
    if known_distance:
        distance = f"{round(node.dist, 2)}"
    else:
        distance = ""
    distance_face = _ClearableTextFace(distance)
    distance_face.margin_bottom = 3
    if node.is_leaf():
        distance_face.margin_right = 0
    else:
        distance_face.margin_right = 1

    # add the branch length label as a "float" so it does not cause
    # dotted extension of the branch length line to make the label fit
    node.add_face(distance_face, column=0, position="float")
    # add another ghost/empty label so the branch length label sits above
    # and not on top of the line
    node.add_face(TextFace(" ", fsize=8), column=0, position="float")

    if not node.type == "species_of_interest" and not node.type == "ancestor":
        # Add an empty face with a _LabelToggleItem at "branch-top" to extract
        # the plot width of the branch line and enable/disable the branch length label
        # depending on whether it fits in this width
        empty_label_toggle_face = StaticItemFace(_LabelToggleItem(distance_face))
        node.add_face(empty_label_toggle_face, column=0, position="branch-top")


def adapt_unknown_branch_length(tree):
    """
    Sets a default length for the branch lines with unknown real length.
    This default values is the average among all known lengths, so that 
    the lines with unknown length are more or less in scale with the others.
    The root length is smaller (1/2 average length) just to save some space in the figure.
    In the code, a branch with unknown length was flagged with a distance attribute set as an impossible number, 10,
    that's why the function ignores all branches with distance equal to 10.
    
    :param tree: the tree object
    """
    distance_list = []
    for node in tree.traverse():
        if node.dist != 10:  # impossible number to flag an unknown length
            distance_list.append(node.dist)
    avg_dist = mean(distance_list)
    for node in tree.traverse():
        if node.is_root():
            # make root branch shorter
            node.dist = avg_dist / 4
            unknown_branch_len_style(node, root=True)   
        if node.dist == 10:
            if node.up.is_root() and node.is_leaf():
                # make branch of single outgroup longer
                node.dist = avg_dist * 2
            else:
                node.dist = avg_dist
            unknown_branch_len_style(node)


def plotting_tree(species, latin_names, original_tree, correction_table, consensus_strategy_for_multi_outgroups, ortholog_db, peak_stats, nextflow_flag):
    """
    Generate a PDF figure of the input tree with branch lengths equal to Ks distances.
    If it is not possible to compute the branch length for a branch, the branch line is dashed. This happens when some\\
    ortholog data to compute the branch-specific Ks contribution are missing.

    :param species: the current focal species
    :param latin_names: a dictionary-like data structure that associates each informal species name to its latin name
    :param original_tree: Newick tree format of the phylogenetic tree among the involved species
    :param correction_table: adjustment results in DataFrame format (contains both possible types of consensus strategy for how to deal with multiple outgroups)
    :param consensus_strategy_for_multi_outgroups: user choice about which consensus strategy to use when dealing with multiple outgroups
    :para ortholog_db: ortholog peak database used to get ortholog data for the relative rate test; if not available, will be ignored
    :param peak_stats: flag to specify whether the ortholog distribution peak is the mode or the median
    :param nextflow_flag: boolean flag to state whether the script is run in the Nextflow pipeline or not
    """
    # Get an equivalent tree where the focal species is the top leaf
    tree = reorder_tree_leaves(original_tree, species)
    node_and_branch_style(tree)
    
    species_node = get_species_node(species, tree)
    
    labeling_internal_nodes(species_node)
    species_history = get_species_history(species_node)
    rate_species_dict, rate_sister_dict = {}, {}

    for ancestor_node in species_history[:-2]:
        # NOTE: at the moment the following function is only used to fill in the dictionaries of branch-specific Ks contributions
        average_peak_of_divergence_event, margin_error_box, error_text = get_branch_length_and_errorbox(species, ancestor_node,
                                                                                                        correction_table, consensus_strategy_for_multi_outgroups,
                                                                                                        latin_names, rate_species_dict, rate_sister_dict)

        # Adding the branch length to the focal species node, otherwise it lacks it
        if ancestor_node.name == species:
            ancestor_node.dist = rate_species_dict[species]
            draw_branch_length_label(ancestor_node, known_distance=True)

        # Adding as TextFaces both the divergent Ks of the node (as mean) and the error range (left-most and right-most boundaries)
        divergence_node = ancestor_node.up # getting parent node, where the current divergence takes place 
        divergence_node.add_feature("rate_species", rate_species_dict[species])
        divergence_node.add_feature("avg_peak", round(average_peak_of_divergence_event, 2))
        divergence_node.add_feature("margins", f"({error_text[0]}, {error_text[1]})")
        ### divergence_node.add_face(AttrFace("margins", fsize=5), column=0, position="branch-right") [ NOT USED FOR NOW ]

    # Setting the branch length of the nodes belonging to the speciation history of the focal species
    for divergence_node in species_history[1:]:
        parent_node = divergence_node.up
        try:
            divergence_node.dist = round(parent_node.rate_species - divergence_node.rate_species, 3)
            draw_branch_length_label(divergence_node, known_distance=True)
        except Exception:
            divergence_node.dist = 10 # impossible number to flag an unknown length
            draw_branch_length_label(divergence_node, known_distance=False)
            unknown_branch_len_style(divergence_node)

    if ortholog_db.empty: # branch-specific Ks contributions can be obtained only from adjustment_tables
        logging.info("Getting branch-specific Ks contributions from rate-adjustment table data")
    else: # if the ortholog DB is available, we can try to compute the branch-specific Ks contributions from there too
        logging.info("Getting branch-specific Ks contributions from rate-adjustment table data")
        logging.info("Computing branch-specific Ks contributions from ortholog peak data in database by applying principles of the relative rate test")

    rate_dict = {}
    get_rates_from_current_analysis(rate_dict, correction_table, species, species_history, latin_names)

    # Setting the branch length of the other remaining nodes
    missing_ortholog_data_from_database = False
    missing_ortholog_data_from_correction_table = False

    for node in species_history[:-1]:
        sister_node = node.get_sisters() # is a list containing the sister NODE (it's only ONE node)

        if not ortholog_db.empty: # if there is an ortholog database that can help with computing the missing branch lengths
            if len(sister_node[0].get_leaves()) > 1:
                missing_ortholog_data_from_database = get_rates_from_ortholog_peak_db(rate_dict, sister_node, latin_names, 
                                                                                      ortholog_db, peak_stats, 
                                                                                      missing_ortholog_data_from_database)
            else:
                if sister_node[0].name in rate_sister_dict.keys(): # if leaf has known length
                    sister_node[0].dist = rate_sister_dict[sister_node[0].name]
                    draw_branch_length_label(sister_node[0], known_distance=True)
                else: # if the leaf has unknown length
                    sister_node[0].dist = 10 # impossible number to flag an unknown length
                    draw_branch_length_label(sister_node[0], known_distance=False)
                    unknown_branch_len_style(sister_node[0])

        else: # if ortholog database not available (the variable was previously set as an empty dataframe)
            if len(sister_node[0].get_leaves()) > 1:
                missing_ortholog_data_from_correction_table = True # correction_tables is not enough to know all branch lengths!
                sister_node[0].dist = 10 # impossible number to flag an unknown length
                draw_branch_length_label(sister_node[0], known_distance=False)
                unknown_branch_len_style(sister_node[0])
                for node in sister_node[0].get_descendants():
                    node.dist = 10 # impossible number to flag an unknown length
                    draw_branch_length_label(node, known_distance=False)
                    unknown_branch_len_style(node)
            else:
                leaf = sister_node[0].get_leaves()[0] # there is only one leaf
                if leaf.name in rate_sister_dict.keys():
                    leaf.dist = rate_sister_dict[leaf.name]
                    draw_branch_length_label(leaf, known_distance=True)
                else: # if the leaf has unknown length
                    leaf.dist = 10 # impossible number to flag an unknown length
                    draw_branch_length_label(leaf, known_distance=False)
                    unknown_branch_len_style(leaf)
    
    # If the ortholog peak database is lacking some required data (must have been deleted by the user) or
    # if the peak database has been deleted and only the correction_table has been used for the branch contributions, gives a warning
    if missing_ortholog_data_from_database or missing_ortholog_data_from_correction_table:
        logging.warning("")
        logging.warning("One or more branch lengths are unknown (dashed line) due to missing ortholog distribution peak data")

    # If in Nextflow mode, tell the user to wait until the pipeline is finished in order to have all branch lengths
    if nextflow_flag:
        if missing_ortholog_data_from_database:
            logging.info(f"As soon as new ortholog data will become available, the tree branch lengths will be updated")
    # If manual mode, tell the user how to get a complete branch tree (probably they deleted some data in the peak database)
    else:
        if missing_ortholog_data_from_database or missing_ortholog_data_from_correction_table:
            logging.warning(f"It's necessary to run a new Nextflow (or manual) pipeline to complete the tree branch length information")
            
    label_leaves_with_latin_names(tree, latin_names)
    adapt_unknown_branch_length(tree)

    ts = TreeStyle()
    # ts.title.add_face(TextFace("  Input tree with branch length equal to Ks distances  ", ftype="Arial", fsize=18), column=0)
    ts.orientation = 1
    ts.branch_vertical_margin = 14
    ts.show_leaf_name = False # because there is a Face showing it
    ts.show_branch_length = False
    ts.margin_left = 25
    ts.margin_right = 25
    ts.margin_top = 25
    ts.scale = 200
    #ts.scale_length =  # to set a fixed scale branch length
    root_of_corrected_tree = species_history[-1]
    root_of_corrected_tree.render(os.path.join("rate_adjustment", f"{species}", f"{_TREE_BRANCH_DISTANCES.format(species)}"), w=4.5, units="in", tree_style=ts)


def plot_uncorrected_phylogeny(tree, species, latin_names, species_history):
    """
    Generates a PDF figure of the input tree with same length for all branches.

    :param tree: input tree from configuration file
    :param species: the current focal species
    :param latin_names: a dictionary-like data structure that associates each informal species name to its latin name
    :param species_history: the list of ancestor nodes of the focal species, including the focal species and going up to the root.
    """
    label_leaves_with_latin_names(tree, latin_names)
    node_and_branch_style(tree)
    ts = TreeStyle()
    # ts.title.add_face(TextFace("  Input phylogenetic tree", ftype="Arial", fsize=18), column=0)
    ts.orientation = 1
    ts.branch_vertical_margin = 14
    ts.show_leaf_name = False # because there is a Face showing it
    ts.show_branch_length = False
    ts.margin_left = 25
    ts.margin_right = 25
    ts.margin_top = 25
    ts.margin_bottom = 25
    ts.scale = 200
    ts.show_scale = False
    tree.render(os.path.join("rate_adjustment", f"{species}", f"{_TREE.format(species)}"), w=4.5, units="in", tree_style=ts)
