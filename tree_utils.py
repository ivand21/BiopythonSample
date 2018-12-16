from Bio import Phylo
import click


def load(file_path, file_format):
    if file_path is None:
        click.echo('No file specified!')
        return None
    else:
        try:
            trees = Phylo.parse(file_path, file_format)
            return trees
        except:
            click.echo('File not found or invalid!')
            return None


def draw_ascii(tree):
    Phylo.draw_ascii(tree)


def get_max_depth(depths):
    max = 0
    for node in depths:
        if not (node.name is None):
            curr = depths[node]
            if max < curr or min == 0:
                max = curr
    return max


def get_min_depth(depths):
    min = 0
    for node in depths:
        if not (node.name is None):
            curr = depths[node]
            if min > curr or min == 0:
                min = curr
    return min


def show_info(tree):
    non_terminals_count = len(tree.get_terminals())
    terminal_count = tree.count_terminals()
    depths = tree.depths()
    if tree.name is not None:
        print("Name: ", tree.name)
    print("Non terminals: ", tree.get_nonterminals())
    print("Non terminals count: ", non_terminals_count)
    print("Terminals: ", tree.get_terminals())
    print("Terminals count: ", terminal_count)
    print("Total branch length: ", tree.total_branch_length())
    print("Total node count: ", non_terminals_count + terminal_count)
    print("Is bifurcating: ", tree.is_bifurcating())
    print("Is preterminal: ", tree.is_preterminal())
    print("Max depth: ", get_max_depth(depths))
    print("Min depth: ", get_min_depth(depths))
