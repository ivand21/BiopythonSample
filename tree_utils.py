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
        
def drawASCII(tree):
    Phylo.draw_ascii(tree)
    
def getMaxDepth(depths):
    max=0
    for node in depths:
        if not (node.name is None):
            curr = depths[node]
            if (max< curr or min==0):
                max = curr
    return max       
    
def getMinDepth(depths):
    min=0
    for node in depths:
        if not (node.name is None):
            curr = depths[node]
            if (min> curr or min==0):
                min = curr
    return min
 
def showInfo(tree):
    nonTerminalsCount = len(tree.get_terminals())
    terminalCount = tree.count_terminals()
    depths = tree.depths()
    print("Non terminals:",tree.get_nonterminals())
    print("Non terminals count:",nonTerminalsCount)
    print("Terminals:",tree.get_terminals())
    print("Terminals count: ",terminalCount)
    print("Total branch length: ",tree.total_branch_length())
    print("Total node count:",nonTerminalsCount+terminalCount )
    print("Is bifurcating: ",tree.is_bifurcating())
    print("Is preterminal",tree.is_preterminal())
    print("Max depth: ", getMaxDepth(depths))
    print("Min depth: ", getMinDepth(depths))
