from Bio import Phylo
from Bio.Phylo.Consensus import *
from Bio.Nexus import Trees
import tree_utils
import click
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
import dendropy
from dendropy import treecalc
from dendropy.calculate import treecompare
from Bio.PDB import *
import sys


@click.group()
def cli():
    pass


@cli.command()
@click.option('-i', '--file-path', help='File with phylo trees to parse')
@click.option('-f', '--file-format', help='Phylo tree file format', default='newick')
@click.option('--draw/--no-draw', default=False, help='If used, tree will be drawn')
@click.option('--ascii/--no-ascii', default=False, help='If used, tree will be drawn- ascii format')
@click.option('--info/--no-info', default=False, help='If used, tree info will be shown')
def display(file_path, file_format, draw, ascii, info):
    trees = tree_utils.load(file_path, file_format)
    if trees is not None:
        for tree in trees:
            is_drawn = False
            if draw:
                Phylo.draw(tree)
                is_drawn = True
            if ascii:
                tree_utils.draw_ascii(tree)
                is_drawn = True
            if info:
                tree_utils.show_info(tree)
                is_drawn = True
            if not is_drawn:
                click.echo(tree)
    else:
        click.echo('Invalid or not existing file')


@cli.command()
@click.option('-i', '--file-path', help='File with phylo trees to parse')
@click.option('-f', '--file-format', help='Phylo tree file format', default='newick')
@click.option('-a', '--algorithm', help='Algorithm used for generating consensus tree', default='strict')
@click.option('-m', '--majority-cutoff', default=0.5, help='Cutoff value for majority consensus algorithm', type=float)
@click.option('-o', '--output-file', help='Output file path - generated tree will be saved to this file')
def consensus(file_path, file_format, algorithm, majority_cutoff, output_file):
    trees = tree_utils.load(file_path, file_format)
    if algorithm.lower() == 'strict':
        consensus_tree = strict_consensus(trees)
    elif algorithm.lower() == 'majority':
        consensus_tree = majority_consensus(trees, majority_cutoff)
    elif algorithm.lower() == 'adam':
        consensus_tree = adam_consensus(trees)
    else:
        click.echo('Invalid algorithm!')
        return
    click.echo('Drawing consensus tree...')
    Phylo.draw(consensus_tree)
    if output_file is not None:
        Phylo.write(consensus_tree, output_file, file_format)


@cli.command()
@click.option('-i', '--file-path', help='File with phylo trees to parse')
@click.option('-f', '--file-format', help='Phylo tree file format', default='newick')
@click.option('-i2', '--file-path2', help='File with phylo trees to parse')
def nearest(file_path, file_format, file_path2):
    # Read the sequences and align
    aln = AlignIO.read(file_path, file_format)

    # Print the alignment
    print(aln)

    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)

    # Print the distance Matrix
    print('\nDistance Matrix\n===================')
    print(dm)

    # Construct the phylogenetic tree using UPGMA algorithm
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    # Draw the phylogenetic tree
    Phylo.draw(tree)

    # Print the phylogenetic tree in the terminal
    print('\nPhylogenetic Tree\n===================')
    Phylo.draw_ascii(tree)


@cli.command()
@click.option('-i', '--file-path', help='File with phylo trees to parse')
@click.option('-f', '--file-format', help='Phylo tree file format', default='newick')
@click.option('-i2', '--file-path2', help='File with phylo trees to parse')
def distance(file_path, file_format, file_path2):
    taxon_namespace = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get_from_path(file_path, file_format,taxon_namespace=taxon_namespace)
    tree2 = dendropy.Tree.get_from_path(file_path2, file_format,taxon_namespace=taxon_namespace)
    sym_diff = treecompare.symmetric_difference(tree1, tree2)
    euc_dis = treecompare.euclidean_distance(tree1, tree2)
    false_pos = treecompare.false_positives_and_negatives(tree1,tree2)
    robinson_dis = treecompare.robinson_foulds_distance(tree1,tree2)
    print("Symetric difference: ",sym_diff)
    print("Robinson Foulds distance: ",robinson_dis)
    print("False positives and negatives: ",false_pos)
    print("Euclidean distance: ",euc_dis)
    
if __name__ == '__main__':
    cli()
