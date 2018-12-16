from Bio import Phylo
from utils import tree_utils
import click


@click.command()
@click.option('--file', help='File with phylo trees to parse')
@click.option('--format', help='Phylo tree file format', default='NEWICK')
@click.option('--draw/--no-draw', default=False, help='If used, tree will be drawn')
def display(file_path, file_format, draw):
    trees = tree_utils.load(file_path, file_format, draw)
    if draw:
        for tree in trees:
            click.echo(tree.name)
            click.echo(tree)
    else:
        for tree in trees:
            tree.ladderize()
            Phylo.draw(tree)
