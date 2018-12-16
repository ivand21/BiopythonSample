from Bio import Phylo
from Bio.Phylo.Consensus import *
import tree_utils
import click


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


if __name__ == '__main__':
    cli()
