from Bio import Phylo
import tree_utils
import click


@click.command()
@click.option('--file-path', help='File with phylo trees to parse')
@click.option('--file-format', help='Phylo tree file format', default='NEWICK')
@click.option('--draw/--no-draw', default=False, help='If used, tree will be drawn')
def display(file_path, file_format, draw):
    trees = tree_utils.load(file_path, file_format)
    if trees is not None:
        if draw:
            for tree in trees:
                click.echo(tree.name)
                click.echo(tree)
        else:
            for tree in trees:
                tree.ladderize()
                Phylo.draw(tree)


if __name__ == '__main__':
    display()
