from Bio import Phylo
import tree_utils
import click


@click.command()
@click.option('--file-path', help='File with phylo trees to parse')
@click.option('--file-format', help='Phylo tree file format', default='NEWICK')
@click.option('--draw/--no-draw', default=False, help='If used, tree will be drawn')
@click.option('--ASCII/--no-ASCII', default=False, help='If used, tree will be drawn- ascii format')
@click.option('--info/--no-info', default=False, help='If used, tree info will be shown')

def display(file_path, file_format, draw, ascii, info):
    trees = tree_utils.load(file_path, file_format)
    if trees is not None:
        for tree in trees:
            if draw:                
                    Phylo.draw(tree)                
            if ascii:
                    tree_utils.drawASCII(tree)                
            if info:
                    tree_utils.showInfo(tree)               
            else:
                    click.echo(tree.name)
                    click.echo(tree)
        


if __name__ == '__main__':
    display()
