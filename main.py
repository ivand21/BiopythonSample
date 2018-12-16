from Bio import Phylo
import tree_utils
import click


@click.group()
def cli():
    pass


@cli.command()
@click.option('--file-path', help='File with phylo trees to parse')
@click.option('--file-format', help='Phylo tree file format', default='newick')
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


if __name__ == '__main__':
    cli()
