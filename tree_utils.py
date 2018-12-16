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
