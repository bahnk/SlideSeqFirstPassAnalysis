#!/usr/bin/python

import click
import nbformat
import re

from jinja2 import Template, BaseLoader, FileSystemLoader
from jinja2.environment import Environment
from json import dump
from os import listdir
from pathlib import Path
from sys import argv


###############################################################################
## Render()
###############################################################################
def Render(name, reference, puck, dge, outdir="tmp", ngenes=2000, test=False):
	"""Produces jupyter notebook."""

	template_dir = Path(argv[0]).parent / "j2"

	env = Environment()
	env.loader = FileSystemLoader(template_dir)

	cells = []

	for filename in sorted(listdir(template_dir)):

		if not filename.endswith("j2"):
			continue

		template = env.get_template(filename)

		template_name = list(template.blocks.keys())[0]
		cell_type = re.sub("cell_\\d+_([a-z]+)\..*", "\\1", filename)

		source = []
		lines = template.render(**locals()).split("\n")
		for i, line in enumerate(lines, start=1):
			if i < len(lines):
				source.append( line.replace("\t", "   ") + "\n" )
			else:
				source.append( line.replace("\t", "   ") )

		cell = {
			"cell_type": cell_type,
			"execution_count": None,
			"metadata": {},
			"outputs": [],
			"source": source
		}

		cells.append(cell)
	
	notebook = {
		"metadata": {
			"kernelspec": {
				"display_name": "Python 3 (ipykernel)",
				"language": "python",
				"name": "python3"
			}
		},
		"nbformat": 4,
		"nbformat_minor": 5,
		"cells": cells
	}

	nbformat.validate(notebook)

	return notebook
	###########################################################################

@click.command()
@click.option('--name', prompt="name", help="Sample name", default="sample")
@click.option('--reference', prompt="reference", help="Path to referene.")
@click.option('--puck', prompt="puck", help="Path to the puck.")
@click.option('--dge', prompt="dge", help="Path to the DGE matrix.")
@click.option('--outdir', prompt="outdir", help="Out directory.", default="tmp")
@click.option('--ngenes', prompt="ngenes", help="Var genes.", default=2000)
@click.option('--test', prompt="test", help="Test mode.", default=False)
@click.option('--nbook', prompt="nbook", help="Out.", default="nbook.ipynb")
def render(name, reference, puck, dge, outdir, ngenes, test, nbook):
	args = {
		"name": name, 
		"reference": reference,
		"puck": puck,
		"dge": dge,
		"outdir": outdir,
		"ngenes": ngenes,
		"test": test
	}
	nbk = Render(**args)
	with open(nbook, "w") as j:
		dump(nbk, j, indent=3)

if __name__ == "__main__":
	render()

