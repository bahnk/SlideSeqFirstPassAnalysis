#!/usr/bin/python

import json
import nbformat
import re

from jinja2 import Template, BaseLoader, FileSystemLoader
from jinja2.environment import Environment
from json import dump
from os import listdir
from pathlib import Path
from sys import argv

###############################################################################
## RenderCell()
###############################################################################
def RenderCell(environment, filename, args):
	"""Renders a cell with jinja2."""

	template = environment.get_template(filename)

	template_name = list(template.blocks.keys())[0]

	regex = re.compile(
		"cell_(?P<num>\\d+)_(?P<type>\\w+)_(?P<hide>\\w+)_(?P<scroll>\\w+).(?P<name>\\w+)\..+"
	)
	m = regex.match(filename)
	assert m, "Cannot parse template filename"

	source = []
	lines = template.render(**args).split("\n")
	for i, line in enumerate(lines, start=1):
		if i < len(lines):
			source.append( line.replace("\t", "   ") + "\n" )
		else:
			source.append( line.replace("\t", "   ") )

	if "markdown" == m.group("type"):

		cell = {
			"cell_type": m.group("type"),
			"metadata": {},
			"source": source
		}

		return cell

	if "code" == m.group("type"):

		tags = []
		if m.group("hide") == "hide":
			tags.append("hide-input")
		if m.group("scroll") == "scroll":
			tags.append("output_scroll")

		cell = {
			"cell_type": m.group("type"),
			"execution_count": None,
			"metadata": {
				"scrolled": True,
				"tags": tags
			},
			"outputs": [],
			"source": source
		}

		return cell
	############################################################################

###############################################################################
## Render()
###############################################################################
def Render(args):
	"""Produces jupyter notebook."""

	template_dir = Path(argv[0]).parent / "j2"

	env = Environment()
	env.loader = FileSystemLoader(template_dir)

	cells = []

	title = {
		"cell_type": "markdown",
		"metadata": {},
		"source": ["# " + args["name"]]
	}
	cells.append(title)


	for filename in sorted(listdir(template_dir)):

		if not filename.endswith("j2"):
			continue

		cell = RenderCell(env, filename, args)
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
	############################################################################

if __name__ == "__main__":

	json_file = open(argv[1], "r")
	args = json.load(json_file)
	json_file.close()

	nbk = Render(args)

	with open(argv[2], "w") as j:
		dump(nbk, j, indent=3)

