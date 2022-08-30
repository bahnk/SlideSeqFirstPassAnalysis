#!/usr/bin/python

from os import makedirs
from os.path import join

import pandas as pd
import re

path = "bin/notebooks/scanpy/scanpy.py"
directory = "bin/notebooks/scanpy/j2"

try:
	makedirs(directory)
except OSError as e:
	pass

# get file content
f = open(path, "r")
lines = f.readlines()
f.close()

# regex to parse the file
var_regex = re.compile("# j2 variable (\\w+): (\\w+)")
cell_regex = re.compile("# cell (\\w+) (\\w+) (\\w+): (\\w+)")

# the future jinja2 variable
variables = []

# the row for the cell bounds dataframe
rows = []

for i, line in enumerate(lines):

	# add future jinja2 variables
	m_var = var_regex.match(line)
	if m_var:
		variables.append( (m_var.group(1), m_var.group(2)) )

	# get bounds of cells
	m_cell = cell_regex.match(line)
	if m_cell:
		rows.append({
			"Type": m_cell.group(1),
			"Name": m_cell.group(4),
			"Hide": m_cell.group(2),
			"Scroll": m_cell.group(3),
			"Line": i
		})

# cells bounds dataframe
df = pd.DataFrame.from_records(rows)

new_rows = []
for info, d in df.groupby(["Type", "Name", "Hide", "Scroll"]):

	cell_type, name, hide, scroll = info

	# we don't want duplicates in cell names because it would produce duplicates
	# in template names
	assert d.shape[0] == 2, f"There are duplicates for {name} ({cell_type})"

	assert d.Hide.unique().size == 1, f"Hide is different for {name} ({cell_type})"
	assert d.Scroll.unique().size == 1, f"Scroll is different for {name} ({cell_type})"

	new_rows.append({
		**{"Type": cell_type, "Name": name},
		**{
			"Hide": hide,
			"Scroll": scroll,
			"Begin": d.Line.min(),
			"End": d.Line.max()
			}
	})

new_df = pd.DataFrame.from_records(new_rows).sort_values("Begin").reset_index()

for i, row in new_df.iterrows():

	cell_type = row.Type
	name = row.Name

	chunk = lines[ row.Begin : min(row.End+1, len(lines)) ][1:-1]

	# remove empty lines at the beginning
	begin = 0
	for j in range(len(chunk)):
		if chunk[j] == "\n":
			begin = j+1
		else:
			break

	# remove empty lines at the end
	end = len(chunk)
	for j in range(1, len(chunk)+1):
		if chunk[-j] == "\n":
			end = -j
		else:
			break

	chunk = chunk[begin:end]

	if cell_type == "markdown":
		chunk = [ line for line in chunk if line != '"""\n']

	src = []
	
	for line in chunk:
		l = line.replace("\t", "   ")
		for var_type, var_name in variables:
			if "string" == var_type:
				l = re.sub(f"j2_{var_name}", f"\"{{{{ {var_name} }}}}\"", l)
			else:
				l = re.sub(f"j2_{var_name}", f"{{{{ {var_name} }}}}", l)
		src.append(l)

	src = f"{{%- block {name} -%}}\n" + "".join(src) + "{%- endblock -%}\n"

	if cell_type == "markdown":
		t = "markdown"
		suf = "md"
	elif cell_type == "python":
		t = "code"
		suf = "py"
	elif cell_type == "r":
		t = "code"
		suf = "r"

	filename = "cell_{0}_{1}_{2}_{3}.{4}.{5}.j2".format(
		str(i).zfill(2),
		t,
		row.Hide,
		row.Scroll,
		name,
		suf
		)
	filepath = join(directory, filename)

	f = open(filepath, "w")
	f.write(src)
	f.close()

