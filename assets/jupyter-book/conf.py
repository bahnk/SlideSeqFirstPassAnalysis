#!/usr/bin/python

from os import listdir
from os.path import join

import json
import pandas as pd
import re
import sys
import yaml

#################################
def main(in_dir, out_dir, args):#
#################################

	tags = [
		("seurat", "Seurat"),
		("rctd", "RCTD"),
		("sparkx", "SPARK-X"),
		("scanpy", "Scanpy"),
		("destvi", "DestVI")
	]
	
	############################################################################
	
	regex = re.compile("^(?P<analysis>[^\\.]+)\\.(?P<sample>.*)\\.ipynb")
	
	rows = []
	
	for filename in listdir(in_dir):
	
		if not re.match(".*\.ipynb", filename):
			continue
	
		m = regex.match(filename)
	
		if m:
			row = {
				"Analysis": m.group("analysis"),
				"Sample": m.group("sample"),
				"File": filename
			}
			rows.append(row)
	
	df = pd.DataFrame.from_records(rows)
	
	############################################################################
	
	toc = {}
	toc["format"] = "jb-book"
	toc["root"] = "index"
	toc["parts"] = []
	
	for tag, title in tags:
	
		d = df.loc[ df.Analysis == tag ].sort_values("Sample")
	
		part = {}
		part["caption"] = title
		part["numbered"] = False
		part["chapters"] = []
	
		for notebook in d.File:
			part["chapters"].append( {"file": notebook} )
	
		if len( part["chapters"] ) > 0:
			toc["parts"].append(part)
	
	f = open(join(out_dir, "_toc.yml"), "w")
	f.write( yaml.dump(toc, sort_keys=False) )
	f.close()
	
	############################################################################
	
	config = {}
	config["title"] = "Project " + args["project"]
	config["author"] = args["scientist"] + " & BABS"
	config["logo"] = "logo.jpg"
	config["execute"] = {"execute_notebooks": "off"}
	config["sphinx"] = {
		"config": {
				"html_js_files": "https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js"
		}
	}
	
	f = open(join(out_dir, "_config.yml"), "w")
	f.write( yaml.dump(config, sort_keys=False) )
	f.close()
	############################################################################

#in_dir = "results/notebooks"
#out_dir = "tmp"
#args = {"project": "Name", "scientist": "Random"}

if __name__ == "__main__":

	in_dir = sys.argv[1]
	out_dir = sys.argv[2]
	j = open(sys.argv[3], "r")
	args = json.load(j)
	j.close()

	main(in_dir, out_dir, args)

