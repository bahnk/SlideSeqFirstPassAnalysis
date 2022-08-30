#!/usr/bin/python2

from os import makedirs
from os.path import join

import json

directory = "bin/seurat"

samples = [
	"210611_03",
	"210611_11",
	"210611_15",
	"210709_19",
	"210709_26"
]

try:
	makedirs( join(directory, "conf") )
except OSError as e:
	pass

for sample in samples:

	conf = {}
	conf["name"] = sample
	conf["path_dge"] = f"../results/{sample}_dge"
	conf["path_spatial"] = f"../results/{sample}.csv"
	conf["path_outdir"] = f"results/seurat/{sample}"
	conf["clusters_resolution"] = f"0.05"

	f = open(join(directory, "conf", f"{sample}.json"), "w")
	json.dump(conf, f, indent=3)
	f.write("\n")
	f.close()

