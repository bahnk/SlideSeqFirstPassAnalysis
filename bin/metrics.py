#!/usr/bin/python

# coding: utf-8

import json
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp
import sys

#name = "sample"
#dge_path = "/camp/stp/babs/working/bahn/projects/slideseq/sequencing/miniseq/211018_MN01566_0015_A000H3KWV7/processing/results/puck_14R"
#spatial_path = "/camp/stp/babs/working/bahn/projects/slideseq/sequencing/miniseq/211018_MN01566_0015_A000H3KWV7/processing/results/puck_14R.csv"
#out_json = "tmp/test.json"

#############################################################
def compute_metrics(dge_path, spatial_path, name, out_json):#
#############################################################

	spatial = pd.read_csv(spatial_path).set_index("Barcode")
	adata = sc.read_10x_mtx(dge_path, var_names="gene_symbols")
	adata.obs = spatial.loc[ adata.obs.index ]
	
	X = adata.X.toarray() if sp.sparse.issparse(adata.X) else adata.X
	
	metrics = {
			"Name": name,
			"Total beads": int(adata.n_obs),
			"Total genes": int(adata.n_vars),
			"Total counts": int(X.sum()),
			"Detected genes": int((np.apply_along_axis(np.sum, 0, X) > 0).sum())
		}
	
	with open(out_json, "w+") as j:
		json.dump(metrics, j, indent=3)

if __name__ == "__main__":
	compute_metrics(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

