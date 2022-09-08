#!/usr/bin/env python

# coding: utf-8

# j2 variable string: name
# j2 variable string: path_dge
# j2 variable string: path_spatial 
# j2 variable string: path_reference
# j2 variable string: param_destvi_gene_identifier
# j2 variable integer: param_destvi_min_counts
# j2 variable integer: param_destvi_n_variable_genes
# j2 variable bool: param_destvi_test

j2_name = "sample1"
j2_path_dge = "test/data/sample1"
j2_path_spatial = "test/data/sample1.csv"
j2_path_reference = "test/data/reference"
j2_param_destvi_gene_identifier = "gene_symbols"
j2_param_destvi_min_counts = 10
j2_param_destvi_n_variable_genes = 2000
j2_param_destvi_test = True

###############################################################################
# cell markdown null null: intro
"""
We use [DestVI](https://www.nature.com/articles/s41587-022-01272-8) and a single cell referenece dataset to infere the cell types present on each spot.
"""
# cell markdown null null: intro

###############################################################################
# cell markdown null null: libraries
"""
Here are the libraries we need.
"""
# cell markdown null null: libraries

###############################################################################
# cell python nohide scroll: libraries

from os import makedirs, remove, rmdir
from os.path import exists, join
from scipy.sparse import csr_matrix, issparse
from scvi.model import CondSCVI, DestVI

import anndata
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import scanpy as sc
import session_info
# cell python nohide scroll: libraries

###############################################################################
# cell markdown null null: functions
"""
Some useful function to save AnnData objects and models.
"""
# cell markdown null null: functions

###############################################################################
# cell python nohide scroll: functions

###############################################
def save_anndata(anndata_object, num, suffix):#
###############################################
	if not issparse(anndata_object.X):
		anndata_object.X = csr_matrix(anndata_object.X)
	anndata_object.write_h5ad( join("output/h5ad", f"anndata.{num}.{suffix}.h5") )
	###########################################################################

#################################################
def save_model(model_object, directory, suffix):#
#################################################

	model_dir = join(directory, f"model_{suffix}")
	model_archive = join(model_dir, "model.pt")

	if exists(model_archive):
		remove(model_archive)

	if exists(model_dir):
		rmdir(model_dir)

	model_object.save(model_dir)
	###########################################################################

############################
def has_zeros(name, adata):#
############################
	genes = np.sum( 0 == np.sum(adata.X, axis=0) )
	cells = np.sum( 0 == np.sum(adata.X, axis=1) )
	print(f"{name}: {genes} zeros for genes and {cells} zeros for cells")

# cell python nohide scroll: functions

###############################################################################
# cell markdown null null: directories
"""
We create the output directory for this noteboook.
Every outputs will save there.
"""
# cell markdown null null: directories

###############################################################################
# cell python nohide scroll: directories

try:
	makedirs("output/h5ad")
except OSError as e:
	pass
# cell python nohide scroll: directories

###############################################################################
# cell markdown null null: scanpy_config
"""
Configuration of scanpy.
"""
# cell markdown null null: scanpy_config

###############################################################################
# cell python nohide scroll: scanpy_config
sc.settings.figdir = "output"
sc.set_figure_params(figsize=(4, 4), frameon=False)
# cell python nohide scroll: scanpy_config

###############################################################################
# cell markdown null null: title_reference
"""
## scRNA-seq reference dataset preprocessing
"""
# cell markdown null null: title_reference

###############################################################################
# cell markdown null null: reference_load
"""
### Loading

We load the single-cell reference dataset.
"""
# cell markdown null null: reference_load

###############################################################################
# cell python nohide scroll: reference_load
sc_adata = sc.read_10x_mtx(j2_path_reference, var_names=j2_param_destvi_gene_identifier)

args = {
	"filepath_or_buffer": join(j2_path_reference, "types.tsv.gz"),
	"header": None,
	"compression": "gzip"
}
types = pd.read_csv(**args).rename(columns={0: "Type"})
types.index = sc_adata.obs_names
sc_adata.obs = types

sc_adata
# cell python nohide scroll: reference_load

###############################################################################
# cell markdown null null: reference_symbols
"""
### Gene symbols

We use the gene symbols as IDs and put convert everything to upper case.
So, we also need to check that are no duplicates in gene symbols.
"""
# cell markdown null null: reference_symbols

###############################################################################
# cell python nohide scroll: reference_symbols
sc_adata.var_names = sc_adata.var_names.str.upper()
sc_adata = sc_adata[ : , ~ sc_adata.var_names.duplicated() ]
# cell python nohide scroll: reference_symbols

###############################################################################
# cell markdown null null: reference_zeros
"""
### Minimum counts

We need to remove the spots and the genes with no count.
"""
# cell markdown null null: reference_zeros

###############################################################################
# cell python nohide scroll: reference_zeros
sc_adata = sc_adata[ : , ~ np.all(sc_adata.X.toarray() == 0, axis=0) ]
sc_adata = sc_adata[ ~ np.all(sc_adata.X.toarray() == 0, axis=1) , : ]

save_anndata(sc_adata, "00", "raw")
# cell python nohide scroll: reference_zeros

###############################################################################
# cell markdown null null: reference_min_counts
"""
We remove the spots with too few UMIs.
"""
# cell markdown null null: reference_min_counts

###############################################################################
# cell python nohide scroll: reference_min_counts
sc.pp.filter_genes(sc_adata, min_counts=j2_param_destvi_min_counts)
sc_adata.layers["counts"] = sc_adata.X.copy()

save_anndata(sc_adata, "01", "min_counts")
# cell python nohide scroll: reference_min_counts

###############################################################################
# cell markdown null null: reference_variable_genes
"""
### Variable genes

We select the most variable genes.
This selection can create spots with no counts.
So, we need to remove these spots.
"""
# cell markdown null null: reference_variable_genes

###############################################################################
# cell python nohide scroll: reference_variable_genes
sc.pp.highly_variable_genes(
    sc_adata,
    n_top_genes=j2_param_destvi_n_variable_genes,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
)

sc_adata.layers["counts"] = sc_adata.X.copy()
sc_adata = sc_adata[ ~ np.all( sc_adata.X.toarray() == 0 , axis=1 ) , : ]
sc_adata.layers["counts"] = sc_adata.X.copy()

save_anndata(sc_adata, "02", "var_genes")
# cell python nohide scroll: reference_variable_genes

###############################################################################
# cell markdown null null: reference_normalization_log_transform
"""
### Normalization and log-transformation

We normalize and log-transform the counts.
"""
# cell markdown null null: reference_normalization_log_transform

###############################################################################
# cell python nohide scroll: reference_normalization_log_transform
sc.pp.normalize_total(sc_adata, target_sum=10e4)
save_anndata(sc_adata, "03", "normalized")

sc.pp.log1p(sc_adata)
sc_adata.raw = sc_adata
save_anndata(sc_adata, "04", "logtransformed")
# cell python nohide scroll: reference_normalization_log_transform

###############################################################################
# cell markdown null null: title_sample
"""
## The sample
"""
# cell markdown null null: title_sample

###############################################################################
# cell markdown null null: sample_load
"""
### Loading

We load the sample count and spatial data and create an AnnData object.
"""
# cell markdown null null: sample_load

###############################################################################
# cell python nohide scroll: sample_load
spatial = pd.read_csv(j2_path_spatial).set_index("Barcode")

st_adata = sc.read_10x_mtx(j2_path_dge, var_names=j2_param_destvi_gene_identifier)
st_adata.obs = spatial.loc[ st_adata.obs.index ]
st_adata = st_adata[ : , ~ np.all(st_adata.X.toarray() == 0, axis=0) ]
# cell python nohide scroll: sample_load

###############################################################################
# cell markdown null null: sample_preprocessing
"""
### Preprocessing

We need to apply the same preprocessing as with the scRNA-seq reference dataset.
So, we convert the gene symbols to upper case.
Then, we normalize, log-transform the counts and remove the duplicates.
"""
# cell markdown null null: sample_preprocessing

###############################################################################
# cell python nohide scroll: sample_preprocessing

st_adata.var_names = st_adata.var_names.str.upper()

sc.pp.normalize_total(st_adata, target_sum=10e4)
sc.pp.log1p(st_adata)

st_adata.raw = st_adata

sc_adata = sc_adata[ : , ~ sc_adata.var.index.duplicated(keep=False) ]
st_adata = st_adata[ : , ~ st_adata.var.index.duplicated(keep=False) ]
# cell python nohide scroll: sample_preprocessing

###############################################################################
# cell markdown null null: intersection
"""
## Common genes

We select only the genes that are common to the sample and the reference.
The number of genes can end up being quite low because we took only the most variable genes.
"""
# cell markdown null null: intersection

###############################################################################
# cell python nohide scroll: intersection
intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)

print("{0} genes in common over {1}".format(len(intersect), sc_adata.shape[1]))

adata = sc_adata[:, intersect].copy()
adata = adata[ ~ np.all( adata.X.toarray() == 0 , axis=1 ) , : ]
adata.layers["counts"] = adata.X.copy()

st_adata = st_adata[:, intersect].copy()
st_adata = st_adata[ ~ np.all(st_adata.X.toarray() == 0, axis=1) , :]
st_adata.layers["counts"] = st_adata.X.copy()
st_adata.obsm["spatial"] = st_adata.obs[["x", "y"]].to_numpy()

has_zeros("refeference", adata)
has_zeros("spatial", st_adata)
# cell python nohide scroll: intersection

###############################################################################
# cell markdown null null: sc_training_title
"""
## Single-cell model training
"""
# cell markdown null null: sc_training_title

###############################################################################
# cell markdown null null: sc_training_setup
"""
We can now set up the single-cell model.
"""
# cell markdown null null: sc_training_setup

###############################################################################
# cell python nohide scroll: sc_training_setup
adata.obs["Type"] = adata.obs.Type.astype(str) # CondSCVI fails otherwise
CondSCVI.setup_anndata(adata, layer="counts", labels_key="Type")
sc_model = CondSCVI(adata, weight_obs=False)
sc_model.view_anndata_setup()
# cell python nohide scroll: sc_training_setup

###############################################################################
# cell markdown null null: sc_training
"""
We train the model here.
"""
# cell markdown null null: sc_training

# cell python nohide scroll: sc_training
if j2_param_destvi_test: # test mode
	sc_model.train(max_epochs=2)
	sc_model.history["elbo_train"].plot()
else: # normal mode
	sc_model.train(max_epochs=150)
	sc_model.history["elbo_train"].iloc[5:].plot()

plt.show()

save_model(sc_model, "output", "singlecell")
save_anndata(adata, "05", "trained")
# cell python nohide scroll: sc_training

###############################################################################
# cell markdown null null: st_training_title
"""
## Annotation model training
"""
# cell markdown null null: st_training_title

###############################################################################
# cell markdown null null: st_training_setup
"""
We can now set up the single-cell model.
"""
# cell markdown null null: st_training_setup

###############################################################################
# cell python nohide scroll: st_training_setup
DestVI.setup_anndata(st_adata, layer="counts")
st_model = DestVI.from_rna_model(st_adata, sc_model)
st_model.view_anndata_setup()
# cell python nohide scroll: st_training_setup

###############################################################################
# cell markdown null null: st_training
"""
We train the model here.
"""
# cell markdown null null: st_training

###############################################################################
# cell python nohide scroll: st_training

if j2_param_destvi_test: # test mode
	st_model.train(max_epochs=2)
	st_model.history["elbo_train"].plot()
else: # normal mode
	st_model.train(max_epochs=1500)
	st_model.history["elbo_train"].iloc[10:].plot() 

plt.show()

st_adata.obsm["proportions"] = st_model.get_proportions()

save_model(st_model, "output", "spatial")
save_anndata(st_adata, "06", "proportions")
# cell python nohide scroll: st_training

###############################################################################
# cell markdown null null: title_results
"""
## Results
"""
# cell markdown null null: title_results

###############################################################################
# cell markdown null null: max_value
"""
As a first approximation, we just select the max probability for each spot.
"""
# cell markdown null null: max_value

###############################################################################
# cell python nohide scroll: max_value
for cell_type in st_adata.obsm["proportions"].columns:
	data = st_adata.obsm["proportions"][cell_type].values
	st_adata.obs[cell_type] = np.clip(data, 0, np.quantile(data, 0.99))

st_adata.obs["CellType"] = st_adata.obsm["proportions"].idxmax(axis=1)
# cell python nohide scroll: max_value

###############################################################################
# cell markdown null null: cell_types_all
"""
"""
# cell markdown null null: cell_types_all

###############################################################################
# cell python nohide scroll: cell_types_all
sc.settings.figdir = "output"
sc.settings.file_format_figs = "png" 

sc.pl.embedding(
	st_adata,
	basis="spatial",
	color="CellType",
	s=3,
	save=f".cell_type_all",
	show=True
	)
# cell python nohide scroll: cell_types_all

###############################################################################
# cell markdown null null: cell_types_grid
"""
We plot the probability values for each cell type and for each spot.
"""
# cell markdown null null: cell_types_grid

###############################################################################
# cell python nohide scroll: cell_types_grid
sc.settings.figdir = "output"
sc.settings.file_format_figs = "png" 

sc.pl.embedding(
	st_adata,
	basis="spatial",
	color=st_adata.obsm["proportions"].columns,
	cmap="Reds",
	s=3,
	save=f".cell_type_grid",
	show=True
	)
# cell python nohide scroll: cell_types_grid

###############################################################################
# cell markdown null null: session
"""
## Session info
"""
# cell markdown null null: session

###############################################################################
# cell python nohide noscroll: session
session_info.show()
# cell python nohide noscroll: session

