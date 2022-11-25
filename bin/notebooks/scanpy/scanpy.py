#!/usr/bin/env python

# coding: utf-8

# j2 variable string: name
# j2 variable string: path_dge
# j2 variable string: path_spatial 
# j2 variable string: param_scanpy_gene_identifier
# j2 variable string: param_scanpy_mitochondrial_gene_symbol_prefix
# j2 variable value: param_scanpy_min_genes
# j2 variable value: param_scanpy_max_genes
# j2 variable value: param_scanpy_max_pct_mitoch
# j2 variable value: param_scanpy_min_cells
# j2 variable value: param_scanpy_clusters_resolution

j2_name = "sample1"
j2_path_dge = "test/sample1"
j2_path_spatial = "test/sample1.csv"
j2_param_scanpy_gene_identifier = "gene_symbols"
j2_param_scanpy_mitochondrial_gene_symbol_prefix = "mt-"
j2_param_scanpy_min_genes = 5
j2_param_scanpy_max_genes = np.infty
j2_param_scanpy_max_pct_mitoch = 30
j2_param_scanpy_min_cells = 10
j2_param_scanpy_clusters_resolution = 0.2

###############################################################################
# cell markdown null null: intro
"""
This is a very a very basic analysis of this sample with scanpy.
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

from IPython.core.display import display, Image
from numpy.random import choice
from os import makedirs
from os.path import join
from plotly.subplots import make_subplots

import ipywidgets as widgets
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
import seaborn as sns
import session_info
# cell python nohide scroll: libraries

###############################################################################
# cell markdown null null: functions
"""
Some useful functions.
"""
# cell markdown null null: functions

###############################################################################
# cell python nohide scroll: functions

##################
def in_jupyter():#
##################
	try:
		__IPYTHON__
	except NameError:
		return False
	else:
		return "IPKernelApp" in get_ipython().config.keys()
	############################################################################

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
	makedirs("output")
except OSError as e:
	pass


# cell python nohide scroll: directories

###############################################################################
# cell markdown null null: config
"""
Configuration of scanpy.
"""
# cell markdown null null: config

###############################################################################
# cell python nohide scroll: config

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3
# cell python nohide scroll: config

###############################################################################
# cell markdown null null: load
"""
We load the digital expression matrix and the spatial information and create an AnnData object.
"""
# cell markdown null null: load

###############################################################################
# cell python nohide scroll: load

spatial = pd.read_csv(j2_path_spatial).set_index("Barcode")

adata = sc.read_10x_mtx(j2_path_dge, var_names=j2_param_scanpy_gene_identifier)
adata.obs = spatial.loc[ adata.obs.index ]
# cell python nohide scroll: load

###############################################################################
# cell markdown null null: qc_metrics
"""
## Quality control

Scanpy will compute basic QC metrics for us.
"""
# cell markdown null null: qc_metrics

###############################################################################
# cell python nohide scroll: qc_metrics

adata.var["mt"] = adata.var_names.str.match(j2_param_scanpy_mitochondrial_gene_symbol_prefix)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
# cell python nohide scroll: qc_metrics


###############################################################################
# cell markdown null null: qc_plots_gene_counts
"""
### Gene counts

Genes detected for each bead.
"""
# cell markdown null null: qc_plots_gene_counts

###############################################################################
# cell python nohide noscroll: qc_plots_gene_counts

plots = [
	("total_counts", "Total counts"),
	("n_genes_by_counts", "Detected genes")
]

tab_content = []
tab_titles = []

n = min(adata.obs.shape[0], 3000)

for column, title in plots:

	fig = px.violin(adata.obs.sample(n), y=column, box=True, points="all")
	fig.write_image( join("output", f"qc_beads.{column}.pdf") )
	fig.write_image( join("output", f"qc_beads.{column}.png") )
	tab_content.append( go.FigureWidget(fig) )
	tab_titles.append(title)

	fig = px.scatter(adata.obs, x="x", y="y", color=column)
	fig.write_image( join("output", f"qc_beads.{column}.spatial.pdf") )
	fig.write_image( join("output", f"qc_beads.{column}.spatial.png") )
	tab_content.append( go.FigureWidget(fig) )
	tab_titles.append(f"{title} (xy)")

fig = px.scatter(adata.obs.sample(n), x="total_counts", y="n_genes_by_counts")
fig.write_image( join("output", "total_counts_VS_n_genes_by_counts.png") )
fig.write_image( join("output", "total_counts_VS_n_genes_by_counts.pdf") )

tab_content += [go.FigureWidget(fig)]
tab_titles += ["Counts VS Genes"]

# display plots
if in_jupyter():
	tab = widgets.Tab(tab_content)
	for i, title in enumerate(tab_titles):
		tab.set_title(i, title)
	display(tab)

# cell python nohide noscroll: qc_plots_gene_counts

###############################################################################
# cell markdown null null: qc_plots_genes
"""
### Genes

The dropouts.
"""
# cell markdown null null: qc_plots_genes

###############################################################################
# cell python nohide noscroll: qc_plots_genes

adata.var["detected"] = np.where(adata.var.n_cells_by_counts > 0,
	"Detected", "Undetected")

df = adata\
	.var\
	.detected\
	.value_counts()\
	.to_frame()\
	.reset_index()\
	.rename(columns={"index": "Gene", "detected": "Number"})

fig = px.bar(df, x="Gene", y="Number", text_auto=True)

tab_content = [go.FigureWidget(fig)]
tab_titles = ["Detected genes"]

plots = [
	("total_counts", "Total counts"),
	("mean_counts", "Mean counts"),
	("n_cells_by_counts", "Beads"),
	("pct_dropout_by_counts", "Dropout")
]

n = min(adata.obs.shape[0], 3000)

for column, title in plots:
	fig = px.violin(adata.var.sample(n), y=column, box=True, points="all")
	fig.write_image( join("output", f"qc_genes.{column}.pdf") )
	fig.write_image( join("output", f"qc_genes.{column}.png") )
	tab_content.append( go.FigureWidget(fig) )
	tab_titles.append(title)

# display plots
if in_jupyter():
	tab = widgets.Tab(tab_content)
	for i, title in enumerate(tab_titles):
		tab.set_title(i, title)
	display(tab)

# cell python nohide noscroll: qc_plots_genes

###############################################################################
# cell markdown null null: qc_plots_mitoch
"""
### Mitochondrial activity

Percentage of mitochondrial genes detected for each bead.
"""
# cell markdown null null: qc_plots_mitoch

###############################################################################
# cell python nohide noscroll: qc_plots_mitoch
	#("pct_counts_mt", "Percent Mitoch")

n = min(adata.obs.shape[0], 3000)

fig1 = px.violin(adata.obs.sample(n), y="pct_counts_mt", box=True, points="all")
fig1.write_image( join("output", "pct_counts_mt.pdf") )
fig1.write_image( join("output", "pct_counts_mt.png") )

fig2 = px.scatter(adata.obs.sample(n), x="total_counts", y="pct_counts_mt")
fig2.write_image( join("output", "total_counts_VS_pct_counts_mt.png") )
fig2.write_image( join("output", "total_counts_VS_pct_counts_mt.pdf") )

# display plots
if in_jupyter():
	tab = widgets.Tab([go.FigureWidget(fig1), go.FigureWidget(fig2)])
	tab.set_title(0, "Percent Mitoch")
	tab.set_title(1, "Counts VS Mitoch")
	display(tab)

# cell python nohide noscroll: qc_plots_mitoch

###############################################################################
# cell markdown null null: filtering
"""
## Filtering

We filter anormal beads.
"""
# cell markdown null null: filtering

###############################################################################
# cell python nohide noscroll: filtering

print("Before filtering:")
print(adata)

adata = adata[ adata.obs.n_genes_by_counts > j2_param_scanpy_min_genes , :]
adata = adata[ adata.obs.n_genes_by_counts < j2_param_scanpy_max_genes , :]
adata = adata[ adata.obs.pct_counts_mt < j2_param_scanpy_max_pct_mitoch , :]

sc.pp.filter_genes(adata, min_cells=j2_param_scanpy_min_cells)

print("\n\n")
print("After filtering:")
print(adata)
# cell python nohide noscroll: filtering

###############################################################################
# cell markdown null null: transform
"""
## Data transformation

Normalization by sequencing depth to remove technical variability.
And nonlinear transformation to stabilize the variance across genes with different expression levels.
"""
# cell markdown null null: transform

###############################################################################
# cell python nohide noscroll: transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata)
# cell python nohide noscroll: transform

###############################################################################
# cell markdown null null: pca
"""
## Principal components analysis

We compute the PCA.
"""
# cell markdown null null: pca

###############################################################################
# cell python nohide noscroll: pca
sc.pp.pca(adata, use_highly_variable=False)
# cell python nohide noscroll: pca

###############################################################################
# cell markdown null null: plot_pca_variance
"""
We plot the variance explained by the different principal components.
"""
# cell markdown null null: plot_pca_variance

###############################################################################
# cell python nohide noscroll: plot_pca_variance
var_ratio = adata.uns["pca"]["variance_ratio"]

fig = px.scatter(x=range(1, len(var_ratio)+1), y=var_ratio)
fig.write_image( join("output", "pca_variance_ratio.png") )
fig.write_image( join("output", "pca_variance_ratio.pdf") )

fig
# cell python nohide noscroll: plot_pca_variance

###############################################################################
# cell markdown null null: plot_pca
"""
We plot the 5 first principal components.
"""
# cell markdown null null: plot_pca

###############################################################################
# cell python nohide noscroll: plot_pca
df = pd.DataFrame(
	np.hstack([adata.obs.index[:,None], adata.obsm["X_pca"][:,:5]]),
	columns=["Bead"] + [f"PC{i+1}" for i in range(5)]
	)\
	.melt(id_vars=["Bead", "PC1"], var_name="PC", value_name="Value")

fig = px.scatter(df, x="PC1", y="Value", facet_col="PC", facet_col_wrap=2)
fig.write_image( join("output", "pca.png") )
fig.write_image( join("output", "pca.pdf") )

fig
# cell python nohide noscroll: plot_pca

###############################################################################
# cell markdown null null: clustering
"""
## Clustering & UMAP

We cluster the beads based on their expression.
"""
# cell markdown null null: clustering

###############################################################################
# cell python nohide noscroll: clustering
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10, use_rep="X_pca")
sc.tl.leiden(adata, key_added="clusters", resolution=j2_param_scanpy_clusters_resolution)
# cell python nohide noscroll: clustering

###############################################################################
# cell markdown null null: umap
"""
We compute the UMAP.
"""
# cell markdown null null: umap

###############################################################################
# cell python nohide noscroll: umap
sc.tl.umap(adata)
# cell python nohide noscroll: umap

###############################################################################
# cell markdown null null: plot_clusters_umap
"""
We plot the clusters on top of the UMAP.
"""
# cell markdown null null: plot_clusters_umap

###############################################################################
# cell python nohide noscroll: plot_clusters_umap
args = {
		"x": np.apply_along_axis(lambda x: x[0], arr=adata.obsm["X_umap"], axis=1),
		"y": np.apply_along_axis(lambda x: x[1], arr=adata.obsm["X_umap"], axis=1),
		"color": adata.obs.clusters
}
fig = px.scatter(**args)
fig.write_image( join("output", "clusters_umap.png") )
fig
# cell python nohide noscroll: plot_clusters_umap

###############################################################################
# cell markdown null null: plot_clusters_spatial
"""
We plot the clusters on top of the spatial coordinates.
"""
# cell markdown null null: plot_clusters_spatial

###############################################################################
# cell python nohide noscroll: plot_clusters_spatial
fig = px.scatter(adata.obs, x="x", y="y", color="clusters")
fig.write_image( join("output", "clusters_spatial.png") )
fig
# cell python nohide noscroll: plot_clusters_spatial

###############################################################################
# cell markdown null null: save
"""
## Save
"""
# cell markdown null null: save

###############################################################################
# cell python nohide noscroll: save
adata.write_h5ad( join("output", "anndata.h5ad") )
# cell python nohide noscroll: save

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

