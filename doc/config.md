
# Pipeline configuration

The pipeline configuration is achieved with a [YAML](https://en.wikipedia.org/wiki/YAML) file.
It should look like that roughly:

```yaml
params: params.csv
out_dir: results
project: Nobel Prize
scientist: Paul Nurse
```

The parameters are the following:

 * `params`: the path of the `CSV` parameters file for the samples
 * `out_dir`: the path of the output directory that will contain the results
 * `project`: the name of the project or the experiment
 * `scientist`: the name of the scientist

A config file example can be found as `params.yml` in the repository, or can be downloaded [here](https://bioinformatics.crick.ac.uk/shiny/users/bahn/slideseq/test_data/params.yml).

The `CSV` paramters file contains the parameters of the analysis for each sample.
An example can be found as `params.csv` in the repository, or can be downloaded [here](https://bioinformatics.crick.ac.uk/shiny/users/bahn/slideseq/test_data/params.csv).
It should look roughly like this, with samlples as rows and parameters as columns:

```
$ cat params.csv | cut -d "," -f 1,7-12
name,param_scanpy_min_genes,param_scanpy_max_genes,param_scanpy_max_pct_mitoch,param_scanpy_min_cells,param_scanpy_clusters_resolution,param_seurat_mitochondrial_gene_symbol_prefix
sample1,5,np.infty,30,10,0.2,^mt-
sample2,5,np.infty,30,10,0.2,^mt-
```

There are quite a few parameters:

```
$ head -n 1 params.csv | tr "," "\n"
name
execute
path_dge
path_spatial
path_reference
param_scanpy_mitochondrial_gene_symbol_prefix
param_scanpy_min_genes
param_scanpy_max_genes
param_scanpy_max_pct_mitoch
param_scanpy_min_cells
param_scanpy_clusters_resolution
param_seurat_mitochondrial_gene_symbol_prefix
param_seurat_max_percentage_mitochondria
param_seurat_min_gene_count
param_seurat_sctransform_n_cells
param_seurat_n_components
param_seurat_clusters_resolution
param_seurat_n_cpus
param_seurat_memory
param_destvi_min_counts
param_destvi_n_variable_genes
param_destvi_test
param_rctd_n_cpus
param_rctd_mode
param_sparkx_n_cpus
```

In the following sections we will give more details about the parameters.

## General parameters

 * `name`: name of the sample
 * `execute`: execute the notebooks or not (`true` or `false`)

## File location parameters

 * `path_dge`: the count matrix of the sample (beads versus genes)
 * `path_spatial`: the spatial coordinates of each bead
 * `path_reference`: a reference single dataset (a count matrix with cell types)

A more detailed description of these files can be found [here](file_formats.md).
The paths can full, relative or as an URL.

## Scanpy parameters

These parameters are related to the [Scanpy notebook](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/integration-scanorama.html):

 * `param_scanpy_mitochondrial_gene_symbol_prefix`: mitochondrial gene prefix used for the `qc_vars` parameter of the [`calculate_qc_metrics`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html#scanpy.pp.calculate_qc_metrics) function (usually it's `MT-` or `mt-`)
 * `param_scanpy_min_genes`: threshold to filter beads expressing too few genes
 * `param_scanpy_max_genes`: threshold to filter beads expressing too many genes
 * `param_scanpy_max_pct_mitoch`: threshold to filter beads expressing a high mitochondrial activity
 * `param_scanpy_min_cells`: threshold to filter genes expressed in too few beads
 * `param_scanpy_clusters_resolution`: `resolution` parameter of the [`leiden`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html#scanpy.tl.leiden.) function

## Seurat parameters

These parameters are related to the [Seurat notebook](https://satijalab.org/seurat/articles/spatial_vignette.html#slide-seq-1):

 * `param_seurat_mitochondrial_gene_symbol_prefix`: mitochondrial gene prefix used for the `pattern` parameter when running the [`PercentageFeatureSet`](https://satijalab.org/seurat/reference/percentagefeatureset) function (usually `MT-` or `mt-`)
 * `param_seurat_max_percentage_mitochondria`: threshold to filter beads expressing a high mitochondrial activity
 * `param_seurat_min_gene_count`: threshold to filter beads expressing too few genes
 * `param_seurat_sctransform_n_cells`: `ncells` parameter of the [`SCTranform`](https://satijalab.org/seurat/reference/sctransform) function
 * `param_seurat_n_components`: `dims` parameter for both [`RunUMAP`](https://satijalab.org/seurat/reference/runumap) function and  [`FindNeighbors`](https://satijalab.org/seurat/reference/findneighbors) function
 * `param_seurat_clusters_resolution`: `resolution` parameter for the [`FindClusters`](https://satijalab.org/seurat/reference/findclusters) function
 * `param_seurat_n_cpus`: number of cpus to use when parallelizing `Seurat` with [`BiocParallel`](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html) and [`future`](https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html)
 * `param_seurat_memory`: value used to set the `future.globals.maxSize` parameters of the [`options`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/options.html) function

## DestVI parameters

These parameters are related to the [DestVI notebook](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/DestVI_tutorial.html):

 * `param_destvi_min_counts`: `min_counts` parameter of the [`filter_genes`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.filter_genes.html#scanpy.pp.filter_genes) function
 * `param_destvi_n_variable_genes`: `n_top_genes` of the [`highly_variable_genes`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html) function
 * `param_destvi_test`: running `DestVI` in test mode (`true` or `false`)

## RCTD parameters

These parameters are related to the [RCTD notebook](https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html):

 * `param_rctd_n_cpus`: `max_cores` parameter of the [`create.RCTD`](https://rdrr.io/github/dmcable/RCTD/man/create.RCTD.html) function
 * `param_rctd_mode`: `doublet_mode` parameter of the [`run.RCTD`](https://rdrr.io/github/dmcable/RCTD/man/run.RCTD.html) function

## SPARK-X parameters

These parameters are related to the [SPARK-X notebook](https://xzhoulab.github.io/SPARK/02_SPARK_Example/):

 * `param_sparkx_n_cpus`: `numCores` parameter of the [`sparkx`](https://xzhoulab.github.io/SPARK/02_SPARK_Example/) function

