#!/usr/bin/Rscript

# j2 variable string: name
# j2 variable string: path_dge
# j2 variable string: path_spatial 
# j2 variable integer: param_seurat_gene_column
# j2 variable string: param_seurat_mitochondrial_gene_symbol_prefix
# j2 variable float: param_seurat_max_percentage_mitochondria
# j2 variable integer: param_seurat_min_gene_count
# j2 variable integer: param_seurat_sctransform_n_cells
# j2 variable integer: param_seurat_n_components
# j2 variable float: param_seurat_clusters_resolution
# j2 variable integer: param_seurat_n_cpus
# j2 variable integer: param_seurat_memory

j2_name <- "210709_26"
j2_path_dge <- "../results/210709_26_dge"
j2_path_spatial  <- "../results/210709_26.csv"
j2_param_seurat_gene_column <- 2
j2_param_seurat_mitochondrial_gene_symbol_prefix <- "^mt-"
j2_param_seurat_max_percentage_mitochondria <- 30
j2_param_seurat_min_gene_count <- 10
j2_param_seurat_sctransform_n_cells <- 3000
j2_param_seurat_n_components <- 10
j2_param_seurat_clusters_resolution <- 0.2
j2_param_seurat_n_cpus <- 12
j2_param_seurat_memory <- 5000

###############################################################################
# cell markdown null null: intro
"#markdown
This is a very a very basic analysis of this sample with Seurat.
"#markdown
# cell markdown null null: intro

###############################################################################
# cell markdown null null: libraries
"#markdown
Here are the libraries we need.
"#markdown
# cell markdown null null: libraries

###############################################################################
# cell r nohide noscroll: libraries

library(BiocParallel)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(future)
library(ggplot2)
library(patchwork)
library(purrr)
library(sceasy)
# cell r nohide noscroll: libraries

###############################################################################
# cell markdown null null: directories
"#markdown
We create the output directory for this noteboook.
Every outputs will save there.
"#markdown
# cell markdown null null: directories

###############################################################################
# cell r nohide noscroll: directories
dir.create("output", showWarnings=F)
# cell r nohide noscroll: directories

###############################################################################
# cell markdown null null: load
"#markdown
We load the digital expression matrix and the spatial information and create an AnnData object.
"#markdown
# cell markdown null null: load

###############################################################################
# cell r nohide noscroll: load

counts <- Seurat::Read10X(j2_path_dge, gene.column=j2_param_seurat_gene_column)

spatial <-
	readr::read_csv(j2_path_spatial) %>%
	dplyr::filter( Barcode %in% colnames(counts) ) %>%
	dplyr::arrange(Barcode) %>%
	tibble::column_to_rownames("Barcode")

stopifnot( sum( ! spatial$Barcode == colnames(counts) ) == 0 )

args <- list("project"=j2_name, "assay"="Spatial", "meta.data"=spatial)
obj <- do.call(SeuratObject::CreateSeuratObject, c(counts, args))

img <- new(
	"SlideSeq",
	coordinates=obj@meta.data[,c("x","y")],
	assay="Spatial",
	key="_images"
	)

obj@images <- list(image=img)

obj$log_nCount_Spatial <- log(obj$nCount_Spatial)

# cell r nohide noscroll: load

###############################################################################
# cell markdown null null: qc_title
"#markdown
## Quality control
"#markdown
# cell markdown null null: qc_title

###############################################################################
# cell markdown null null: gene_counts
"#markdown
### Gene counts

The distribution of the gene counts.
"#markdown
# cell markdown null null: gene_counts

###############################################################################
# cell r nohide noscroll: violin_plot_count

obj %>%
	Seurat::VlnPlot(., features="nCount_Spatial", pt.size=0, log=TRUE) +
	Seurat::NoLegend() -> vln_count_plot

ggplot2::ggsave("output/violin_count.png", vln_count_plot, width=9, height=6)
ggplot2::ggsave("output/violin_count.pdf", vln_count_plot, width=9, height=6)

vln_count_plot

# cell r nohide noscroll: violin_plot_count

###############################################################################
# cell markdown null null: spatial_umis
"#markdown
### UMIs

UMIs per bead.
"#markdown
# cell markdown null null: spatial_umis

###############################################################################
# cell r nohide noscroll: spatial_umi_plot

obj %>%
	Seurat::SpatialFeaturePlot(., features="log_nCount_Spatial") +
	ggplot2::theme(legend.position="right") -> spatial_umi_plot

ggplot2::ggsave("output/spatial_umi.png", spatial_umi_plot, width=9, height=6)
ggplot2::ggsave("output/spatial_umi.pdf", spatial_umi_plot, width=9, height=6)

spatial_umi_plot
# cell r nohide noscroll: spatial_umi_plot

###############################################################################
# cell markdown null null: violin_mitoch
"#markdown
### Percentage of mitochondria

Distribution of percentage of detected mitochondrial genes.
"#markdown
# cell markdown null null: violin_mitoch

###############################################################################
# cell r nohide noscroll: violin_plot_mitoch

obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern=j2_param_seurat_mitochondrial_gene_symbol_prefix)

obj %>%
	Seurat::VlnPlot(., features="percent.mt", pt.size=0, log=TRUE) +
	Seurat::NoLegend() -> vln_mitoch_plot

ggplot2::ggsave("output/violin_mitoch.png", vln_mitoch_plot, width=9, height=6)
ggplot2::ggsave("output/violin_mitoch.pdf", vln_mitoch_plot, width=9, height=6)

vln_mitoch_plot

# cell r nohide noscroll: violin_plot_mitoch

###############################################################################
# cell markdown null null: scatter_count_mitoch
"#markdown
### Abnormal beads

Gene counts vers mitochondrial percentage
"#markdown
# cell markdown null null: scatter_count_mitoch

###############################################################################
# cell r nohide noscroll: scatter_plot_count_mitoch

scatter_plot <-
	Seurat::FeatureScatter(
		obj,
		feature1="nCount_Spatial",
		feature2="percent.mt"
	)

ggplot2::ggsave("output/scatter_count_mitoch.png", scatter_plot, width=9, height=6)
ggplot2::ggsave("output/scatter_count_mitoch.pdf", scatter_plot, width=9, height=6)

scatter_plot

# cell r nohide noscroll: scatter_plot_count_mitoch

###############################################################################
# cell markdown null null: filter
"#markdown
## Filtering

We remove abnormal cells.
"#markdown
# cell markdown null null: filter

################################################################################
# cell r nohide noscroll: filter

cat("Original\n")
print(obj)

cat("Filter mitochondria\n")
obj <- base::subset(obj, percent.mt < j2_param_seurat_max_percentage_mitochondria)
print(obj)

cat("Filter UMIs\n")
obj <- base::subset(obj, nCount_Spatial >= j2_param_seurat_min_gene_count)
print(obj)

# cell r nohide noscroll: filter

###############################################################################
# cell markdown null null: analysis_title
"#markdown
## Clustering
"#markdown
# cell markdown null null: analysis_title

###############################################################################
# cell r nohide scroll: sctransform
obj <- Seurat::SCTransform(obj, assay="Spatial", ncells=j2_param_seurat_sctransform_n_cells, verbose=T)
# cell r nohide scroll: sctransform

###############################################################################
# cell r nohide scroll: pca
obj <- Seurat::RunPCA(obj)
# cell r nohide scroll: pca

###############################################################################
# cell r nohide scroll: umap
obj <- Seurat::RunUMAP(obj, dims=1:j2_param_seurat_n_components)
# cell r nohide scroll: umap

###############################################################################
# cell r nohide noscroll: neighbors
obj <- Seurat::FindNeighbors(obj, dims=1:j2_param_seurat_n_components)
# cell r nohide noscroll: neighbors

###############################################################################
# cell r nohide noscroll: clusters
BiocParallel::register(MulticoreParam(j2_param_seurat_n_cpus))
base::options(future.globals.maxSize = j2_param_seurat_memory * 1024^2)
future::plan("multiprocess", workers=j2_param_seurat_n_cpus)

obj <- Seurat::FindClusters(obj, resolution=j2_param_seurat_clusters_resolution, verbose=T)
# cell r nohide noscroll: clusters

###############################################################################
# cell r nohide noscroll: plot_umap_clusters
umap_plot <- Seurat::DimPlot(obj, reduction="umap", label=TRUE)

ggplot2::ggsave("output/umap.png", umap_plot, width=9, height=6)
ggplot2::ggsave("output/umap.pdf", umap_plot, width=9, height=6)

umap_plot
# cell r nohide noscroll: plot_umap_clusters

###############################################################################
# cell r nohide noscroll: plot_spatial_clusters
spatial_cluster_plot <- Seurat::SpatialDimPlot(obj, stroke=0)

ggplot2::ggsave("output/spatial.png", spatial_cluster_plot, width=9, height=6)
ggplot2::ggsave("output/spatial.pdf", spatial_cluster_plot, width=9, height=6)

spatial_cluster_plot
# cell r nohide noscroll: plot_spatial_clusters

###############################################################################
# cell r nohide noscroll: plot_spatial_individual_clusters

obj %>%
	slot(., "meta.data") %>%
	dplyr::pull(paste0("SCT_snn_res.", j2_param_seurat_clusters_resolution)) %>%
	levels() %>%
	as.numeric(.) %>%
	sort() -> clusters


ncol <- 3
1:ceiling( length(clusters) / ncol) %>%
	purrr::map(function(x) ncol * (x-1) + 1:ncol) %>%
	purrr::map(function(x) clusters[x]) %>%
	purrr::map(function(x) {

		plot <- Seurat::SpatialDimPlot(
			obj,
			cells.highlight=Seurat::CellsByIdentities(object=obj, idents=x),
			facet.highlight=T
		)

		basename <- paste0("output/individual_cluster_", paste(na.omit(x), collapse="-"))

		ggplot2::ggsave(
			filename=paste0(basename, ".png"),
			plot=plot,
			widt=9, height=6, dpi=300
			)

		ggplot2::ggsave(
			filename=paste0(basename, ".pdf"),
			plot=plot,
			widt=9, height=6, dpi=300
			)

		print(plot)

		list("clusters"=x, "plot"=plot)

	}) -> clusters_plots

# cell r nohide noscroll: plot_spatial_individual_clusters

###############################################################################
# cell markdown null null: markers_title
"#markdown
## Marker genes
"#markdown
# cell markdown null null: markers_title

###############################################################################
# cell r nohide noscroll: find_all_markers
markers <- Seurat::FindAllMarkers(obj, assay="SCT", only.pos=T)
write.csv(markers, "output/markers.csv", row.names=F)
markers
# cell r nohide noscroll: find_all_markers

# cell r nohide noscroll: plot_all_markers

if ( nrow(markers) > 0 )
{
	markers %>%
		dplyr::group_by(cluster) %>%
		dplyr::group_map(function(x, n) {
	
			df <-
				x %>%
				dplyr::arrange(dplyr::desc(avg_log2FC) ) %>%
				dplyr::top_n( max(6, nrow(x)) )
	
			genes <- df$gene
	
			for ( gene in genes ) {
	
				g <- Seurat::FeaturePlot(obj, features=gene, ncol=1)
				g <- g + patchwork::plot_annotation(title=paste0("Cluster ", n))
	
				ggplot2::ggsave(
					filename=sprintf("output/markers_cluster_%s_gene_%s.png", n, gene),
					plot=g,
					height=3, width=3.5, dpi=200
					)
	
				ggplot2::ggsave(
					filename=sprintf("output/markers_cluster_%s_gene_%s.pdf", n, gene),
					plot=g,
					height=3, width=3.5, dpi=200
					)
	
				print(g)
			}
		})
}

# cell r nohide noscroll: plot_all_markers

###############################################################################
# cell markdown null null: save
"#markdown
## Save
"#markdown
# cell markdown null null: save

###############################################################################
# cell r nohide noscroll: save
saveRDS(obj, "output/seurat_object.rds")

sceasy::convertFormat(
	obj,
	from="seurat",
	to="anndata",
	assay="Spatial",
	outFile="output/seurat_object.h5ad"
	)
# cell r nohide noscroll: save

###############################################################################
# cell markdown null null: session
"#markdown
## Session info
"#markdown
# cell markdown null null: session

###############################################################################
# cell r nohide noscroll: session
sessionInfo()
# cell r nohide noscroll: session


################################################################################
#cat("Variable\n")
#
#Seurat::DefaultAssay(obj) <- "SCT"
#
#obj <- Seurat::FindSpatiallyVariableFeatures(
#	obj,
#	assay="SCT",
#	slot="scale.data",
#	features=Seurat::VariableFeatures(obj)[1:1000],
#	selection.method="moransi",
#	x.cuts=100,
#	y.cuts=100
#	)
#
#features <- Seurat::SpatiallyVariableFeatures(obj, selection.method="moransi")
#
#spatial_feat_plot <- SpatialFeaturePlot(
#	obj,
#	features=head(features, 6),
#	ncol=3,
#	alpha=c(0.1, 1),
#	max.cutoff="q95"
#	)
#
#basepath <- file.path(path_outdir, sprintf("%s_plot_spatial_feat", name))
#ggplot2::ggsave(filename=sprintf("%s.png", basepath), plot=spatial_feat_plot)
#ggplot2::ggsave(filename=sprintf("%s.pdf", basepath), plot=spatial_feat_plot)

