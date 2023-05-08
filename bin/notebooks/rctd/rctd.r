#!/usr/bin/Rscript

# j2 variable string: path_reference
# j2 variable string: path_dge
# j2 variable string: path_spatial 
# j2 variable integer: param_rctd_gene_column
# j2 variable integer: param_rctd_n_cpus
# j2 variable string: param_rctd_mode

j2_path_reference <- "results/sample1/reference"
j2_path_dge <- "test/sample1"
j2_path_spatial  <- "test/sample1.csv"
j2_param_rctd_gene_column <- 2
j2_param_rctd_n_cpus <- 12
j2_param_rctd_mode <- "doublet"

###############################################################################
# cell markdown null null: intro
"#markdown
We run [RCTD](https://www.nature.com/articles/s41587-021-00830-w) on the sample.
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

library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(plotly)
library(readr)
library(spacexr)
library(tibble)
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
# cell markdown null null: title_reference
"#markdown
## Create a reference object
"#markdown
# cell markdown null null: title_reference

###############################################################################
# cell markdown null null: load_reference
"#markdown
We load the digital expression matrix and the cell types for the reference.
"#markdown
# cell markdown null null: load_reference

###############################################################################
# cell r nohide noscroll: load_reference

ref_counts <- Seurat::Read10X(j2_path_reference, gene.column=j2_param_rctd_gene_column)

cell_types <-
	readr::read_tsv(file.path(j2_path_reference, "types.tsv.gz"), col_names=F) %>%
	dplyr::pull(X1) %>%
	setNames(colnames(ref_counts)) %>%
	as.factor()

ref_nUMIs <- colSums(ref_counts)
# cell r nohide noscroll: load_reference

###############################################################################
# cell markdown null null: create_reference
"#markdown
We create a RCTD reference object.
"#markdown
# cell markdown null null: create_reference

###############################################################################
# cell r nohide noscroll: reference_reference
ref <- spacexr::Reference(round(ref_counts), cell_types, ref_nUMIs)
saveRDS(ref, "output/reference.rds")
# cell r nohide noscroll: reference_reference

###############################################################################
# cell markdown null null: title_sample
"#markdown
## Create a sample object
"#markdown
# cell markdown null null: title_sample

###############################################################################
# cell markdown null null: load_sample
"#markdown
We load the digital expression matrix and the spatial information of the sample.
"#markdown
# cell markdown null null: load_sample

###############################################################################
# cell r nohide noscroll: load_sample

counts <- Seurat::Read10X(j2_path_dge, gene.column=j2_param_rctd_gene_column)

nUMIs <- colSums(counts)

spatial <-
    readr::read_csv(j2_path_spatial) %>%
    tibble::column_to_rownames("Barcode")

spatial <- spatial[colnames(counts),]

spots <- spacexr::SpatialRNA(spatial, counts, nUMIs)

# cell r nohide noscroll: load_sample

###############################################################################
# cell markdown null null: create_sample
"#markdown
We create a RCTD object.
"#markdown
# cell markdown null null: create_sample

###############################################################################
# cell r nohide noscroll: create_sample
plan("multicore", workers=j2_param_rctd_n_cpus)
rctd <- spacexr::create.RCTD(spots, ref, max_cores=j2_param_rctd_n_cpus)
# cell r nohide noscroll: create_sample

###############################################################################
# cell markdown null null: title_rctd
"#markdown
## Run RCTD
"#markdown
# cell markdown null null: title_rctd

###############################################################################
# cell r nohide noscroll: run_rctd
rctd <- spacexr::run.RCTD(rctd, doublet_mode=j2_param_rctd_mode)

saveRDS(rctd, "output/rctd.rds")

readr::write_csv(
	rctd@results$results_df %>% tibble::rownames_to_column("cell"),
	"output/rctd.csv"
)
# cell r nohide noscroll: run_rctd

###############################################################################
# cell markdown null null: title_results
"#markdown
## Results
"#markdown
# cell markdown null null: title_results

###############################################################################
# cell markdown null null: results_first_type
"#markdown
### First type
"#markdown
# cell markdown null null: results_first_type

###############################################################################
# cell r nohide noscroll: plot_first_type

if ( "first_type" %in% colnames(rctd@results$results_df) )
{
	df <-
		rctd@results$results_df %>%
		tibble::rownames_to_column("cell") %>%
		dplyr::inner_join( spatial %>% tibble::rownames_to_column("cell") )
	
	fig <- plotly::plot_ly(data=df, x=~x, y=~y, color=~first_type)
	
	plotly::save_image(fig, "output/first_type.png")
	plotly::save_image(fig, "output/first_type.pdf")
	
	#fig

	g <- ggplot2::ggplot(df, aes(x=x, y=y, color=first_type)) + geom_point()
	g
}
# cell r nohide noscroll: plot_first_type

###############################################################################
# cell markdown null null: results_second_type
"#markdown
### Second type
"#markdown
# cell markdown null null: results_second_type

###############################################################################
# cell r nohide noscroll: plot_second_type

if ( "second_type" %in% colnames(rctd@results$results_df) )
{
	df <-
		rctd@results$results_df %>%
		tibble::rownames_to_column("cell") %>%
		dplyr::inner_join( spatial %>% tibble::rownames_to_column("cell") )
	
	fig <- plotly::plot_ly(data=df, x=~x, y=~y, color=~second_type)
	
	plotly::save_image(fig, "output/second_type.png")
	plotly::save_image(fig, "output/second_type.pdf")
	
	#fig

	g <- ggplot2::ggplot(df, aes(x=x, y=y, color=second_type)) + geom_point()
	g
}
# cell r nohide noscroll: plot_second_type

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

