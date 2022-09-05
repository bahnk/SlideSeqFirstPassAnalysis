#!/usr/bin/Rscript

# j2 variable string: path_dge
# j2 variable string: path_spatial 
# j2 variable integer: param_sparkx_n_cpus

j2_name <- "sample1"
j2_path_dge <- "test/sample1"
j2_path_spatial  <- "test/sample1.csv"
j2_param_sparkx_n_cpus <- 12

###############################################################################
# cell markdown null null: intro
"#markdown
We run [SPARK-X](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02404-0) on the sample.
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

library(SPARK)
library(Seurat)
library(dplyr)
library(future)
library(readr)
library(tibble)
# cell r nohide noscroll: libraries

##############################################################################
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
# cell markdown null null: load_sample
"#markdown
We load the digital expression matrix and the spatial information of the sample.
"#markdown
# cell markdown null null: load_sample

###############################################################################
# cell r nohide noscroll: load_sample

counts <- Seurat::Read10X(j2_path_dge)

spatial <-
    readr::read_csv(j2_path_spatial) %>%
    dplyr::select(SeqBarcode, x, y) %>%
    tibble::column_to_rownames("SeqBarcode")

spatial <- spatial[colnames(counts),]
# cell r nohide noscroll: load_sample

###############################################################################
# cell markdown null null: remove_mitoch
"#markdown
We remove mitochondrial genes if present.
"#markdown
# cell markdown null null: remove_mitoch

###############################################################################
# cell r nohide noscroll: remove_mitoch

mt_idx <- grep("mt-", rownames(counts), ignore.case=T)
if ( length(mt_idx) != 0 )
{
	counts <- counts[-mt_idx,]
}
# cell r nohide noscroll: remove_mitoch

###############################################################################
# cell markdown null null: run_sparkx
"#markdown
We run SPARK-X.
"#markdown
# cell markdown null null: run_sparkx

###############################################################################
# cell r nohide noscroll: run_sparkx

spark <- 
	SPARK::sparkx(
		count_in=counts,
		locus_in=as.matrix(spatial),
		numCores=j2_param_sparkx_n_cpus,
		option="mixture",
		verbose=T
	)
# cell r nohide noscroll: run_sparkx

###############################################################################
# cell markdown null null: export_results
"#markdown
We export the results.
"#markdown
# cell markdown null null: export_results

###############################################################################
# cell r nohide noscroll: export_results

spark$res_mtest %>%
	tibble::rownames_to_column("gene") %>%
	readr::write_csv("output/results.csv")
# cell r nohide noscroll: export_results

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

