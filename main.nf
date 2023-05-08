#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.nio.file.Paths
import nextflow.util.ArrayBag

///////////////////////////////////////////////////////////////////////////////
//// METHODS //////////////////////////////////////////////////////////////////

include {
	addValue;
	removeKeys;
	toBoolean;
	getDataPaths;
	checkAllParams;
	addMetrics
	} from "./modules/utils.nf"

///////////////////////////////////////////////////////////////////////////////
//// PROCESSES ////////////////////////////////////////////////////////////////

////////////
// nbconvert
nbconvert_template = Channel.fromPath("$workflow.projectDir/assets/jupyter-nbconvert/template.tpl")

/////////////
// parameters
include { params_json } from "./modules/process/params.nf"

/////////
// seurat
include { notebook_r as seurat } from "./modules/process/notebooks.nf"
render_seurat = Channel.fromPath("$workflow.projectDir/bin/notebooks/seurat/render.py")
j2_seurat = Channel.fromPath("$workflow.projectDir/bin/notebooks/seurat/j2")
prefix_seurat = Channel.value("seurat")

///////
// rctd
include { notebook_r as rctd } from "./modules/process/notebooks.nf"
render_rctd = Channel.fromPath("$workflow.projectDir/bin/notebooks/rctd/render.py")
j2_rctd = Channel.fromPath("$workflow.projectDir/bin/notebooks/rctd/j2")
prefix_rctd = Channel.value("rctd")

/////////
// sparkx
include { notebook_r as sparkx } from "./modules/process/notebooks.nf"
render_sparkx = Channel.fromPath("$workflow.projectDir/bin/notebooks/sparkx/render.py")
j2_sparkx = Channel.fromPath("$workflow.projectDir/bin/notebooks/sparkx/j2")
prefix_sparkx = Channel.value("sparkx")

/////////
// scanpy
include { notebook_py as scanpy } from "./modules/process/notebooks.nf"
render_scanpy = Channel.fromPath("$workflow.projectDir/bin/notebooks/scanpy/render.py")
j2_scanpy = Channel.fromPath("$workflow.projectDir/bin/notebooks/scanpy/j2")
prefix_scanpy = Channel.value("scanpy")

/////////
// destvi
include { notebook_py_gpu as destvi } from "./modules/process/notebooks.nf"
render_destvi = Channel.fromPath("$workflow.projectDir/bin/notebooks/destvi/render.py")
j2_destvi = Channel.fromPath("$workflow.projectDir/bin/notebooks/destvi/j2")
prefix_destvi = Channel.value("destvi")

///////////////
// jupyter_book
include { jupyter_book } from "./modules/process/report.nf"
jupyter_book_conf = Channel.fromPath("$workflow.projectDir/assets/jupyter-book/conf.py")
jupyter_book_logo = Channel.fromPath("$workflow.projectDir/assets/jupyter-book/logo.jpg")
jupyter_book_index = Channel.fromPath("$workflow.projectDir/assets/jupyter-book/index.md")

//////////
// metrics
include { metrics } from "./modules/local/metrics/main"
bin_metrics = Channel.fromPath("$workflow.projectDir/bin/metrics.py")

///////////////////////////////////////////////////////////////////////////////
//// DESIGN ///////////////////////////////////////////////////////////////////

scanpy_meta = "$workflow.projectDir/modules/local/scanpy/meta.yml"
seurat_meta = "$workflow.projectDir/modules/local/seurat/meta.yml"
rctd_meta = "$workflow.projectDir/modules/local/rctd/meta.yml"
destvi_meta = "$workflow.projectDir/modules/local/destvi/meta.yml"
sparkx_meta = "$workflow.projectDir/modules/local/sparkx/meta.yml"
metas = [scanpy_meta, seurat_meta, rctd_meta, destvi_meta, sparkx_meta]

Channel
	.fromPath(params.params)
	.splitCsv(header: true)
	.map{ toBoolean(it) }
	.map{ getDataPaths(it) }
	.map{ [ addValue(it[0], "project", params.project) , it[1] ] }
	.map{ [ addValue(it[0], "scientist", params.scientist) , it[1] ] }
	.map{ [
		[
			"name": it[0].name,
			"execute": it[0].execute,
			"execution": [:],
			"parameters": removeKeys(it[0], ["execute"])
		],
		it[1],
		metas
	] }
	.map{ checkAllParams(it) }
	.set{SAMPLES}

///////////////////////////////////////////////////////////////////////////////
//// MAIN WORKFLOW ////////////////////////////////////////////////////////////

workflow {

	params_json(SAMPLES)

	metrics(params_json.out.combine(bin_metrics))

	// filters
	metrics
		.out
		.map{[ addMetrics(it[0], it[3].toString()) , *it[1..2] ]}
		.filter{ it[0]["metrics"]["Total beads"] >= params.min_beads }
		.set{ SELECTED_SAMPLES }

	seurat(
		SELECTED_SAMPLES
			.filter{ it[0]["execution"]["seurat"] }
			.combine(prefix_seurat)
			.combine(render_seurat)
			.combine(j2_seurat)
			.combine(nbconvert_template)
	)

	rctd(
		SELECTED_SAMPLES
			.filter{ it[0]["execution"]["rctd"] }
			.combine(prefix_rctd)
			.combine(render_rctd)
			.combine(j2_rctd)
			.combine(nbconvert_template)
	)

	sparkx(
		SELECTED_SAMPLES
			.filter{ it[0]["execution"]["sparkx"] }
			.combine(prefix_sparkx)
			.combine(render_sparkx)
			.combine(j2_sparkx)
			.combine(nbconvert_template)
	)

	scanpy(
		SELECTED_SAMPLES
			.filter{ it[0]["execution"]["scanpy"] }
			.combine(prefix_scanpy)
			.combine(render_scanpy)
			.combine(j2_scanpy)
			.combine(nbconvert_template)
	)

	destvi(
		SELECTED_SAMPLES
			.filter{ it[0]["execution"]["destvi"] }
			.combine(prefix_destvi)
			.combine(render_destvi)
			.combine(j2_destvi)
			.combine(nbconvert_template)
	)

	seurat
		.out.ipynb
		.concat(rctd.out.ipynb)
		.concat(sparkx.out.ipynb)
		.concat(scanpy.out.ipynb)
		.concat(destvi.out.ipynb)
		.collect{ it[3] }
		.map{[
			["project": params.project, "scientist": params.scientist],
			it
		]}
		.combine(jupyter_book_conf)
		.combine(jupyter_book_index)
		.combine(jupyter_book_logo)
		.set{TO_JUPYTER_BOOK}
	
	jupyter_book(TO_JUPYTER_BOOK)
}

