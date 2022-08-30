#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.nio.file.Paths

///////////////////////////////////////////////////////////////////////////////
//// METHODS //////////////////////////////////////////////////////////////////

/////////////////////////////////
def addValue(map, key, value) {//
/////////////////////////////////

	def new_map = map.clone()

	new_map.put(key, value)

	return new_map
}

/////////////////////////////
def removeKeys(map, keys) {//
/////////////////////////////

	def new_map = [:]

	map.each{
		if ( ! keys.contains(it.key) )
		{
			new_map.put(it.key, it.value)
		}
	}

	return new_map
}

/////////////////////
def getPaths(map) {//
/////////////////////

	def new_map = [:]
	def files = []

	map.each{
		if ( it.key.startsWith("path_") )
		{
			if ( it.key == "path_reference" )
			{
				if ( it.value.startsWith("http") )
				{
					def path_dir = Paths.get(params.out_dir, map["name"], "reference")
					def dir = new File( path_dir.toString() )
					dir.mkdirs()

					def matrix = new File(Paths.get(dir.toString(), "matrix.mtx.gz").toString())
					matrix.bytes = new byte[0]
					def url_matrix = new URL(it.value + "/matrix.mtx.gz")
					matrix << url_matrix.getBytes()

					def barcodes = new File(Paths.get(dir.toString(), "barcodes.tsv.gz").toString())
					barcodes.bytes = new byte[0]
					def url_barcodes = new URL(it.value + "/barcodes.tsv.gz")
					barcodes << url_barcodes.getBytes()

					def features = new File(Paths.get(dir.toString(), "features.tsv.gz").toString())
					features.bytes = new byte[0]
					def url_features = new URL(it.value + "/features.tsv.gz")
					features << url_features.getBytes()

					def types = new File(Paths.get(dir.toString(), "types.tsv.gz").toString())
					types.bytes = new byte[0]
					def url_types = new URL(it.value + "/types.tsv.gz")
					types << url_types.getBytes()

					new_map.put( it.key , dir.getName() )
					files.add(dir)
				}
				else
				{
					new_map.put( it.key , new File(it.value).getName() )
					files.add( Paths.get(it.value) )

				}
			}
			else
			{
				new_map.put( it.key , new File(it.value).getName() )
				files.add( Paths.get(it.value) )
			}
		}
		else
		{
			new_map.put(it.key, it.value)
		}
	}

	return [ new_map , files ]
}

///////////////////////////////////////////////////////////////////////////////
//// PROCESSES ////////////////////////////////////////////////////////////////

////////////
// nbconvert
nbconvert_template = Channel.fromPath("assets/jupyter-nbconvert/template.tpl")

/////////////
// parameters
include { params_json } from "./modules/process/params.nf"

/////////
// seurat
include { seurat } from "./modules/process/notebooks.nf"
render_seurat = Channel.fromPath("bin/notebooks/seurat/render.py")
j2_seurat = Channel.fromPath("bin/notebooks/seurat/j2")

/////////
// scanpy
include { scanpy } from "./modules/process/notebooks.nf"
render_scanpy = Channel.fromPath("bin/notebooks/scanpy/render.py")
j2_scanpy = Channel.fromPath("bin/notebooks/scanpy/j2")

///////////////
// jupyter_book
include { jupyter_book } from "./modules/process/report.nf"
jupyter_book_conf = Channel.fromPath("assets/jupyter-book/conf.py")
jupyter_book_logo = Channel.fromPath("assets/jupyter-book/logo.jpg")
jupyter_book_index = Channel.fromPath("assets/jupyter-book/index.md")

///////////////////////////////////////////////////////////////////////////////
//// DESIGN ///////////////////////////////////////////////////////////////////

Channel
	.fromPath(params.params)
	.splitCsv(header: true)
	.map{ getPaths(it) }
	.map{ [ addValue(it[0], "project", params.project) , it[1] ] }
	.map{ [ addValue(it[0], "scientist", params.scientist) , it[1] ] }
	.set{SAMPLES}

///////////////////////////////////////////////////////////////////////////////
//// MAIN WORKFLOW ////////////////////////////////////////////////////////////

workflow {

	SAMPLES.view()

	//params_json(SAMPLES)

	//params_json.out.view()

	//seurat( params_json.out.combine(render_seurat).combine(j2_seurat) )
	//scanpy(
	//	params_json.out
	//		.combine(render_scanpy)
	//		.combine(j2_scanpy)
	//		.combine(nbconvert_template)
	//)

	//seurat
	//	.out.ipynb
	//	.concat(scanpy.out.ipynb)
	//	.collect{ it[3] }
	//	.map{[
	//		["project": params.project, "scientist": params.scientist],
	//		it
	//	]}
	//	.combine(jupyter_book_conf)
	//	.combine(jupyter_book_index)
	//	.combine(jupyter_book_logo)
	//	.set{TO_JUPYTER_BOOK}

	//jupyter_book(TO_JUPYTER_BOOK)
	
}

