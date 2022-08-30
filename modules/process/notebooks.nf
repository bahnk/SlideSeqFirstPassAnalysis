import java.nio.file.Paths

process seurat {

	label "python"

	tag { "${metadata.name}" }

	publishDir Paths.get( params.out_dir ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${metadata.name}/seurat/${filename}" }
		//saveAs: { filename ->
		//	if ( filename.indexOf(".ipynb") != -1 )
		//	{
		//		"notebooks/${filename}"
		//	}
		//	else if ( filename.indexOf(".html") != -1 )
		//	{
		//		"html/${filename}"
		//	}
		//	else if ( filename.indexOf("seurat") != -1 )
		//	{
		//		"${metadata.name}/${filename}"
		//	}
		//}

	input:
		tuple val(metadata), path(files), path(conf), path(render), path(j2)
	
	output:
		tuple val(metadata), path(files), path(conf), path("seurat.${metadata.name}.ipynb"), emit: ipynb
		tuple val(metadata), path(files), path(conf), path("seurat.${metadata.name}.html"), emit: html
		tuple val(metadata), path(files), path(conf), path("output"), emit: output

	script:		

		metadata["notebook"] = "seurat"

		"""
		python3 "${render}" "${conf}" "seurat.${metadata.name}.ipynb"

		#xvfb-run --auto-servernum --server-num=1 \
		#jupyter-nbconvert \
		#	--execute \
		#	--ExecutePreprocessor.timeout=36000 \
		#	--inplace \
		#	"seurat.${metadata.name}.ipynb"
		"""

		if ( metadata["execute"] )
		{
			"""
			python3 "${render}" "${conf}" "seurat.${metadata.name}.ipynb"

			xvfb-run --auto-servernum --server-num=1 \
				jupyter-nbconvert \
					--execute \
					--ExecutePreprocessor.timeout=36000 \
					--inplace \
					"seurat.${metadata.name}.ipynb"

			jupyter-nbconvert \
				--to html \
				"seurat.${metadata.name}.ipynb"
			"""
		}

		else
		{
			"""
			python3 "${render}" "${conf}" "seurat.${metadata.name}.ipynb"
			jupyter-nbconvert --to html "seurat.${metadata.name}.ipynb"
			"""
		}
}

process scanpy {

	label "python"

	tag { "${metadata.name}" }

	publishDir Paths.get( params.out_dir ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${metadata.name}/scanpy/${filename}" }
		//saveAs: { filename ->
		//	if ( filename =~ /ipynb$/ )
		//	{
		//		"notebooks/${filename}"
		//	}
		//	else if ( filename =~ /html$/ )
		//	{
		//		"html/${filename}"
		//	}
		//	else if ( filename ==~ /^scanpy$/ )
		//	{
		//		"${metadata.name}/${filename}"
		//	}
		//}

	input:
		tuple val(metadata), path(files), path(conf), path(render), path(j2), path(nbconvert_tpl)
	
	output:
		tuple val(metadata), path(files), path(conf), path("scanpy.${metadata.name}.ipynb"), emit: ipynb
		tuple val(metadata), path(files), path(conf), path("scanpy.${metadata.name}.html"), emit: html
		tuple val(metadata), path(files), path(conf), path("output"), emit: output

	script:		

		metadata["notebook"] = "scanpy"

		if ( metadata["execute"] )
		{
			"""
			python3 "${render}" "${conf}" "scanpy.${metadata.name}.ipynb"

			jupyter-nbconvert \
				--execute \
				--ExecutePreprocessor.timeout=36000 \
				--inplace \
				"scanpy.${metadata.name}.ipynb"

			jupyter-nbconvert \
				--to html \
				--template-file "${nbconvert_tpl}" \
				"scanpy.${metadata.name}.ipynb"
			"""
		}

		else
		{
			"""
			python3 "${render}" "${conf}" "scanpy.${metadata.name}.ipynb"
			jupyter-nbconvert --to html "scanpy.${metadata.name}.ipynb"
			"""
		}
}

