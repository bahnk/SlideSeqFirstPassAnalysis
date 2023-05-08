import java.nio.file.Paths

process notebook_py {

	label "notebook"
	label "python"

	tag { "${metadata.name}" }

	publishDir Paths.get( params.out_dir ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${metadata.name}/${prefix}/${filename}" }

	input:
		tuple \
			val(metadata),
			path(files),
			path(conf),
			val(prefix),
			path(render),
			path(j2),
			path(nbconvert_tpl)
	
	output:
		tuple val(metadata), path(files), path(conf), path("${basename}.ipynb"), emit: ipynb
		tuple val(metadata), path(files), path(conf), path("${basename}.rmd"), emit: rmd
		tuple val(metadata), path(files), path(conf), path("${basename}.html"), emit: html
		tuple val(metadata), path(files), path(conf), path("output"), emit: output

	script:		

		basename = "${prefix}.${metadata.name}"

		r_cmd = 'library(rmarkdown);rmarkdown::convert_ipynb("'
		r_cmd += "${prefix}.${metadata.name}.ipynb"
		r_cmd += '","'
		r_cmd += "${prefix}.${metadata.name}.rmd"
		r_cmd += '")'

		if ( metadata.execute )
		{
			"""
			python3 "${render}" "${conf}" "${basename}.ipynb"

			jupyter-nbconvert \
				--execute \
				--ExecutePreprocessor.timeout=36000 \
				--inplace \
				"${basename}.ipynb"

			jupyter-nbconvert \
				--to html \
				--template-file "${nbconvert_tpl}" \
				"${basename}.ipynb"

			R -e '${r_cmd}'
			"""
		}

		else
		{
			"""
			python3 "${render}" "${conf}" "${basename}.ipynb"
			jupyter-nbconvert --to html "${basename}.ipynb"
			mkdir output
			R -e '${r_cmd}'
			"""
		}
}

process notebook_py_gpu {

	label "notebook"
	label "python"
	label "gpu"

	tag { "${metadata.name}" }

	publishDir Paths.get( params.out_dir ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${metadata.name}/${prefix}/${filename}" }

	input:
		tuple \
			val(metadata),
			path(files),
			path(conf),
			val(prefix),
			path(render),
			path(j2),
			path(nbconvert_tpl)
	
	output:
		tuple val(metadata), path(files), path(conf), path("${basename}.ipynb"), emit: ipynb
		tuple val(metadata), path(files), path(conf), path("${basename}.rmd"), emit: rmd
		tuple val(metadata), path(files), path(conf), path("${basename}.html"), emit: html
		tuple val(metadata), path(files), path(conf), path("output"), emit: output

	script:		

		basename = "${prefix}.${metadata.name}"

		r_cmd = 'library(rmarkdown);rmarkdown::convert_ipynb("'
		r_cmd += "${basename}.ipynb"
		r_cmd += '","'
		r_cmd += "${basename}.rmd"
		r_cmd += '")'

		if ( metadata.execute )
		{
			"""
			python3 "${render}" "${conf}" "${basename}.ipynb"

			jupyter-nbconvert \
				--execute \
				--ExecutePreprocessor.timeout=36000 \
				--inplace \
				"${basename}.ipynb"

			jupyter-nbconvert \
				--to html \
				--template-file "${nbconvert_tpl}" \
				"${basename}.ipynb"

			R -e '${r_cmd}'
			"""
		}

		else
		{
			"""
			python3 "${render}" "${conf}" "${basename}.ipynb"
			jupyter-nbconvert --to html "${basename}.ipynb"
			mkdir output
			R -e '${r_cmd}'
			"""
		}
}

process notebook_r {

	label "notebook"
	label "python_and_r"

	tag { "${metadata.name}" }

	publishDir Paths.get( params.out_dir ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${metadata.name}/${prefix}/${filename}" }

	input:
		tuple \
			val(metadata),
			path(files),
			path(conf),
			val(prefix),
			path(render),
			path(j2),
			path(nbconvert_tpl)
	
	
	output:
		tuple val(metadata), path(files), path(conf), path("${basename}.ipynb"), emit: ipynb
		tuple val(metadata), path(files), path(conf), path("${basename}.rmd"), emit: rmd
		tuple val(metadata), path(files), path(conf), path("${basename}.html"), emit: html
		tuple val(metadata), path(files), path(conf), path("output"), emit: output

	script:		

		basename = "${prefix}.${metadata.name}"

		r_cmd = 'library(rmarkdown);rmarkdown::convert_ipynb("'
		r_cmd += "${basename}.ipynb"
		r_cmd += '","'
		r_cmd += "${basename}.rmd"
		r_cmd += '")'

		if ( metadata.execute )
		{
			"""
			python3 "${render}" "${conf}" "${basename}.ipynb"

			xvfb-run --auto-servernum --server-num=1 \
				jupyter-nbconvert \
					--execute \
					--ExecutePreprocessor.timeout=36000 \
					--inplace \
					"${basename}.ipynb"

			jupyter-nbconvert \
				--to html \
				--template-file "${nbconvert_tpl}" \
				"${basename}.ipynb"

			R -e '${r_cmd}'
			"""
		}

		else
		{
			"""
			python3 "${render}" "${conf}" "${basename}.ipynb"
			jupyter-nbconvert --to html "${basename}.ipynb"
			mkdir output
			R -e '${r_cmd}'
			"""
		}
}

