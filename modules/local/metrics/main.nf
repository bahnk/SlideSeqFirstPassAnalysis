
import java.nio.file.Paths

process metrics {

	label "python"

	tag { "${metadata.name}" }

	publishDir Paths.get( params.out_dir , "metrics" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename ->
			if ( filename.indexOf("metrics.json") != -1 )
			{
				"${metadata.name}.json"
			}
			else
			{
				null
			}
		}

	input:
		tuple val(metadata), path(files), path(conf), path(script)
	
	output:
		tuple val(metadata), path(files), path(conf), path("metrics.json")

	script:		

		name = metadata.name
		dge = files[0]
		spatial = files[1]

		"""
		python3 $script $dge $spatial "${name}" metrics.json
		"""
}

