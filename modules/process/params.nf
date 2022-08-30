import java.nio.file.Paths
import groovy.json.JsonOutput

process params_json {

	tag { "${metadata.name}" }

	publishDir Paths.get( params.out_dir ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename ->
			if ( filename.indexOf(".json") != -1 )
			{
				"params/${filename}"
			}
		}

	input:
		tuple val(metadata), path(files)
	
	output:
		tuple val(metadata), path(files), path("${metadata.name}.json")

	script:		
		
		def map = [:]
		metadata.each{
			map[it.key] = it.value.toString()
		}
		def json = JsonOutput.toJson(map)
		def pretty = JsonOutput.prettyPrint(json)

		"""
		echo -E '${pretty}' > "${metadata.name}.json"
		"""
}
