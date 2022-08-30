import java.nio.file.Paths
import groovy.json.JsonOutput

process jupyter_book {

	label "python"

	publishDir Paths.get( params.out_dir ), mode: "copy", overwrite: "true"

	input:
		tuple val(metadata), path(notebooks), path(conf), path(index), path(logo)
	
	output:
		path "_build"

	script:		

		def map = [:]
		map["project"] = metadata["project"]
		map["scientist"] = metadata["scientist"]
		def json = JsonOutput.toJson(map)
		def pretty = JsonOutput.prettyPrint(json)

		"""
		echo -E '${pretty}' > "args.json"
		python3 "${conf}" "." "." "args.json"
		jupyter-book build .
		"""
}

