
import groovy.json.JsonSlurper
import java.nio.file.Files
import java.nio.file.Paths
import nextflow.util.ArrayBag
import org.yaml.snakeyaml.Yaml

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

//////////////////////
def toBoolean(map) {//
//////////////////////

	def new_map = [:]

	map.each{
		if ( ["true", "false"].contains( it.value.toString().toLowerCase() ) )
		{
			if ( "true" == it.value.toString().toLowerCase() )
			{
				new_map.put(it.key, true)
			}
			else
			{
				new_map.put(it.key, false)
			}
		}
		else
		{
			new_map.put(it.key, it.value)
		}
	}

	return new_map
}

/////////////////////////
def getDataPaths(map) {//
/////////////////////////

	def new_map = [:]
	def files = []

	map.each{

		// count matrix
		if ( it.key == "path_dge" )
		{
			if ( it.value.startsWith("http") )
			{
				def path_dir = Paths.get(params.out_dir, map["name"], "counts")
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

				new_map.put( it.key , dir.getName() )
				files.add( Paths.get(dir.getAbsolutePath()) )
			}

			else
			{
				new_map.put( it.key , new File(it.value).getName() )
				files.add( Paths.get(it.value) )
			}
		}

		else if ( it.key == "path_spatial" )
		{
			if ( it.value.startsWith("http") )
			{
				def spatial = new File(Paths.get(params.out_dir, map["name"], "spatial.csv").toString())
				spatial.bytes = new byte[0]
				def url = new URL(it.value)
				spatial << url.getBytes()

				new_map.put( it.key , spatial.getName() )
				files.add( Paths.get(spatial.getAbsolutePath()) )
			}

			else
			{
				new_map.put( it.key , new File(it.value).getName() )
				files.add( Paths.get(it.value) )
			}
		}

		// reference
		else if ( it.key == "path_reference" )
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
				files.add( Paths.get(dir.getAbsolutePath()) )
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

	return [ new_map , new ArrayBag(files) ]
}

//////////////////////
def readYaml(path) {//
//////////////////////
	return new Yaml().load( new File(path).text )
}

/////////////////////////////
def checkAllParams(item) { //
/////////////////////////////

	/*
	 * we want to:
	 *  1. check if the parameters of a process were defined
	 *  2. if none of them were defined we don't execute the process
	 *  3. if only some of them are missing we set them to default values
	 */

	def metadata = item[0]
	def paths = item[2]

	// the first parameter is the metadata
	// the second parameter is the array of paths (DGE, spatial and reference)
	// the last parameter is the array of the paths of the process info
	// So, we'll return the first and second parameters

	// the general params
	metadata = checkGeneralParams(metadata, paths)

	// the specific params
	for ( path in paths )
	{
		metadata = checkSpecificParams(metadata, path.toString())
	}
	
	return [ metadata , item[1] ]
}

def checkGeneralParams(metadata, paths) {

	def general_params = [] as Set

	for ( path in paths )
	{
		def info = readYaml(path.toString())
		def params = info.parameters.findAll{ ! it.specific }.collect{ it.name }
		general_params.addAll(params)
	}

	def missing = [] as Set

	for ( param in general_params )
	{
		if ( metadata.parameters.containsKey(param) )
		{
			if ( ! metadata.parameters[param] )
			{
				missing = missing + param
			}
		}
		
		else
		{
			missing = missing + param
		}
	}

	if ( metadata.containsKey("missing_general_parameters") )
	{
		metadata.missing_parameters = metadata.missing_parameters.addAll(missing)
	}

	else
	{
		metadata["missing_general_parameters"] = missing
	}

	return metadata
}

///////////////////////////////////////////
def checkSpecificParams(metadata, path) {//
///////////////////////////////////////////

	def process_info = readYaml(path)
	def process_name = process_info.name
	def process_params = process_info.parameters

	def defined_params = metadata.parameters
	def execution = metadata.execution

	def params = [:]
	def status = [:]

	// we attribute each parameter with a flag:
	// valued if it has a value set by the user
	// empty if it doesn't have a value set by the user
	// no key if there was no column for this parameter in the csv file

	process_params.each {

		if ( it.specific )
		{

			def param_name = it.name

			if ( defined_params.containsKey(param_name) )
			{
				if ( defined_params[param_name] )
				{
					status[param_name] = "valued"
					params[param_name] = defined_params[param_name]
				}

				else
				{
					status[param_name] = "empty"
					params[param_name] = it.default
				}
			}

			else
			{
				status[param_name] = "no key"
				params[param_name] = it.default
			}
		}
	}

	// check the parameters defined by the user
	specific_params = process_params.findAll{ it.specific }.collect{ it.name }
	defined_specific_params = status.findAll{ k, v -> v == "valued" }.keySet()

	def process_specific_params =
		process_params.findAll{ it.specific }.collect{ it.name }

	// if 0 parameter is defined we don't execute the process and remove their
	// key
	if ( defined_specific_params.size() == 0 )
	{
		execution[process_name] = false
		metadata.execution << execution

		new_params = [:]

		for ( entry in metadata.parameters )
		{
			if ( ! process_specific_params.contains(entry.key) )
			{
				new_params[entry.key] = entry.value
			}
		}

		metadata.parameters = new_params
	}

	else
	{
		execution[process_name] = true
		metadata.execution << execution
		metadata.parameters << params
	}

	return metadata
}

////////////////////////////////
def addMetrics(metadata, path)//
////////////////////////////////
{
	def slurper = new JsonSlurper()
	def metrics = slurper.parseText( new File(path).text )

	def new_map = metadata.clone()
	new_map["metrics"] = metrics

	return new_map
}

