import java.nio.file.Paths
import nextflow.util.ArrayBag

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

