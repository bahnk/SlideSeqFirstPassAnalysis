
# Details about input file formats

For each sample, the pipeline needs 3 files whose the locations are passed as parameters:

 * `path_dge`: the count matrix of the sample (beads versus genes)
 * `path_spatial`: the spatial coordinates of each bead
 * `path_reference`: a reference single dataset (a count matrix with cell types)

An example of each of these can be found in the `test/data` folder:

```
$ ls test/data/sample1/ # path_dge
barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz

$ ls test/data/sample1.csv # path_spatial
test/data/sample1.csv

$ ls test/data/reference/
barcodes.tsv.gz  features.tsv.gz  types.tsv.gz
```

The `test/data/reference` normally contain a `matrix.mtz.gz`.
However, it is too big to be put on GitHub.
So, it can be found [here](https://bioinformatics.crick.ac.uk/shiny/users/bahn/slideseq/test_data/reference/matrix.mtx.gz) instead.

The sample count matrix (`path_dge`) and the reference count matrix (`path_reference`) are both under the [10X file format](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices).
So, it means there are just a folder containing 3 files:

 * `barcodes.tsv.gz`: the list of the barcodes (for cells or beads)
 * `features.tsv.gz`: the list of the genes
 * `matrix.mtx.gz`: the count matrix (barcodes versus genes), under the [matrix market exchange format](https://math.nist.gov/MatrixMarket/formats.html)

Here is what these files look like:

```
$ zcat test/data/sample1/barcodes.tsv.gz | head -n 4
AAAAAAATCTTAGT
AAAAAAATTCGGTG
AAAAAACATCTTTC
AAAAAACGAAATAG

$ zcat test/data/sample1/features.tsv.gz | head -n 4
ENSMUSG00000109644      0610005C13Rik   Gene Expression
ENSMUSG00000108652      0610006L08Rik   Gene Expression
ENSMUSG00000007777      0610009B22Rik   Gene Expression
ENSMUSG00000086714      0610009E02Rik   Gene Expression

$ zcat test/data/sample1/matrix.mtx.gz | head -n 6
%%MatrixMarket matrix coordinate integer general
%metadata_json: {"institute": "The Crick Institute"}
54838 42236 2929086
1 2257 1
1 2830 1
1 11275 1
```

Please note that the first column of `features.tsv.gz` is the gene ID and the second one is the gene symbol.
It has to be this way because this is what the [`read_10x_mtx`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_10x_mtx.html) and [`Read10X`](https://satijalab.org/seurat/reference/read10x) functions are expecting.
The same comment applies to the reference dataset.

The single-cell reference dataset has an additional file named `types.tsv.gz` which stores the cell types associated with each barcodes:

```
$ zcat test/data/reference/types.tsv.gz | head
Ependymal
CA3
CA3
CA3
```

Finally, the `path_spatial` file stores the coordinates of each bead barcode:

```
$ head -n 4 test/data/sample1.csv
Barcode,x,y
AAAACAATAAAAGG,2451.4,1092.5
AAAATATAAAAACC,3789.6,3876.1
AAATCAAAAAAATG,3905.5,3863.3
```

