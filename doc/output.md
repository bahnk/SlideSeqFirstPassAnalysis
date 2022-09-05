
# Pipeline output

The output of the pipeline should be in the directory specified by the `out_dir` parameter in the [parameter file](config.md).

Each sample has its own directory based on its name.
Additionally, there are a `params` directory containing the parameters as a `JSON` for each sample, and `_build` directory containing the output of `jupyter-book`:

```bash
$ ls -l results 
_build  params  sample1  sample2

$ ls -l results/params
sample1.json  sample2.json
```

The `_build` directory is self-contained html book.
It can be open with the `results/_build/html/index.html` file.

In each sample directory, you should find the output files of each no detailed [here](steps.md) organised by order.

```
$ ls results/sample1
counts  destvi  rctd  reference  scanpy  seurat  sparkx  spatial.csv
```

 * `counts`, `spatial.csv` and `reference`: the files passed as input in the [CSV configuration file](config.md) (`path_dge`, `path_spatial` and `path_reference` respectively)
 * `destvi`, `rctd`, `scanpy`, `seurat` and `sparkx`: the notebook files with their input and output for each analysis as detailed [here](steps.md).

