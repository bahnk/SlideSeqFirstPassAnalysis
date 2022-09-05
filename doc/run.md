
# Running the pipeline

Just `ssh` to camp and open terminal demultiplexer (`tmux` or `screen`) if you use one.

Then, be sure that your singularity config directory is not in your home.
For example:

```bash
$ ls -l ~/.singularity
lrwxrwxrwx 1 username domain users 40 Aug 26  2021 /camp/home/username/.singularity -> /camp/stp/babs/working/username/.singularity
```

You can now download parameters file examples:

```
wget https://bioinformatics.crick.ac.uk/shiny/users/bahn/slideseq/test_data/params.yml
wget https://bioinformatics.crick.ac.uk/shiny/users/bahn/slideseq/test_data/params.csv
```

Now, you need to parametrise the pipeline as detailed [here](config.md).
You can either edit the [`params.yml`](../params.yml) or you can overwrite the parameters with the command line as explained later.
Normally, you should only have to change the `params` parameter.

We need to load these two modules before running the pipeline:

```bash
module load Nextflow/22.04.0 Singularity/3.6.4
```

Then, pull the latest version of the pipeline:
```bash
nextflow pull bahnk/SlideSeqFirstPassAnalysis -r main
```

Finally, you can run the pipeline this way:

```bash
nextflow run bahnk/SlideSeqFirstPassAnalysis -r main -params-file params.yml
```

Alternatively, if you don't want to edit the `params.yml` file, then you can overwrite the parameters this way:

```bash
nextflow run bahnk/SlideSeqFirstPassAnalysis -r main -params-file params.yml --out_dir /path/to/output_directory --params /path/to/csv_params_file
```

