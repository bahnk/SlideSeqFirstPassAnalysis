
# Slide-seq first pass analysis

Here is the documentation:

 1. [Pipeline steps](doc/steps.md)
 2. [Pipeline configuration](doc/config.md)
 3. [Running the pipeline](doc/run.md)
 4. [Pipeline output](doc/output.md)

In a nutshell, first be sure that your singularity config directory is not in your home.
For example:

```
$ ls -l ~/.singularity
lrwxrwxrwx 1 username domain users 40 Aug 26  2021 /camp/home/username/.singularity -> /camp/stp/babs/working/username/.singularity
```

Then, to test the pipeline it you can just run:

```
# create a new directory and cd in it
mkdir test
cd test

# download the example parameters files
wget https://bioinformatics.crick.ac.uk/shiny/users/bahn/slideseq/test_data/params.yml
wget https://bioinformatics.crick.ac.uk/shiny/users/bahn/slideseq/test_data/params.csv

# load nextflow and singularity
module load Nextflow/22.04.0 Singularity/3.6.4

# pull the latest version
nextflow pull bahnk/SlideSeqFirstPassAnalysis -r main

# run the pipeline and pray
nextflow run bahnk/SlideSeqFirstPassAnalysis -r main -params-file params.yml
```

