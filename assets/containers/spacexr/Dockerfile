
FROM ubuntu:22.04

RUN apt-get update
RUN apt-get --yes upgrade
RUN apt-get install --yes apt-utils

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London
RUN apt-get install --yes tzdata

RUN apt-get install --yes less
RUN apt-get install --yes vim

RUN apt-get install --yes "pandoc*"
RUN apt-get install --yes texlive
RUN apt-get install --yes snakemake

###############################################################################
## PYTHON

RUN apt-get install --yes \
	python3-louvain \
	python3-networkx \
	python3-numpy \
	python3-pandas \
	python3-scikit-fmm \
	python3-scikit-rf \
	python3-scipy \
	python3-seaborn \
	python3-simpy \
	python3-simpy3 \
	python3-sklearn \
	python3-sklearn-lib \
	python3-sklearn-pandas \
	python3.9-numpy

RUN apt-get install --yes \
  python3-jinja-vanish \
  python3-jinja2

RUN apt-get install --yes ipython3

RUN apt-get install --yes \
	python3-ipython \
	python3-ipython-genutils \
	python3-ipyparallel

RUN apt-get install --yes "jupyter*"

RUN apt-get install --yes \
	python3-ipykernel \
	python3-ipywidgets \
	python3-jupyter-client \
	python3-jupyter-console \
	python3-jupyter-core \
	python3-jupyter-server \
	python3-jupyter-sphinx \
	python3-jupyter-sphinx-theme

RUN apt-get install --yes \
	python3-jupyterlab-pygments \
	python3-jupyterlab-server \
  python3-jupyterlab-server \
  python3-jupyterlab-pygments

RUN apt-get install --yes \
	python3-notebook \
	python3-nb2plots \
	python3-nbclient \
	python3-nbconvert \
	python3-nbformat \
	python3-nbsphinx

RUN apt-get install --yes \
	python3-oauth \
	python3-oauth2client \
	python3-oauthlib \
	python3-urllib3

RUN apt-get install --yes python3-termcolor
RUN apt-get install --yes python3-yaml

RUN apt-get install --yes python3-pip
RUN pip install upsetplot
RUN pip install umap-learn
#RUN pip install leiden
#RUN pip install leidenalg
RUN pip install anndata
RUN pip install scanpy
RUN pip install jupyter-book

###############################################################################
# R

RUN apt-get install --yes \
	r-base \
	r-base-core \
	r-base-core-dbg \
	r-base-dev \
	r-base-html \
	r-base-latex

RUN apt-get install --yes \
	r-bioc-hdf5array \
	r-bioc-rhdf5 \
	r-bioc-rhdf5filters \
	r-bioc-rhdf5lib \
	r-bioc-s4vectors \
	r-cran-argparse \
	r-cran-argparser \
	r-cran-data.table \
	r-cran-dbplyr \
	r-cran-dplyr \
	r-cran-hdf5r \
	r-cran-itertools \
	r-cran-jsonlite \
	r-cran-openxlsx \
	r-cran-optparse \
	r-cran-plyr \
	r-cran-purrr \
	r-cran-purrrlyr \
	r-cran-purrrogress \
	r-cran-r.utils \
	r-cran-rcurl \
	r-cran-readr \
	r-cran-readxl \
	r-cran-reshape \
	r-cran-reshape2 \
	r-cran-rjava \
	r-cran-rjson \
	r-cran-rlang \
	r-cran-rsqlite \
	r-cran-stringi \
	r-cran-stringr \
	r-cran-tibble \
	r-cran-tidygraph \
	r-cran-tidyr \
	r-cran-tidyselect \
	r-cran-tidyverse \
	r-cran-udunits2 \
	r-cran-yaml

RUN apt-get install --yes \
	r-cran-bookdown \
	r-cran-htmltools \
	r-cran-htmlwidgets \
	r-cran-knitr \
	r-cran-markdown \
	r-cran-rmarkdown

RUN apt-get install --yes \
	r-cran-rcpp \
	r-cran-rcppannoy \
	r-cran-rcpparmadillo

RUN apt-get install --yes \
	r-cran-leiden \
	r-cran-rtsne \
	r-cran-tsne

RUN apt-get install --yes \
	r-bioc-annotationdbi \
	r-bioc-biocinstaller \
	r-bioc-biocneighbors \
	r-bioc-biocparallel \
	r-bioc-biomart \
	r-bioc-biostrings \
	r-bioc-genomicranges \
	r-bioc-go.db \
	r-bioc-keggrest \
	r-bioc-org.hs.eg.db \
	r-bioc-rtracklayer \
	r-bioc-tximport

RUN apt-get install --yes \
	r-bioc-deseq \
	r-bioc-deseq2 \
	r-bioc-dexseq \
	r-bioc-edger \
	r-bioc-limma \
	r-bioc-summarizedexperiment

RUN apt-get install --yes \
	r-bioc-glmgampoi \
	r-bioc-monocle \
	r-bioc-singlecellexperiment \
	r-bioc-singler \
	r-cran-seurat \
	r-cran-seuratobject

RUN apt-get install --yes \
	r-bioc-complexheatmap \
	r-cran-cowplot \
	r-cran-ggalluvial \
	r-cran-ggally \
	r-cran-gganimate \
	r-cran-ggbeeswarm \
	r-cran-ggdendro \
	r-cran-ggdendro \
	r-cran-ggeffects \
	r-cran-ggforce \
	r-cran-ggfortify \
	r-cran-ggm \
	r-cran-ggplot.multistats \
	r-cran-ggplot2 \
	r-cran-ggplot2 \
	r-cran-ggpubr \
	r-cran-ggpubr \
	r-cran-ggraph \
	r-cran-ggrepel \
	r-cran-ggrepel \
	r-cran-ggridges \
	r-cran-ggsci \
	r-cran-ggseqlogo \
	r-cran-ggsignif \
	r-cran-ggtext \
	r-cran-ggthemes \
	r-cran-ggvis \
	r-cran-gridbase \
	r-cran-gridextra \
	r-cran-lattice \
	r-cran-latticeextra \
	r-cran-pheatmap \
	r-cran-plotly \
	r-cran-rcolorbrewer \
	r-cran-upsetr \
	r-cran-venndiagram \
	r-cran-vioplot \
	r-cran-viridis \
	r-cran-viridislite

RUN apt-get install --yes \
	r-cran-dendextend \
	r-cran-densityclust \
	r-cran-mclust

RUN apt-get install --yes \
	r-cran-fit.models \
	r-cran-irlba \
	r-cran-markovchain \
	r-cran-matrixmodels \
	r-cran-matrixstats

RUN apt-get install --yes r-cran-devtools

RUN apt-get install --yes r-cran-irkernel

RUN R -e "BiocManager::install('org.Mm.eg.db')"
RUN R -e "IRkernel::installspec(prefix='/usr/local')"
RUN R -e "devtools::install_github('dmcable/spacexr', build_vignettes=FALSE)"

###############################################################################

RUN pip install -U jupyterlab

