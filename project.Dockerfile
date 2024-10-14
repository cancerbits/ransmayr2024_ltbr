# pull base image
FROM rocker/tidyverse:4.1.2

# who maintains this image
LABEL maintainer Christoph Hafemeister "christoph.hafemeister@ccri.at"
LABEL version 4.1.2-v1

# Change the default CRAN mirror
RUN echo "options(repos = c(CRAN = 'https://mran.microsoft.com/snapshot/2022-02-01'), download.file.method = 'libcurl')" >> ${R_HOME}/etc/Rprofile.site

RUN apt-get update -qq \
  && apt-get -y --no-install-recommends install \
  htop \
  nano \
  libigraph-dev \
  libcairo2-dev \
  libxt-dev \
  libcurl4-openssl-dev \
  libcurl4 \
  libxml2-dev \
  libxt-dev \
  openssl \
  libssl-dev \
  wget \
  curl \
  bzip2 \
  libbz2-dev \
  libpng-dev \
  libhdf5-dev \
  pigz \
  libudunits2-dev \
  libgdal-dev \
  libgeos-dev \
  libboost-dev \
  libboost-iostreams-dev \
  libboost-log-dev \
  libboost-system-dev \
  libboost-test-dev \
  libz-dev \
  libarmadillo-dev \
  libglpk-dev \
  jags \
  libgsl-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  cmake \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# add nlopt (for muscat R package)
RUN git clone -b v2.7.1 --single-branch https://github.com/stevengj/nlopt.git \
  && cd nlopt \
  && mkdir build \
  && cd build \
  && cmake .. \
  && make \
  && make install \
  && cd ../.. \
  && rm -rf nlopt

RUN R -e "update.packages(ask = FALSE)"

RUN install2.r --error \
    DT \
    profvis \
    tictoc \
    markdown \
    plotly \
    bench \
    mclust 

RUN R -e "BiocManager::install(version = "3.14")"
RUN R -e "BiocManager::install(c('multtest', 'limma', 'Biobase', 'monocle', 'rtracklayer', 'IRanges', 'GenomeInfoDb', 'GenomicRanges', 'BiocGenerics', 'DESeq2', 'MAST', 'SingleCellExperiment', 'SummarizedExperiment', 'S4Vectors', 'splatter', 'GEOquery', 'infercnv', 'glmGamPoi', 'hypeR', 'fgsea', 'ComplexHeatmap'), update = FALSE)"
RUN R -e "remotes::install_github(repo = 'satijalab/sctransform', ref = 'e9e52a4', dependencies = TRUE)"
RUN R -e "remotes::install_github(repo = 'satijalab/seurat', ref = 'v4.1.0', dependencies = TRUE)"
RUN R -e "remotes::install_github(repo = 'satijalab/azimuth', ref = 'v0.4.3', dependencies = TRUE)"

RUN chmod -R a+rw ${R_HOME}/site-library # so that everyone can dynamically install more libraries within container
RUN chmod -R a+rw ${R_HOME}/library

# add custom options for rstudio sessions
# make sure sessions stay alive forever
RUN echo "session-timeout-minutes=0" >> /etc/rstudio/rsession.conf
# make sure we get rstudio server logs in the container
# RUN echo $'[*]\nlog-level=warn\nlogger-type=file\n' > /etc/rstudio/logging.conf
# make sure authentication is not needed so often
RUN echo "auth-timeout-minutes=0" >> /etc/rstudio/rserver.conf
RUN echo "auth-stay-signed-in-days=365" >> /etc/rstudio/rserver.conf

# make sure all R related binaries are in PATH in case we want to call them directly
ENV PATH ${R_HOME}/bin:$PATH

# Install Mambaforge (for snakemake) with conda 
RUN mkdir -p /opt/mambaforge \
  && wget -O /opt/mambaforge/mambaforge.sh https://github.com/conda-forge/miniforge/releases/download/4.10.1-5/Mambaforge-4.10.1-5-Linux-x86_64.sh \
  && chmod +x /opt/mambaforge/mambaforge.sh \
  && /opt/mambaforge/mambaforge.sh -u -b -p /opt/mambaforge \
  && rm /opt/mambaforge/mambaforge.sh
ENV PATH /opt/mambaforge/condabin:$PATH
ENV CONDA_AUTO_UPDATE_CONDA=false
RUN conda init bash

RUN R -e "install.packages('treemap', repo = 'https://packagemanager.posit.co/cran/__linux__/jammy/2023-08-14')"
