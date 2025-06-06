# nf_scLT image

################################################################

# Start from Debian basic image: NBBB remember --platform linux/x86_64 in docker build
FROM python:3.11-bullseye
LABEL mantainer = "andrea.cossa@ieo.it"

# Compilers, utilities, git
RUN apt-get --quiet update --yes \
    && apt-get --quiet upgrade --yes \
    && apt-get --quiet install --yes \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      build-essential \
      libz-dev \
      wget \
      procps \
      git \
      cmake \
    && rm -rf /var/lib/apt/lists/*

# Mamba
RUN wget https://github.com/conda-forge/miniforge/releases/download/24.11.3-0/Miniforge3-24.11.3-0-Linux-x86_64.sh \
    && bash Miniforge3-24.11.3-0-Linux-x86_64.sh -b -p /mambaforge
ENV PATH="/mambaforge/bin/:$PATH"

# Python env
RUN mamba install -c conda-forge pigz=2.8 python=3.10 numpy=1.22.4 pandas=1.5.3 scipy=1.9.0
RUN pip install plotting_utils

# Genomic tools
RUN mamba install -c bioconda --yes samtools=1.20 pysam=0.22.1 bcftools=1.20 fgbio=2.1.0 cellsnp-lite bwa picard freebayes

# STAR, from source
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.tar.gz \
    && tar -xzf 2.7.9a.tar.gz \
    && cd STAR-2.7.9a/source \
    && make
ENV PATH="/STAR-2.7.9a/bin/Linux_x86_64/:$PATH"

# Default command to initialize the container
CMD [ "bash" ]

#################################################################

