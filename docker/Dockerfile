# MiTo image

################################################################

## Start Debian basic image: NBBB ricordarsi --platform linux/x86_64 in docker build
FROM python:3.9-bullseye

## Change as u want
LABEL mantainer = "andrea.cossa@ieo.it"

## Compilers, utilities, git
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

## Mamba
RUN wget https://github.com/conda-forge/miniforge/releases/download/24.11.3-0/Miniforge3-24.11.3-0-Linux-x86_64.sh \
    && bash Miniforge3-24.11.3-0-Linux-x86_64.sh -b -p /mambaforge
ENV PATH="/mambaforge/bin/:$PATH"

## Clone repo and install code
RUN git clone https://github.com/andrecossa5/MiTo.git \
    && cd MiTo \
    && mamba env update -n base -f envs/environment.yml \
    && cd .. \
    && pip install . 

## Default command to initialize the container
CMD [ "bash" ]

#################################################################
