FROM openjdk:7

# build-essential
RUN apt-get update
RUN apt-get -y install build-essential \
    && apt-get -y install git \
    && apt-get -y install curl \
    && apt-get -y install software-properties-common

# install miniconda3 and python3
RUN apt-get update --fix-missing \
    && apt-get install -y wget bzip2 ca-certificates \
    && wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /opt/conda \
    && rm /tmp/miniconda.sh
ENV PATH /opt/conda/bin:$PATH
RUN conda update conda

# install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/cutrun_env/bin:$PATH
RUN touch .Rprofile
RUN touch .Renviron

# install CUT&RUNTools
RUN git clone https://bitbucket.org/qzhudfci/cutruntools.git
WORKDIR /cutruntools
RUN sh ucsc-tools.install
RUN apt-get install -y zlib1g-dev
RUN sh make_kseq_test.sh
ENV PATH /root/.local/bin:$PATH
ENV PATH /cutruntools:$PATH

# install SEACR
WORKDIR /
RUN git clone https://github.com/FredHutch/SEACR.git