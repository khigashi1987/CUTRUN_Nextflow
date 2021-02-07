FROM openjdk:7

# build-essential
RUN apt-get update
RUN apt-get -y install build-essential \
    && apt-get -y install git \
    && apt-get -y install curl \
    && apt-get -y install software-properties-common

# install miniconda2, python2.7
RUN apt-get -qq update && apt-get -qq -y install bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=2 \
    && conda update conda
ENV PATH /opt/conda/bin:$PATH

# install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /usr/local/envs/cutrun_env/bin:$PATH
RUN touch .Rprofile
RUN touch .Renviron

# install CUT&RUNTools
RUN git clone https://bitbucket.org/qzhudfci/cutruntools.git
WORKDIR /cutruntools
RUN sh atactk.install.sh
RUN sh ucsc-tools.install
RUN sh make_kseq_test.sh
ENV PATH /root/.local/bin:$PATH
ENV PATH /cutruntools:$PATH

# install SEACR
WORKDIR /
RUN git clone https://github.com/FredHutch/SEACR.git