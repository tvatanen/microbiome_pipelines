FROM ubuntu:18.04

# Note that you'll need to have Bowtie2 indexed human genome in the local build directory
# (see the bottom of the Dockerfile)
# $ docker build . -t kneaddata:v0.10.0
# $ docker run --rm -it -v `pwd`:/data kneaddata:v0.10.0 /bin/bash

ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
    zlib1g-dev \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git mercurial subversion

# install conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Make RUN commands use `bash --login` (instead of '/bin/sh -c'):
SHELL ["/bin/bash", "--login", "-c"]

RUN conda init bash
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --add channels biobakery
RUN conda install -c anaconda zlib
RUN conda install kneaddata -c biobakery

RUN mkdir /hg37
RUN kneaddata_database --download human_genome bowtie2 /hg37/

