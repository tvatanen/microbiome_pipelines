FROM ubuntu:18.04

# $ docker build . -t humann:v3.6
# *OR*
# $ gcloud builds submit . --timeout 12h --tag gcr.io/osullivan-lab/humann:v3.6
# $ docker run --rm -it -v `pwd`:/data humann:v3.0.0 /bin/bash

ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
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
RUN conda activate base
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --add channels biobakery
RUN conda install humann=3.6 -c biobakery
RUN mkdir /humann_db
RUN humann_databases --download chocophlan full /humann_db --update-config yes
RUN humann_databases --download uniref uniref90_diamond /humann_db --update-config yes
RUN humann_databases --download utility_mapping full /humann_db --update-config yes
