# modified from https://github.com/biobakery/biobakery/blob/master/docker/metaphlan/Dockerfile
FROM ubuntu:18.04

# also install python version 2 used by bowtie2
RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y python python-dev python3 python3-dev python3-pip apt-transport-https openjdk-8-jre wget zip
RUN pip3 install boto3 cloudpickle awscli
RUN apt-get install -y bowtie2
RUN pip3 install numpy
RUN pip3 install cython
RUN pip3 install biom-format

# install cmseq
RUN apt-get install -y git 
RUN git clone https://github.com/SegataLab/cmseq.git
RUN cd cmseq && \
    python3 setup.py install && \
    cd ../ && \
    rm -r cmseq

RUN pip3 install metaphlan==3.0.13
RUN metaphlan --install

WORKDIR /tmp