FROM ubuntu:20.04
LABEL author="Rodrigo Martin <rodrigo.martin@bsc.es>"

ENV PATH=$PATH:$HOME/bin
ARG DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    python3 \
    python3-pip \
    git \
    wget \
    build-essential \
    libz-dev \
    libglib2.0-dev \
    libbz2-dev \
    liblzma-dev \
    default-jre \
    autoconf \
    samtools \
    biobambam2 \
    bwa


RUN mkdir $HOME/bin

# Python dependencies
RUN pip install pysam

# art_illumina
RUN wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
RUN tar -xvzf artbinmountrainier2016.06.05linux64.tgz
RUN cp art_bin_MountRainier/art_illumina $HOME/bin

# Velvet
RUN wget https://github.com/dzerbino/velvet/archive/refs/tags/v1.2.10.tar.gz && tar -xvzf v1.2.10.tar.gz
RUN make -C velvet-1.2.10
RUN cp velvet-1.2.10/velvetg $HOME/bin && cp velvet-1.2.10/velveth $HOME/bin

# Exonerate
RUN git clone https://github.com/adamewing/exonerate.git
RUN cd exonerate && autoreconf -fi  && ./configure && make && make install

RUN export PATH=$PATH:$HOME/bin

RUN git clone --recurse-submodules --remote-submodules https://gitlab.bsc.es/rmarti1/genome-variator
