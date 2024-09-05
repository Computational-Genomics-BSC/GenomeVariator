FROM ubuntu:20.04
LABEL author="Rodrigo Martin <rodrigo.martin@bsc.es>"

# Install dependencies
RUN apt-get update && apt-get upgrade -y && DEBIAN_FRONTEND=noninteractive apt-get install -y \
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
RUN pip install pysam variant-extractor

# Velvet
RUN wget https://github.com/dzerbino/velvet/archive/refs/tags/v1.2.10.tar.gz && tar -xvzf v1.2.10.tar.gz
RUN make -C velvet-1.2.10
RUN cp velvet-1.2.10/velvetg $HOME/bin && cp velvet-1.2.10/velveth $HOME/bin

# Exonerate
RUN git clone https://github.com/adamewing/exonerate.git
RUN cd exonerate && autoreconf -fi  && ./configure && make && make install

RUN export PATH=$PATH:$HOME/bin

# Copy the content of the repository to the folder GenomeVariator
COPY dependencies/ /GenomeVariator/
COPY src/ /GenomeVariator/
COPY utils/ /GenomeVariator/
