FROM ubuntu:22.04
MAINTAINER Mikhail Kolmogorov, mikolmogorov@gmail.com

# update and install dependencies
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
    apt-get -y install cmake git make gcc g++ autoconf bzip2 lzma-dev zlib1g-dev tabix libbz2-dev && \
	apt-get -y install libcurl4-openssl-dev libpthread-stubs0-dev liblzma-dev libhdf5-dev && \
	apt-get -y install python3-pip python3-virtualenv virtualenv python3-dev && \
	apt-get -y install wget libz-dev libncurses5-dev libgsl-dev && \
	apt-get -y install graphviz libgraphviz-dev pkg-config && \
	apt-get clean && \
	apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python3 --version && \
	python3 -m pip install --upgrade pip

ARG MM_VER=2.28
RUN wget https://github.com/lh3/minimap2/releases/download/v$MM_VER/minimap2-$MM_VER.tar.bz2 && \
	tar xvf minimap2-$MM_VER.tar.bz2 && \
	rm minimap2-$MM_VER.tar.bz2 && \
	cd minimap2-$MM_VER && \
	make && \
	cp minimap2 /usr/bin/

### samtools
# 1.15
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2 && \
    tar xvf samtools-1.15.tar.bz2 && \
	rm -r /opt/samtools/samtools-1.15.tar.bz2 && \
	cd samtools-1.15/ && \
	autoheader && \
	autoconf -Wno-header && \
	./configure && \
	make && \
	cp samtools /usr/bin/samtools

### bcftools
WORKDIR /opt/bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.15/bcftools-1.15.tar.bz2 && \
    tar xvf bcftools-1.15.tar.bz2 && \
	cd bcftools-1.15/ && \
	autoheader && \
	autoconf -Wno-header && \
	./configure --enable-libgsl && \
	make && \
	cp bcftools /usr/bin/bcftools

#build and install Flye
WORKDIR /opt/flye
RUN git clone https://github.com/fenderglass/Flye && \
	cd Flye && \
	python setup.py install

#Install strainy & dependencies
WORKDIR /opt/strainy
RUN git clone -b docker https://github.com/katerinakazantseva/strainy && \
	cd strainy && \
	python3 -m pip install -r requirements.txt && \
	python3 setup.py install


ENV PYTHONUNBUFFERED "1"
