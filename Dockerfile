# Dockerfile for testing scCloud

FROM continuumio/miniconda3:4.6.14
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y \
		build-essential \
		automake \
		zlib1g-dev \
		python3-igraph \
		libxml2-dev \
		cmake \
		libfftw3-dev \
		default-jre \
		git

RUN apt-get clean && \
     rm -rf /var/lib/apt/lists/*

RUN pip install cython pybind11 numpy
RUN pip install fitsne xlrd

COPY . /scCloudPy/
WORKDIR /scCloudPy/tests
RUN git clone https://github.com/klarman-cell-observatory/scCloud-test-data.git
WORKDIR /scCloudPy/
RUN pip install -e .


