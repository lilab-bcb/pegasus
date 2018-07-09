FROM continuumio/miniconda3:latest	
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install --no-install-recommends -y \
		build-essential \
		automake \
		zlib1g-dev \
		python3-igraph \
		libxml2-dev \
		cmake \
		gnupg \
		lsb-release \
		libfftw3-dev \
		default-jdk && \
	export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
	echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
	curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
	apt-get update && apt-get install -y google-cloud-sdk

RUN pip install --upgrade pip && \
	pip install numpy && \
	pip install Cython && \
	pip install pybind11 && \
	mkdir /software && \
	git clone https://github.com/bli25broad/hnsw.git /software/hnswlib && \
	cd /software/hnswlib/python_bindings && python setup.py install && cd ../../..

COPY scrtools-0.3.0.tar.gz /software
RUN tar -xzf /software/scrtools-0.3.0.tar.gz -C /software && \
	cd /software/scrtools-0.3.0 && pip install . && cd ../..
