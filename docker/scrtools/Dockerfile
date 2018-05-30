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
		lsb-release && \
	export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
	echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
	curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
	apt-get update && apt-get install -y google-cloud-sdk

RUN pip install --upgrade pip && \
	pip install numpy==1.13.3 && \
	pip install louvain && \
	pip install scanpy==0.4.4

RUN mkdir /software

RUN git clone https://github.com/DmitryUlyanov/Multicore-TSNE.git /software/Multicore-TSNE && \
	cd /software/Multicore-TSNE && pip install . && cd ../..

COPY scrtools-0.2.0.tar.gz /software
RUN tar -xzf /software/scrtools-0.2.0.tar.gz -C /software && \
	cd /software/scrtools-0.2.0 && pip install . && cd ../..


pip install numpy
pip install Cython
conda install -c conda-forge fftw
export CPATH=$CPATH:/ahg/regevdata/users/libo/miniconda3/envs/scrtools/include