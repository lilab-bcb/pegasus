FROM continuumio/miniconda3:latest	
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install --no-install-recommends -y \
		build-essential \
		automake \
		zlib1g-dev \
		python3-igraph \
		libxml2-dev \
		cmake

RUN pip install --upgrade pip
RUN pip install numpy==1.13.3
RUN pip install louvain
RUN pip install scanpy==0.4.4

RUN mkdir /software

RUN git clone https://github.com/DmitryUlyanov/Multicore-TSNE.git /software/Multicore-TSNE
RUN cd /software/Multicore-TSNE && pip install . && cd ../..

COPY scrtools-0.1.0.tar.gz /software
RUN tar -xzf /software/scrtools-0.1.0.tar.gz -C /software
RUN cd /software/scrtools-0.1.0 && pip install . && cd ../..
