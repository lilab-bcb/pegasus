# Dockerfile for testing pegasus

FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN mkdir -p /usr/share/man/man1 && \
    apt-get update && apt-get install -y \
		build-essential \
		zlib1g-dev \
		python3-igraph \
		libxml2-dev \
		libfftw3-dev \
		default-jdk \
		git \
		python3 \
		python3-dev \
		python3-pip

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install setuptools --no-cache-dir && \
    python -m pip install cython --no-cache-dir && \
    python -m pip install leidenalg --no-cache-dir && \
    python -m pip install fitsne --no-cache-dir

RUN apt-get -qq -y autoremove && \
    apt-get clean && \
	rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

COPY . /pegasus/
WORKDIR /pegasus/tests
RUN git clone https://github.com/klarman-cell-observatory/pegasus-test-data.git
WORKDIR /pegasus/
RUN python -m pip install -e .

WORKDIR /pegasus/tests
