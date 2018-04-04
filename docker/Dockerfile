FROM openjdk:8
RUN apt-get update && apt-get install --no-install-recommends -y \
        automake \
        autotools-dev \
	    ca-certificates \
        make \
		zlib1g-dev \
		libbz2-dev \
		liblzma-dev \
		gcc \
		r-base \
		unzip \
		wget

RUN curl -k -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
RUN bash ./miniconda.sh -b -p /miniconda
RUN rm -f ./miniconda.sh
ENV PATH=/miniconda/bin/:$PATH
RUN conda install -y matplotlib numpy pysam
ENV MPLBACKEND="Agg"
RUN mkdir /scripts
COPY scripts/pyadapter_trim.py /scripts/
COPY scripts/pyMakeVplot.py /scripts/
COPY scripts/single_stats.py /scripts/a
RUN wget --quiet -O "/picard.jar" "https://github.com/broadinstitute/picard/releases/download/2.17.0/picard.jar"

RUN wget --quiet -O "/bowtie2.zip" "https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip"
RUN unzip /bowtie2.zip
RUN rm -f /bowtie2.zip

RUN pip install python-Levenshtein MACS2
RUN wget --quiet -O "/samtools-1.6.tar.bz2" "https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2"
RUN tar xf /samtools-1.6.tar.bz2
RUN rm -f /samtools-1.6.tar.bz2
RUN mkdir /samtools
RUN cd /samtools-1.6/ ; ./configure --without-curses --prefix=/samtools/ ; make ; make install
RUN rm -rf /samtools-1.6
ENV PATH=/bowtie2-2.3.4-linux-x86_64:/samtools/bin:${PATH}



