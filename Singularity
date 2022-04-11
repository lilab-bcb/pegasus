Bootstrap: docker
From: continuumio/miniconda3:4.10.3

%post
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda install -y -c conda-forge mamba
    mamba install -y -c bioconda pegasuspy>=1.5.0
    mamba install -y -c bioconda rpy2 bioconductor-deseq2 bioconductor-fgsea harmony-pytorch forceatlas2-python scanorama
    mamba install -y -c conda-forge scvi-tools
    pip install nmf-torch
