Installation
------------

Linux
+++++
This installation instruction has been tested on Ubuntu Linux 18.04 and 19.04.

Suppose your Linux user directory is ``/home/foo``. We will create two folders ``/home/foo/miniconda3`` and ``/home/foo/software``.

Please use the commands below to install **scCloud** locally via Miniconda_::

	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh .
	export CONDA_PATH=/home/foo
	bash Miniconda3-latest-Linux-x86_64.sh -p $CONDA_PATH/miniconda3
	mv Miniconda3-latest-Linux-x86_64.sh $CONDA_PATH/miniconda3
	source ~/.bashrc
	conda create -n scCloud -y pip
	conda activate scCloud
	conda install -y -c anaconda numpy
	conda install -y -c anaconda cython
	conda install -y -c anaconda pybind11 
	conda install -y -c conda-forge fftw
	conda install -y -c anaconda pytables
	conda install -y -c anaconda cmake
	conda install -y -c anaconda libxml2
	export CPATH=$CPATH:$CONDA_PATH/miniconda3/envs/scCloud/include
	mkdir -p $CONDA_PATH/software
	git clone https://github.com/nmslib/hnsw.git $CONDA_PATH/software/hnswlib
	cd $CONDA_PATH/software/hnswlib/python_bindings
	python setup.py install
	cd $OLDPWD
	git clone https://github.com/bli25broad/fishers_exact_test.git $CONDA_PATH/software/fisher_test
	cd $CONDA_PATH/software/fisher_test
	pip install .
	cd $OLDPWD
	git clone https://github.com/lilab-cbb/anndata.git $CONDA_PATH/software/anndata
	cd $CONDA_PATH/software/anndata
	pip install .
	cd $OLDPWD
	git clone https://github.com/bli25broad/louvain-igraph.git $CONDA_PATH/software/louvain
	cd $CONDA_PATH/software/louvain
	pip install .
	cd $OLDPWD
	git clone https://github.com/bli25broad/Multicore-TSNE.git $CONDA_PATH/software/MulticoreTSNE
	cd $CONDA_PATH/software/MulticoreTSNE
	pip install .
	cd $OLDPWD
	git clone https://github.com/bli25broad/umap.git $CONDA_PATH/software/umap
	cd $CONDA_PATH/software/umap
	pip install .
	cd $OLDPWD
	git clone https://github.com/broadinstitute/scRNA-Seq.git $CONDA_PATH/software/scRNA-Seq
	cd $CONDA_PATH/software/scRNA-Seq/scCloud
	pip install -e .
	cd $OLDPWD

Mac OS
++++++

Suppose your MacOS user directory is ``/Users/foo``. We will create two folders ``/Users/foo/miniconda3`` and ``/Users/foo/software``.

Please use the commands below to install **scCloud** locally via Miniconda_::

	curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
	export CONDA_PATH=/Users/foo
	bash Miniconda3-latest-MacOSX-x86_64.sh -p $CONDA_PATH/miniconda3
	mv Miniconda3-latest-MacOSX-x86_64.sh $CONDA_PATH/miniconda3
	conda create -n scCloud -y pip
	mkdir -p $CONDA_PATH/miniconda3/envs/scCloud/etc/conda/activate.d
	mkdir -p $CONDA_PATH/miniconda3/envs/scCloud/etc/conda/deactivate.d
	printf '#!/bin/sh\n\nexport KMP_DUPLICATE_LIB_OK=true\n' > $CONDA_PATH/miniconda3/envs/scCloud/etc/conda/activate.d/env_vars.sh
	printf '#!/bin/sh\n\nunset KMP_DUPLICATE_LIB_OK' > $CONDA_PATH/miniconda3/envs/scCloud/etc/conda/deactivate.d/env_vars.sh
	conda activate scCloud
	conda install -y -c anaconda numpy
	conda install -y -c anaconda cython
	conda install -y -c anaconda pybind11
	conda install -y -c anaconda lightgbm
	conda install -y -c anaconda cmake
	conda install -y -c anaconda pytables
	conda install -y -c conda-forge fftw
	export CPATH=$CPATH:$CONDA_PATH/miniconda3/envs/scCloud/include
	mkdir -p $CONDA_PATH/software
	git clone https://github.com/nmslib/hnsw.git $CONDA_PATH/software/hnswlib
	cd $CONDA_PATH/software/hnswlib/python_bindings
	python setup.py install
	cd $OLDPWD
	git clone https://github.com/bli25broad/fishers_exact_test.git $CONDA_PATH/software/fisher_test
	cd $CONDA_PATH/software/fisher_test
	pip install .
	cd $OLDPWD
	git clone https://github.com/lilab-cbb/anndata.git $CONDA_PATH/software/anndata
	cd $CONDA_PATH/software/anndata
	pip install .
	cd $OLDPWD
	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
	brew install gcc
	sudo ln -s /usr/local/bin/gcc-9 /usr/local/bin/gcc
	sudo ln -s /usr/local/bin/g++-9 /usr/local/bin/g++
	export PATH=$CONDA_PATH/miniconda3/envs/scCloud/bin:/usr/local/bin:$PATH
	git clone https://github.com/bli25broad/louvain-igraph.git $CONDA_PATH/software/louvain
	cd $CONDA_PATH/software/louvain
	pip install .
	cd $OLDPWD
	git clone https://github.com/bli25broad/Multicore-TSNE.git $CONDA_PATH/software/MulticoreTSNE
	cd $CONDA_PATH/software/MulticoreTSNE
	pip install .
	cd $OLDPWD
	git clone https://github.com/bli25broad/umap.git $CONDA_PATH/software/umap
	cd $CONDA_PATH/software/umap
	pip install .
	cd $OLDPWD
	git clone https://github.com/broadinstitute/scRNA-Seq.git $CONDA_PATH/software/scRNA-Seq
	cd $CONDA_PATH/software/scRNA-Seq/scCloud
	pip install -e .
	cd $OLDPWD

Use **scCloud** in UGER
++++++++++++++++++++++++

First, you need to request a RedHat7 server::

	qrsh -q interactive -l h_vmem=4g -l os=RedHat7 -P regevlab

Then, if you have installed **scCloud**, you could activate the virtual environment::

	source activate scCloud

Or, you can use an installed version by typing::

	source /ahg/regevdata/users/libo/miniconda3/bin/activate scCloud

.. _Miniconda: http://conda.pydata.org/miniconda.html
