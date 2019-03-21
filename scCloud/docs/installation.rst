Installation
------------

Linux
+++++
This installation instruction has been tested on Ubuntu Linux 18.04.

Suppose your Linux user directory is ``/users/foo``. We will create two folders ``/users/foo/miniconda3`` and ``/users/foo/software``.

Please use the commands below to install **scCloud** locally via Miniconda_::

	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh .
	bash Miniconda3-latest-Linux-x86_64.sh -p /users/foo/miniconda3
	mv Miniconda3-latest-Linux-x86_64.sh /users/foo/miniconda3
	source ~/.bashrc
	conda create -n scCloud -y pip
	conda activate scCloud
	conda install -y -c anaconda numpy
	conda install -y -c anaconda cython
	conda install -y -c anaconda pybind11 
	conda install -y -c conda-forge fftw
	conda install -y -c anaconda pytables
	export CPATH=$CPATH:/users/foo/miniconda3/envs/scCloud/include
	mkdir -p /users/foo/software
	git clone https://github.com/nmslib/hnsw.git /users/foo/software/hnswlib
	cd /users/foo/software/hnswlib/python_bindings
	python setup.py install
	cd $OLDPWD
	git clone https://github.com/bli25broad/fishers_exact_test.git /users/foo/software/fisher_test
	cd /users/foo/software/fisher_test
	pip install .
	cd $OLDPWD
	git clone https://github.com/broadinstitute/scRNA-Seq.git /users/foo/software/scRNA-Seq
	cd /users/foo/software/scRNA-Seq/scCloud
	pip install -e .
	cd $OLDPWD

Mac OS
++++++

Suppose your Linux user directory is ``/Users/foo``. We will create two folders ``/Users/foo/miniconda3`` and ``/Users/foo/software``.

Please use the commands below to install **scCloud** locally via Miniconda_::

	curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -p /Users/foo/miniconda3
	mv Miniconda3-latest-Linux-x86_64.sh /Users/foo/miniconda3
	conda create -n scCloud -y pip
	mkdir -p /Users/foo/miniconda3/scCloud/etc/conda/activate.d
	mkdir -p /Users/foo/miniconda3/scCloud/etc/conda/deactivate.d
	printf '#!/bin/sh\n\nexport KMP_DUPLICATE_LIB_OK=true\n' > /Users/foo/miniconda3/scCloud/etc/conda/activate.d/env_vars.sh
	printf '#!/bin/sh\n\nunset KMP_DUPLICATE_LIB_OK' > /Users/foo/miniconda3/scCloud/etc/conda/deactivate.d/env_vars.sh
	conda activate scCloud
	conda install -y -c anaconda numpy
	conda install -y -c anaconda cython
	conda install -y -c anaconda pybind11
	conda install -y -c anaconda lightgbm
	conda install -y -c anaconda cmake
	conda install -y -c anaconda pytables
	conda install -y -c conda-forge fftw
	export CPATH=$CPATH:/Users/foo/miniconda3/envs/scCloud/include
	mkdir -p /Users/foo/software
	git clone https://github.com/nmslib/hnsw.git /Users/foo/software/hnswlib
	cd /Users/foo/software/hnswlib/python_bindings
	python setup.py install
	cd $OLDPWD
	git clone https://github.com/bli25broad/fishers_exact_test.git /Users/foo/software/fisher_test
	cd /Users/foo/software/fisher_test
	pip install .
	cd $OLDPWD
	git clone https://github.com/bli25broad/louvain-igraph.git /Users/foo/software/louvain
	cd /Users/foo/software/louvain
	pip install .
	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
	brew install gcc
	sudo ln -s /usr/local/bin/gcc-8 /usr/local/bin/gcc
	sudo ln -s /usr/local/bin/g++-8 /usr/local/bin/g++
	export PATH=/Users/foo/miniconda3/envs/scCloud/bin:/usr/local/bin:$PATH
	git clone https://github.com/broadinstitute/scRNA-Seq.git /Users/foo/software/scRNA-Seq
	cd /Users/foo/software/scRNA-Seq/scCloud
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
