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
	mkdir -p $CONDA_PATH/software
	git clone https://github.com/klarman-cell-observatory/scCloudPy.git $CONDA_PATH/software/scCloudPy
	cd $CONDA_PATH/software/scCloudPy/scCloud
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
	mkdir -p $CONDA_PATH/software
	export MACOSX_DEPLOYMENT_TARGET=10.9
	git clone https://github.com/klarman-cell-observatory/scCloudPy.git $CONDA_PATH/software/scCloudPy
	cd $CONDA_PATH/software/scCloudPy/scCloud
	pip install -e .
	cd $OLDPWD
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
