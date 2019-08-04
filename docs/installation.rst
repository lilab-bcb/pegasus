Installation
------------

Prerequisites
^^^^^^^^^^^^^

Platform
########

This installation instruction has been tested on Linux (Ubuntu Linux 18.04 and 19.04) and macOS (macOS 10.14 Mojave).

Python
######

scCloud works with Python 3 (Python 3.5 or greater) only. 

For Ubuntu, install Python 3 by::

	sudo apt install python

For macOS, Python 3 can be installed by either Miniconda package manger (see below_), Homebrew_, or the `Python website`_.

Others
######

For Ubuntu, install the following packages::

	sudo apt install build-essential cmake libxml2-dev wget

For macOS, install by::

	brew install cmake libxml2 curl


Install via PyPI
^^^^^^^^^^^^^^^^

Linux
#####

Assuming Python 3 is installed on your OS. Then install ``pip`` for Python 3, along with some extra libraries for building scCloud::

	sudo apt install python3-pip


Now install scCloud via ``pip``::

	pip3 install scCloud

Alternatively, if you want to use ``mkl`` package for speed improvement, type::

	pip3 install scCloud[mkl]

If you want to use scCloud's FIt-SNE feature. First, install ``fftw`` library::

	sudo apt install libfftw3-dev

Then type::

	pip3 install scCloud[fitsne]

Or if you want both extra features, type::

	pip3 install scCloud[mkl,fitsne]

.. _below: ./installation.html#install-via-miniconda

.. _Homebrew: https://brew.sh

.. _Python website: https://www.python.org/downloads/mac-osx/


macOS
######

``pip3`` is already installed. Follow the same steps as above for Linux, except that when installing ``fftw``, you have to go to its website_ for instruction.

.. _website: http://www.fftw.org/

------------------------

Install via Miniconda
^^^^^^^^^^^^^^^^^^^^^

You can also choose to install scCloud in a dedicated Conda environment without affecting your OS.


Package Manager
###############

Linux
*****

Use the following commands to install a Miniconda on your system::

	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh .
	export CONDA_PATH=/home/foo
	bash Miniconda3-latest-Linux-x86_64.sh -p $CONDA_PATH/miniconda3
	mv Miniconda3-latest-Linux-x86_64.sh $CONDA_PATH/miniconda3
	source ~/.bashrc

Feel free to change ``/home/foo`` to your own directory on handling Miniconda.

macOS
*****

Use the following commands::

	curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
	export CONDA_PATH=/Users/foo
	bash Miniconda3-latest-MacOSX-x86_64.sh -p $CONDA_PATH/miniconda3
	mv Miniconda3-latest-MacOSX-x86_64.sh $CONDA_PATH/miniconda3

Feel free to change ``/Users/foo`` to your own directory on handling Miniconda.

Install scCloud
###############

Both Linux and macOS share this installation step.

1. Create a conda environment for scCloud. Let ``scCloud`` be its name, but you are free to choose your own::

	conda create -n scCloud -y pip
For macOS, you may need to do the following extra commands after creation::

	mkdir -p $CONDA_PATH/miniconda3/envs/scCloud/etc/conda/activate.d
	mkdir -p $CONDA_PATH/miniconda3/envs/scCloud/etc/conda/deactivate.d
	printf '#!/bin/sh\n\nexport KMP_DUPLICATE_LIB_OK=true\n' > $CONDA_PATH/miniconda3/envs/scCloud/etc/conda/activate.d/env_vars.sh
	printf '#!/bin/sh\n\nunset KMP_DUPLICATE_LIB_OK' > $CONDA_PATH/miniconda3/envs/scCloud/etc/conda/deactivate.d/env_vars.sh
where ``$CONDA_PATH`` is set in the previous step.

2. Enter conda environment by activating::

	conda activate scCloud
or::

	source activate scCloud

3. (Optional) If you want to use the Intel ``mkl`` package for speed improvement, type::

	conda install -y -c anaconda numpy
Also, if you want to use scCloud's FIt-SNE feature, which depends on ``fftw`` package, type::

	conda install -y -c conda-forge fftw

4. Install scCloud::

	pip install scCloud
If you want to use scCloud's FIt-SNE feature, type::

	pip install scCloud[fitsne]

-----------------------------------

Use **scCloud** in UGER
^^^^^^^^^^^^^^^^^^^^^^^

First, you need to request a RedHat7 server::

	qrsh -q interactive -l h_vmem=4g -l os=RedHat7 -P regevlab

Then, if you have installed **scCloud**, you could activate the virtual environment::

	source activate scCloud

Or, you can use an installed version by typing::

	source /ahg/regevdata/users/libo/miniconda3/bin/activate scCloud

.. _Miniconda: http://conda.pydata.org/miniconda.html
