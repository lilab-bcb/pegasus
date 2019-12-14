Installation
------------

Pegasus works with Python 3 only. **Note:** Pegasus doesn't work with Python 3.8 yet, because Cython has not been updated for this Python version.

Linux
^^^^^

The installation has been tested on Ubuntu (18.04 and 19.04).

Prerequisites
#############

First, install the following dependencies by::

	sudo apt install build-essential cmake libxml2-dev zlib1g-dev wget

Next, you can install Pegasus system-wide by PyPI (see `Linux Install via PyPI`_), or within a Miniconda environment (see `Linux Install via Miniconda`_).

.. _Linux Install via PyPI: ./installation.html#install-via-pypi
.. _Linux Install via Miniconda: ./installation.html#install-via-miniconda

Install via PyPI
################

First, install Python 3 and *pip* tool for Python 3::

	sudo apt install python3 python3-pip

Now install Pegasus via *pip*::

	pip3 install pegasuspy

There are optional packages that you can install:

- **mkl**: This package improves math routines for science and engineering applications::

	pip3 install mkl

- **fitsne**: This package is to calculate t-SNE plots using a faster algorithm FIt-SNE::

	sudo apt install libfftw3-dev
	pip3 install fitsne

- **leiden**: This package provides Leiden clustering algorithm, besides the default Louvain algorithm in Pegasus::

	pip3 install leidenalg

Install via Miniconda
#####################

You can also choose to install pegasus in a dedicated Conda environment without affecting your OS settings.

1. Use the following commands to install a Miniconda on your system::

	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh .
	export CONDA_PATH=/home/foo
	bash Miniconda3-latest-Linux-x86_64.sh -p $CONDA_PATH/miniconda3
	mv Miniconda3-latest-Linux-x86_64.sh $CONDA_PATH/miniconda3
	source ~/.bashrc

Feel free to change ``/home/foo`` to your own directory on handling Miniconda.

In addition, remember to type ``yes`` when asked if you wish the installer to initialize Miniconda3 by running conda init.

2. Create a conda environment for pegasus. This tutorial uses ``pegasus`` as the environment name, but you are free to choose your own::

	conda create -n pegasus -y python=3.7

3. Enter ``pegasus`` environment by activating::

	source activate pegasus

4. (Optional) Install the following dependency if you want *mkl* to for optimized math routines::

	conda install -y -c anaconda numpy

5. Install pegasus::

	pip install pegasuspy

6. (Optional) You can install the following optional features:

	- **fitsne**: Generate t-SNE plot by a faster algorithm FIt-SNE::

		sudo apt install libfftw3-dev
		pip install fitsne

	- **leiden**: Leiden clustering algorithm::

		pip install leidenalg

---------------

macOS
^^^^^

Prerequisites
#############

First, install Homebrew by following the instruction on its website: https://brew.sh/. Then install the following dependencies::

	brew install cmake libxml2 curl libomp

Next, you can install Pegasus system-wide by PyPI (see `macOS Installation via PyPI`_), or within a Miniconda environment (see `macOS Installation via Miniconda`_).

.. _macOS Installation via PyPI: ./installation.html#id2
.. _macOS Installation via Miniconda: ./installation.html#id3

Install via PyPI
################

1. You need to install Python first::

	brew install python3

2. Starting from macOS Mojave (i.e. 10.14), *python-igraph*, one of the dependencies of Pegasus, needs to set the following environment variable before installation::

	export MACOSX_DEPLOYMENT_TARGET=10.14
	pip3 install python-igraph

You should change ``10.14`` to your macOS version number. For example, ``10.15`` is the number for Catalina.

3. Now install Pegasus::

	pip3 install pegasuspy

There are optional packages that you can install:

- **mkl**: This package improves math routines for science and engineering applications::

	pip3 install mkl

- **fitsne**: This package is to calculate t-SNE plots using a faster algorithm FIt-SNE. First, you need to install its dependency *fftw*::

	brew install fftw

Then install *fitsne* by::

	pip3 install fitsne

- **leiden**: This package provides Leiden clustering algorithm, besides the default Louvain algorithm in Pegasus::

	pip3 install leidenalg

Install via Miniconda
#####################

1. Use the following commands to install a Miniconda on your system::

	curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
	export CONDA_PATH=/Users/foo
	bash Miniconda3-latest-MacOSX-x86_64.sh -p $CONDA_PATH/miniconda3
	mv Miniconda3-latest-MacOSX-x86_64.sh $CONDA_PATH/miniconda3

Feel free to change ``/Users/foo`` to your own directory on handling Miniconda.


2. Create a conda environment for pegasus. This tutorial uses ``pegasus`` as the environment name, but you are free to choose your own::

	conda create -n pegasus -y python=3.7

3. Enter ``pegasus`` environment by activating::

	conda activate pegasus

4. (Optional) Install the following dependency if you want *mkl* to for optimized math routines::

	conda install -y -c anaconda numpy

5. **For macOS 10.14 or later:** for these macOS versions, you need to set the following environment variable before installing Pegasus::

	export MACOSX_DEPLOYMENT_TARGET=10.14

where ``10.14`` is the version number for macOS Mojave. You should change it to your own OS version. For example, ``10.15`` is for macOS Catalina.

5. Install pegasus::

	pip install pegasuspy

6. (Optional) You can install the following optional features:

	- **fitsne**: Generate t-SNE plot by a faster algorithm FIt-SNE::

		conda install -y -c conda-forge fftw
		pip install fitsne

	- **leiden**: Leiden clustering algorithm::

		pip install leidenalg

---------------

Development Version
^^^^^^^^^^^^^^^^^^^^^^

To install Pegasus development version directly from `its GitHub respository <https://github.com/klarman-cell-observatory/pegasus>`_, please do the following steps:

1. Install prerequisite libraries as mentioned in above sections. 

2. Install Git. See `here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_ for how to install Git.

3. Use git to fetch repository source code, and install from it::

	git clone https://github.com/klarman-cell-observatory/pegasus.git
	cd pegasus
	pip install -e .

where ``-e`` option of ``pip`` means to install in editing mode, so that your Pegasus installation will be automatically updated upon modifications in source code.