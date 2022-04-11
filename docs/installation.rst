Installation
------------

Pegasus works with Python 3.7, 3.8 and 3.9.

Linux
^^^^^

Ubuntu/Debian
###############

Prerequisites
+++++++++++++++

On Ubuntu/Debian Linux, first install the following dependency by::

	sudo apt install build-essential

Next, you can install Pegasus system-wide by PyPI (see `Ubuntu/Debian install via PyPI`_), or within a Miniconda environment (see `Install via Conda`_).

To use the Force-directed-layout (FLE) embedding feature, you'll need Java. You can either install `Oracle JDK`_, or install OpenJDK which is included in Ubuntu official repository::

	sudo apt install default-jdk

Ubuntu/Debian install via PyPI
+++++++++++++++++++++++++++++++++

First, install Python 3, *pip* tool for Python 3 and *Cython* package::

	sudo apt install python3 python3-pip
	python3 -m pip install --upgrade pip
	python3 -m pip install cython

Now install Pegasus with the required dependencies via *pip*::

	python3 -m pip install pegasuspy

or install Pegasus with all dependencies::

	python3 -m pip install pegasuspy[all]

Alternatively, you can install Pegasus with some of the additional optional dependencies as below:

- **torch**: This includes ``harmony-pytorch`` for data integration and ``nmf-torch`` for NMF and iNMF data integration, both of which uses PyTorch::

	python3 -m pip install pegasuspy[torch]

- **louvain**: This includes ``louvain`` package, which provides Louvain clustering algorithm, besides the default Leiden algorithm in Pegasus::

	python3 -m pip install pegasuspy[louvain]

.. note::

	If installing from Python 3.9, to install ``louvain``, you'll need to install the following packages system-wide first in order to locally compile it::

		sudo apt install flex bison libtool

- **tsne**: This package is to calculate t-SNE plots using a fast algorithm FIt-SNE::

	sudo apt install libfftw3-dev
	python3 -m pip install pegasuspy[tsne]

- **forceatlas**: This includes ``forceatlas2-python`` package, a multi-thread Force Atlas 2 implementation for trajectory analysis::

	python3 -m pip install pegasuspy[forceatlas]

- **scanorama**: This includes ``scanorama`` package, a widely-used method for batch correction::

	python3 -m pip install pegasuspy[scanorama]

- **mkl**: This includes ``mkl`` package, which improve math routines for science and engineering applications. Notice that *mkl* not included in **pegasuspy[all]** above::

	python3 -m pip install pegasuspy[mkl]

- **rpy2**: This includes ``rpy2`` package, which is used by Pegasus wrapper on R functions, such as ``fgsea`` and ``DESeq2``::

	python3 -m pip install pegasuspy[rpy2]

- **scvi**: This includes ``scvi-tools`` package for data integration::

	python3 -m pip install pegasuspy[scvi]



--------------------------

Fedora
########

Prerequisites
++++++++++++++

On Fedora Linux, first install the following dependency by::

	sudo dnf install gcc gcc-c++

Next, you can install Pegasus system-wide by PyPI (see `Fedora install via PyPI`_), or within a Miniconda environment (see `Install via Conda`_).

To use the Force-directed-layout (FLE) embedding feature, you'll need Java. You can either install `Oracle JDK`_, or install OpenJDK which is included in Fedora official repository (e.g. ``java-latest-openjdk``)::

	sudo dnf install java-latest-openjdk

or other OpenJDK version chosen from the searching result of command::

	dnf search openjdk

Fedora install via PyPI
+++++++++++++++++++++++++

We'll use Python 3.8 in this tutorial.

First, install Python 3 and *pip* tool for Python 3::

	sudo dnf install python3.8
	python3.8 -m ensurepip --user
	python3.8 -m pip install --upgrade pip

Now install Pegasus with the required dependencies via *pip*::

	python3.8 -m pip install pegasuspy

or install Pegasus with all dependencies::

	python3.8 -m pip install pegasuspy[all]

Alternatively, you can install Pegasus with some of the additional optional dependencies as below:

- **torch**: This includes ``harmony-pytorch`` for data integration and ``nmf-torch`` for NMF and iNMF data integration, both of which uses PyTorch::

	python3.8 -m pip install pegasuspy[torch]

- **louvain**: This includes ``louvain`` package, which provides Louvain clustering algorithm, besides the default Leiden algorithm in Pegasus::

	python3.8 -m pip install pegasuspy[louvain]

.. note::

	If installing from Python 3.9, to install ``louvain``, you'll need to install the following packages system-wide first in order to locally compile it::

		sudo dnf install flex bison libtool

- **tsne**: This package is to calculate t-SNE plots using a fast algorithm FIt-SNE::

	sudo apt install libfftw3-dev
	python3.8 -m pip install pegasuspy[tsne]

- **forceatlas**: This includes ``forceatlas2-python`` package, a multi-thread Force Atlas 2 implementation for trajectory analysis::

	python3.8 -m pip install pegasuspy[forceatlas]

- **scanorama**: This includes ``scanorama`` package, a widely-used method for batch correction::

	python3.8 -m pip install pegasuspy[scanorama]

- **mkl**: This includes ``mkl`` package, which improve math routines for science and engineering applications. Notice that *mkl* not included in **pegasuspy[all]** above::

	python3.8 -m pip install pegasuspy[mkl]

- **rpy2**: This includes ``rpy2`` package, which is used by Pegasus wrapper on R functions, such as ``fgsea`` and ``DESeq2``::

	python3.8 -m pip install pegasuspy[rpy2]

- **scvi**: This includes ``scvi-tools`` package for data integration::

	python3.8 -m pip install pegasuspy[scvi]


.. _Ubuntu/Debian install via PyPI: ./installation.html#ubuntu-debian-install-via-pypi
.. _Fedora install via PyPI: ./installation.html#fedora-install-via-pypi

---------------

macOS
^^^^^

Prerequisites
#############

First, install Homebrew by following the instruction on its website: https://brew.sh/. Then install the following dependencies::

	brew install libomp

And install macOS command line tools::

	xcode-select --install

Next, you can install Pegasus system-wide by PyPI (see `macOS installation via PyPI`_), or within a Miniconda environment (see `Install via Conda`_).

To use the Force-directed-layout (FLE) embedding feature, you'll need Java. You can either install `Oracle JDK`_, or install OpenJDK via Homebrew::

	brew install java

.. _macOS installation via PyPI: ./installation.html#macos-install-via-pypi

macOS install via PyPI
#######################

1. You need to install Python and *pip* tool first::

	brew install python3
	python3 -m pip install --upgrade pip

2. Now install Pegasus with required dependencies via *pip*::

	python3 -m pip install pegasuspy

or install Pegasus with all dependencies::

	python3 -m pip install pegasuspy[all]

Alternatively, you can install Pegasus with some of the additional optional dependencies as below:

- **torch**: This includes ``harmony-pytorch`` for data integration and ``nmf-torch`` for NMF and iNMF data integration, both of which uses PyTorch::

	python3 -m pip install pegasuspy[torch]

- **louvain**: This includes ``louvain`` package, which provides Louvain clustering algorithm, besides the default Leiden algorithm in Pegasus::

	python3 -m pip install pegasuspy[louvain]

- **tsne**: This package is to calculate t-SNE plots using a fast algorithm FIt-SNE::

	sudo apt install libfftw3-dev
	python3 -m pip install pegasuspy[tsne]

- **forceatlas**: This includes ``forceatlas2-python`` package, a multi-thread Force Atlas 2 implementation for trajectory analysis::

	python3 -m pip install pegasuspy[forceatlas]

- **scanorama**: This includes ``scanorama`` package, a widely-used method for batch correction::

	python3 -m pip install pegasuspy[scanorama]

- **mkl**: This includes ``mkl`` packages, which improve math routines for science and engineering applications. Notice that *mkl* not included in **pegasuspy[all]** above::

	python3 -m pip install pegasuspy[mkl]

- **rpy2**: This includes ``rpy2`` package, which is used by Pegasus wrapper on R functions, such as ``fgsea`` and ``DESeq2``::

	python3 -m pip install pegasuspy[rpy2]

- **scvi**: This includes ``scvi-tools`` package for data integration::

	python3 -m pip install pegasuspy[scvi]


----------------------

Install via Conda
^^^^^^^^^^^^^^^^^^

Alternatively, you can install Pegasus via Conda, which is a separate virtual environment without touching your system-wide packages and settings.

You can install Anaconda_, or Miniconda_ (a minimal installer of conda). In this tutorial, we'll use Miniconda.

1. Download `Miniconda installer`_ for your OS. For example, if on 64-bit Linux, then use the following commands to install Miniconda::

	export CONDA_PATH=/home/foo
	bash Miniconda3-latest-Linux-x86_64.sh -p $CONDA_PATH/miniconda3
	mv Miniconda3-latest-Linux-x86_64.sh $CONDA_PATH/miniconda3
	source ~/.bashrc

where ``/home/foo`` should be replaced by the directory to which you want to install Miniconda. Similarly for macOS.

2. Create a conda environment for pegasus. This tutorial uses ``pegasus`` as the environment name, but you are free to choose your own::

	conda create -n pegasus -y python=3.8

Also notice that Python ``3.8`` is used in this tutorial. To choose a different version of Python, simply change the version number in the command above.

3. Enter ``pegasus`` environment by activating::

	conda activate pegasus

4. Install Pegasus via conda::

	conda install -y -c bioconda pegasuspy

5. (Optional) Use the following command to add support on ``harmony-pytorch`` and ``nmf-torch``::

	conda install -y -c bioconda harmony-pytorch
	pip install nmf-torch

Enalbe Force Atlas 2 for trajectory analysis::

	conda install -y -c bioconda forceatlas2-python

Enable support on ``scanorama``::

	conda install -y -c bioconda scanorama

Enable support on ``fgsea`` and ``deseq2`` packages::

	conda install -y -c bioconda rpy2 bioconductor-fgsea bioconductor-deseq2

Enable support on ``scvi-tools``::

	conda install -y -c conda-forge scvi-tools

--------------------------

Install via Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^

Singularity_ is a container engine similar to Docker. Its main difference from Docker is that Singularity can be used with unprivileged permissions.

.. note::

	Please notice that Singularity Hub has been offline since April 26th, 2021 (see `blog post`_). All existing containers held there are in archive, and we can no longer push new builds.

	So if you fetch the container from Singularity Hub using the following command::

		singularity pull shub://klarman-cell-observatory/pegasus

	it will just give you a Singularity container of Pegasus v1.2.0 running on Ubuntu Linux 20.04 base with Python 3.8, in the name ``pegasus_latest.sif`` of about 2.4 GB.

On your local machine, first `install Singularity`_, then you can use our `Singularity spec file`_ to build a Singularity container by yourself::

	singularity build pegasus.sif Singularity

After that, you can interact with it by running the following command::

	singularity run pegasus.sif

Please refer to `Singularity image interaction guide`_ for details.


--------------------------

Development Version
^^^^^^^^^^^^^^^^^^^^^^

To install Pegasus development version directly from `its GitHub respository <https://github.com/lilab-bcb/pegasus>`_, please do the following steps:

1. Install prerequisite libraries as mentioned in above sections.

2. Install Git. See `here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_ for how to install Git.

3. Use git to fetch repository source code, and install from it::

	git clone https://github.com/lilab-bcb/pegasus.git
	cd pegasus
	pip install -e .[all]

where ``-e`` option of ``pip`` means to install in editing mode, so that your Pegasus installation will be automatically updated upon modifications in source code.


.. _Oracle JDK: https://www.oracle.com/java/
.. _Anaconda: https://www.anaconda.com/products/individual#Downloads
.. _Miniconda: https://docs.conda.io/en/latest/index.html
.. _Miniconda installer: https://docs.conda.io/en/latest/miniconda.html
.. _Singularity: http://singularity.lbl.gov/
.. _blog post: https://vsoch.github.io//2021/singularity-hub-archive/
.. _install Singularity: https://sylabs.io/guides/3.8/user-guide/quick_start.html#quick-installation-steps
.. _Singularity spec file: https://raw.githubusercontent.com/lilab-bcb/pegasus/master/Singularity
.. _Singularity image interaction guide: https://singularityhub.github.io/singularityhub-docs/docs/interact
