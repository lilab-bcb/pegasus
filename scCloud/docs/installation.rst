Installation
------------

Install **scrtools** locally via Miniconda_::

	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh .
	bash Miniconda3-latest-Linux-x86_64.sh -p /users/foo/miniconda3
	mv Miniconda3-latest-Linux-x86_64.sh /users/foo/miniconda3
	conda create -n scrtools -y pip
	source activate scrtools
	conda install -y -c anaconda numpy
	conda install -y -c anaconda cython 
	conda install -y -c conda-forge fftw
	export CPATH=$CPATH:/users/foo/miniconda3/envs/scrtools/include
	pip install pybind11
	git clone https://github.com/nmslib/hnsw.git hnswlib
	cd hnswlib/python_bindings
	python setup.py install
	cd ../..
	git clone https://github.com/broadinstitute/scRNA-Seq.git scRNA-Seq
	cd scRNA-Seq/scrtools
	pip install .
	cd ../..

Use **scrtools** in UGER
++++++++++++++++++++++++

First, you need to request a RedHat7 server::

	qrsh -q interactive -l h_vmem=4g -l os=RedHat7 -P regevlab

Then, if you have installed **scrtools**, you could activate the virtual environment::

	source activate scrtools

Or, you can use an installed version by typing::

	source /ahg/regevdata/users/libo/miniconda3/bin/activate scrtools

.. _Miniconda: http://conda.pydata.org/miniconda.html
