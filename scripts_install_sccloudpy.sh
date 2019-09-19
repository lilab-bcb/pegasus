#!/bin/bash

set -e

apt-get update
apt-get install --no-install-recommends -y libfftw3-dev
pip3 install --upgrade pip
pip3 install --upgrade bokeh
pip3 install cmake
pip3 install Cython
pip3 install --upgrade numpy
pip3 install MulticoreTSNE-modified==0.1.post2
pip3 install sccloud[fitsne]==0.14.0

echo "sccloud done"
