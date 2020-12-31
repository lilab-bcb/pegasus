Bootstrap: docker
From: ubuntu:20.04

%post
    apt update
    apt install -y --no-install-recommends build-essential libfftw3-dev default-jdk python3 python3-pip python3-dev
    ln -s /usr/bin/python3 /usr/bin/python
    python -m pip install --upgrade pip
    python -m pip install pegasuspy
    python -m pip install fitsne leidenalg
 