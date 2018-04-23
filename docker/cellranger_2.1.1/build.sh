#!/usr/bin/env bash

# Download cellranger-2.1.1.tar into this directory from cell ranger website before building
docker build -t cellranger-2.1.1 .
docker tag cellranger-2.1.1 regevlab/cellranger-2.1.1
docker push regevlab/cellranger-2.1.1
