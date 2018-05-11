#!/usr/bin/env bash

docker build -t scrtools_graphlayout .
docker tag scrtools_graphlayout regevlab/scrtools_graphlayout
docker push regevlab/scrtools_graphlayout
