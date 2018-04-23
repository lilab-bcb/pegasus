#!/usr/bin/env bash

docker build -t scrtools .
docker tag scrtools regevlab/scrtools
docker push regevlab/scrtools
