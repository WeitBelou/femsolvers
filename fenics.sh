#!/usr/bin/env bash

# Cache for form instant compiler
docker volume create --name instant-cache > /dev/null 2>&1

docker run -t --rm \
           -v instant-cache:/home/fenics/.instant \
           -v $(pwd):/home/fenics/shared \
           -w /home/fenics/shared quay.io/fenicsproject/stable "$@"