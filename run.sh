#!/usr/bin/env bash

docker run --rm \
    --env HOST_UID=$(id -u) \
    --env HOST_GID=$(id -g) \
    -v $(pwd):/home/fenics/shared \
    -w /home/fenics/shared \
    quay.io/fenicsproject/stable \
    "python3 -m femsolvers"
