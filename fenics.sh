#!/usr/bin/env bash

# Cache for form instant compiler
docker volume create --name instant-cache > /dev/null 2>&1

CUSTOM_IMAGE_NAME="ice_island/custom_fenics_image"
docker build -t ${CUSTOM_IMAGE_NAME} .

docker run -t --rm \
           -v instant-cache:/home/fenics/.instant \
           -v $(pwd):/home/fenics/shared ${CUSTOM_IMAGE_NAME} "$@"