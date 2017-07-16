#!/usr/bin/env bash

TEMP=`getopt -o i: --long image: -- "$@"`
eval set -- "${TEMP}"

while true ; do
    case "$1" in
        -i|--image)
            case "$2" in
                "") shift 2 ;;
                *) IMAGE=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

echo Image: ${IMAGE}

fenicsproject run ${IMAGE} "pytest tests"