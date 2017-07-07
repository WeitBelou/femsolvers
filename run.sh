#!/usr/bin/env bash


TEMP=`getopt -o i:c: --long image:,config: -- "$@"`
eval set -- "${TEMP}"

while true ; do
    case "$1" in
        -i|--image)
            case "$2" in
                "") shift 2 ;;
                *) IMAGE=$2 ; shift 2 ;;
            esac ;;
        -c|--config)
            case "$2" in
                "") shift 2 ;;
                *) CONFIG=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

echo Image: ${IMAGE}
echo Config: ${CONFIG}

fenicsproject run ${IMAGE} "python3 src/runner.py ${CONFIG}"