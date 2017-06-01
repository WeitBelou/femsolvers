FROM quay.io/fenicsproject/stable

WORKDIR /home/fenics/shared
ADD . ${WORKDIR}

USER root
RUN pip install -r requirements.txt
USER root