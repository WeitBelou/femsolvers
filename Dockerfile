FROM quay.io/fenicsproject/stable

WORKDIR /home/fenics/shared
ADD requirements.txt .

USER root
RUN pip install -r requirements.txt

ENTRYPOINT ["python", "src/dispatch.py"]
CMD [""]