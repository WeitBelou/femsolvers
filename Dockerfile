FROM quay.io/fenicsproject/stable

WORKDIR /home/fenics/shared
ADD . .

USER root
RUN pip install -r requirements.txt

ENTRYPOINT ["python", "src/main.py"]
CMD [""]