FROM python:3.8.1-buster

LABEL maintainer="Lantern Team <lantern@broadinstitute.org>"

COPY requirements.txt .
RUN pip install --upgrade pip
RUN pip3 install numpy
RUN pip3 install -r requirements.txt

RUN mkdir /tools
WORKDIR /tools

COPY create_zarr_optimus.py .
