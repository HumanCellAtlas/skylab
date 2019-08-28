FROM python:3.7.2

LABEL maintainer="Lantern Team <lantern@broadinstitute.org>"

COPY requirements.txt .
RUN pip3 install numpy==1.17.0
RUN pip3 install -r requirements.txt

RUN mkdir /tools
WORKDIR /tools

COPY create_zarr_optimus.py .
