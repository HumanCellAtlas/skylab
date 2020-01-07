FROM python:3.7.4

COPY requirements.txt .
RUN pip install -r requirements.txt

RUN mkdir /tools
COPY optimus_zarr_to_loom.py /tools/
COPY unpackZARR.sh /tools/
COPY ss2_plate_zarr_to_loom.py /tools/

ENV PATH=${PATH}:/tools/
