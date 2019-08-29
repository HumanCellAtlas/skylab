FROM python:3.7.4

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY optimus_zarr_to_loom.py .
COPY unpackZARR.sh .