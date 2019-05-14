FROM python:2.7

LABEL contributor="Jaeyoung Chun <chunj@mskcc.org>" \
    software="Barcode Attacher for scATAC-seq" \
    version="0.0.1" \
    description="Barcode Attacher for scATAC-seq"

RUN pip install --upgrade pip && \
    pip install biopython==1.73

RUN mkdir /tools
WORKDIR /tools

COPY bc_attacher.py /tools

# ENTRYPOINT ["python", "bc_attacher.py"]
