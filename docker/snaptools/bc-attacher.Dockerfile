FROM python:3.6.2

LABEL contributor="Jaeyoung Chun <chunj@mskcc.org>"

# RUN apt update && apt install -y

RUN pip install --upgrade pip && \
    pip3 install biopython==1.73

RUN mkdir /tools
WORKDIR /tools

COPY bc_attacher.py /tools

ENTRYPOINT ["python", "bc_attacher.py"]
