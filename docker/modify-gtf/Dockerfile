FROM python:3.7.2

LABEL maintainer="Lantern Team <lantern@broadinstitute.org>" \
    software="modify_gtf" \
    version="0.1.0" \
    description="Tools to modify gtf files."


## Append tools to path
ENV PATH=/tools/:${PATH}

## Copy Scripts
COPY SetGeneNameToId.py /tools/