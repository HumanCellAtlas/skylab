FROM python:3.6.2

LABEL maintainer="Nick Barkas <nbarkas@broadinstitute.org>" \
    software="umi_tools" \
    version="0.5.5" \
    description="tools for extraction correction, deduplication and counting of UMIs" \
    website="https://github.com/CGATOxford/UMI-tools"

RUN git clone https://github.com/CGATOxford/UMI-tools.git
WORKDIR UMI-tools
RUN git checkout tags/0.5.5
RUN pip install .

