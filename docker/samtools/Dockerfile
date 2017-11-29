FROM ubuntu:zesty

LABEL maintainer="Ambrose J. Carr <acarr@broadinstitute.org>" \
  software="samtools" \
  version="1.3.1" \
  description="processing sequence alignments in SAM and BAM formats" \
  website="https://samtools.github.io"

RUN apt update && apt install -y \
  samtools=0.1.3-3
