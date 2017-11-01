FROM ubuntu:zesty

LABEL maintainer="Ambrose J. Carr <acarr@broadinstitute.org>" \
  software="STAR" \
  version="2.5.3a-40ead6e" \
  description="RNA-seq aligner, patch of version 2.5.3a to allow uBam alignment" \
  website="https://github.com/alexdobin/STAR"

RUN apt update && apt install -y \
  git \
  samtools

RUN git clone https://github.com/alexdobin/STAR.git && \
  cd STAR && \
  git reset --hard 40ead6efcd9ee2ebf77a510283b69b1a7aad1183 && \
  cp ./bin/Linux_x86_64_static/STAR /usr/local/bin && \
  cp ./LICENSE .. && \
  cd .. && \
  rm -rf STAR
