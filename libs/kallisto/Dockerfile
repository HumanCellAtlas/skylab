FROM ubuntu:zesty

LABEL maintainer="Ambrose J. Carr <acarr@broadinstitute.org>" \
  software="kallisto" \
  version="0.43.1" \
  description="kallisto RNA-seq pseudoaligner" \
  website="https://pachterlab.github.io/kallisto"

RUN apt update && apt install -y \
  curl

# download and install the kallisto binary, keep the license
RUN curl -L -o kallisto_linux-v0.43.1.tar.gz \
  https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz && \
  tar -xzf kallisto_linux-v0.43.1.tar.gz && \
  cp kallisto_linux-v0.43.1/kallisto /usr/local/bin/ && \
  cp kallisto_linux-v0.43.1/license.txt . && \
  rm kallisto_linux-v0.43.1.tar.gz && \
  rm -r kallisto_linux-v0.43.1
