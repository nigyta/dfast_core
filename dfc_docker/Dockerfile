# Pin to linux/amd64: the bundled binaries under bin/Linux are x86_64,
# so the image must be amd64 (also keeps the libidn.so path below valid).
FROM --platform=linux/amd64 python:3.13

# 環境変数の設定
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV INFERNAL_NCPU 1

# Create working directory
RUN mkdir /work && chmod 777 /work 

# Install dependency
RUN pip install biopython && \
    apt-get update && \
    apt install -y default-jre zip prodigal infernal ncbi-blast+ && \
    ln -s /usr/bin/cmscan /usr/local/bin/cmscan && \
    ln -s /usr/bin/cmsearch /usr/local/bin/cmsearch && \
    apt-get clean && rm -rf /var/lib/apt/lists/* 

RUN cd /tmp  && \
    curl -LO https://github.com/UCSC-LoweLab/tRNAscan-SE/archive/v2.0.12.tar.gz && \
    tar xfz v2.0.12.tar.gz && \
    cd tRNAscan-SE-2.0.12 && \
    ./configure && make && make install && \
    cd .. && \
    rm -r /tmp/tRNAscan-SE-2.0.12 /tmp/v2.0.12.tar.gz && \
    cd /work

# For compatibility (required for older version of BLAST binaries that link libidn.so.11).
# Debian 13 (trixie) only ships libidn12, so install it and provide a libidn.so.11 symlink.
RUN apt-get update && apt-get install -y libidn12 && \
    ln -sf /usr/lib/x86_64-linux-gnu/libidn.so.12 /usr/lib/x86_64-linux-gnu/libidn.so.11 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Prepare reference data (currently disabled)
# RUN dfast_file_downloader.py --protein dfast bifido cyanobase ecoli lab --cdd Cog --hmm TIGR

# PlasmidFinder (v3.x) and KMA
RUN pip install plasmidfinder && \
    printf '#!/bin/bash\nexec python3 -m plasmidfinder "$@"\n' > /usr/local/bin/plasmidfinder.py && \
    chmod +x /usr/local/bin/plasmidfinder.py

# MobileElementFinder (installed into the same env; biopython pin relaxed)
RUN cd /opt && git clone https://bitbucket.org/mhkj/mge_finder.git && \
    sed -i 's/biopython<=1.80/biopython/' mge_finder/setup.cfg && \
    pip install ./mge_finder "setuptools<81" && \
    rm -rf mge_finder

# KMA is used for PlasmidFinder, but blastn will be used instead of KMA for the FASTA-formated nucleotide sequence search.
# RUN cd / && \
#     git clone --branch 1.0.1 --depth 1 https://bitbucket.org/genomicepidemiology/kma.git && \
#     cd kma && make && \
#     mv kma* /usr/local/bin/ && \
#     cd / && rm -r /kma

WORKDIR /work

CMD /bin/bash

ENV INCREMENT_THIS_TO_BUMP 1
# Install dfast_core
RUN cd / && \
    git clone https://github.com/nigyta/dfast_core && \
    ln -s /dfast_core/dfast /usr/local/bin/ && \
    ln -s /dfast_core/scripts/dfast_file_downloader.py /usr/local/bin/ && \
    ln -s /dfast_core/scripts/reference_util.py /usr/local/bin/ && \
    ln -s /dfast_core/scripts/reference_util_for_nucl.py /usr/local/bin/

# Install DFAST Record tools
RUN pip install "git+https://github.com/ddbj/dr_tools.git"
