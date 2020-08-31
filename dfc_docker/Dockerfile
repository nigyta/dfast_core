FROM python:3.8.2

# 環境変数の設定
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# Create working directory
RUN mkdir /work

# Install dependency
RUN pip install biopython && \
    apt-get update && \
    apt install -y default-jre zip prodigal infernal && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install dfast_core
RUN cd / && \
    git clone https://github.com/nigyta/dfast_core && \
    ln -s /dfast_core/dfast /usr/local/bin/ && \
    ln -s /dfast_core/scripts/dfast_file_downloader.py /usr/local/bin/ && \
    ln -s /dfast_core/scripts/reference_util.py /usr/local/bin/

RUN cd /tmp  && \
	curl -LO https://github.com/UCSC-LoweLab/tRNAscan-SE/archive/v2.0.6.tar.gz && \
	tar xfz v2.0.6.tar.gz && \
	cd tRNAscan-SE-2.0.6 && \
	./configure && make && make install && \
    cd .. && \
	rm -r /tmp/tRNAscan-SE-2.0.6 && \
	cd /work

# Prepare reference data (currently disabled)
# RUN dfast_file_downloader.py --protein dfast bifido cyanobase ecoli lab --cdd Cog --hmm TIGR


WORKDIR /work

CMD /bin/bash
