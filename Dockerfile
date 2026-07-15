FROM vibpsb/i-adhore:3.1

MAINTAINER Cecilia Sensalari

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

# Install Python 3.9 and the last pip release line that supports it.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        software-properties-common && \
    add-apt-repository -y ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y --no-install-recommends python3.9 python3.9-distutils && \
    curl --fail --location --silent --show-error \
        https://bootstrap.pypa.io/pip/3.9/get-pip.py \
        --output /tmp/get-pip.py && \
    python3.9 /tmp/get-pip.py && \
    python3.9 -m pip --version && \
    rm /tmp/get-pip.py && \
    rm -rf /var/lib/apt/lists/*

# Install non-python wgd dependencies
RUN apt-get update && apt-get install -yq \
    wget \
    git \
    curl \
    default-jdk \
    build-essential \
    mcl \
    ncbi-blast+ \
    muscle \
    fasttree

# Install PAML from source
COPY /vendor/paml4.9j.tgz /paml4.9j.tgz
RUN tar -xzf paml4.9j.tgz && cd paml4.9j/src && make -f Makefile && \
    mv codeml /bin && cd /

# Install DIAMOND
RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz && \
    tar -xzf diamond-linux64.tar.gz && mv diamond /bin

# Install OrthoMCLight
RUN wget https://raw.githubusercontent.com/VIB-PSB/OrthoMCLight/main/orthomclight.pl -P /bin && \
    wget https://raw.githubusercontent.com/VIB-PSB/OrthoMCLight/main/orthomclight_module.pm -P /bin && \
    chmod a+rx /usr/bin/orthomclight*

# Copy ksrates files
ADD /requirements.txt /ksrates/requirements.txt
ADD /setup.py /ksrates/setup.py
ADD /ksrates /ksrates/ksrates
ADD /wgd_ksrates /ksrates/wgd_ksrates
ADD /README.md /ksrates/README.md
ADD /ksrates_cli.py /ksrates/ksrates_cli.py

# Download the 37 angiosperm sequence zipped file from Zenodo for the reciprocal retention pipeline
RUN wget https://zenodo.org/records/15225340/files/original_angiosperm_sequences.tar.gz \
    -P /ksrates/ksrates/reciprocal_retention

# Install ksrates and requirements from requirements.txt
RUN python3.9 -m pip install -r /ksrates/requirements.txt && \
    python3.9 -m pip install /ksrates && \
    rm -r /ksrates
