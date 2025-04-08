FROM vibpsb/i-adhore:3.1

MAINTAINER Cecilia Sensalari

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

# Install Python3.8
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y python3.8 python3.8-distutils && \
	curl -sS https://bootstrap.pypa.io/get-pip.py | python3.8 && \
	python3.8 -m pip install --upgrade pip

# Install non-python wgd dependencies
RUN apt-get install -yq git curl default-jdk build-essential mcl ncbi-blast+ muscle fasttree 

# Install PAML from source
RUN apt-get install -y wget && wget https://gitlab.renkulab.io/tzzteresa/polyploid-summer-school-2023day4-2/-/raw/rimjhim.choudhury-master-patch-51949/src/paml4.9j.tgz && \
	tar -xzf paml4.9j.tgz && cd paml4.9j/src && make -f Makefile && mv codeml /bin && cd /

# Install DIAMOND
RUN	wget http://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz && \
	tar -xzf diamond-linux64.tar.gz && mv diamond /bin

# Install OrthoMCLight
RUN wget https://raw.githubusercontent.com/VIB-PSB/OrthoMCLight/main/orthomclight.pl -P /bin && \
	wget https://raw.githubusercontent.com/VIB-PSB/OrthoMCLight/main/orthomclight_module.pm -P /bin && \
	chmod a+rx /usr/bin/orthomclight*

# Copy ksrates files
ADD /requirements.txt /install/requirements.txt
ADD /setup.py /install/setup.py
ADD /ksrates /install/ksrates
ADD /wgd_ksrates /install/wgd_ksrates
ADD /README.md /install/README.md
ADD /ksrates_cli.py /install/ksrates_cli.py

# Install ksrates and requirements from requirements.txt
RUN python3.8 -m pip install /install && \
	rm -r /install
