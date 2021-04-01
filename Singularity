Bootstrap: docker
From: vibpsb/i-adhore

%labels
	AUTHOR cecilia.sensalari@psb.vib-ugent.be

%environment
	export LC_ALL=C.UTF-8
	export LANG=C.UTF-8
    export DEBIAN_FRONTEND=noninteractive

%files
	requirements.txt install/requirements.txt
	setup.py install/setup.py
	ksrates install/ksrates
	wgd_ksrates install/wgd_ksrates
	README.md install/README.md
	ksrates_cli.py install/ksrates_cli.py

%post

	# Install PAML from source

	apt-get install -y wget && wget http://abacus.gene.ucl.ac.uk/software/paml4.9j.tgz && \
	tar -xzf paml4.9j.tgz && cd paml4.9j/src && make -f Makefile && mv codeml /bin && cd /

	# Install Python3, wgd dependencies...

	apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -yq install python3-pip python3-tk git curl \
	default-jdk build-essential mcl ncbi-blast+ muscle mafft prank fasttree phyml

	# Install ksrates and requirements from requirements.txt

	python3 -m pip install /install
