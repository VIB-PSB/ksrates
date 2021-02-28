Installation
************

*ksrates* can be either installed locally on the machine that is going to pipeline (local computer or a remote compute cluster) or can be shipped in a Singularity and Docker container. Such containers are an isolated portable environment in which the package and the relative dependencies are already installed, which comes in handy when local installation is not possible for example due to permission reasons.

Container availability
======================

Singularity runs natively only on Linux; on Windows it requires either WSL2 or a virtual machine (VM); on macOS it is available as a beta version or it requires a VM. Docker runs natively on Linux and Windows, while on macOS it can be installed as an app that makes use of a VM under the hood.

.. note::
   WSL2 (Windows Subsystem for Linux 2) is a native Windows 10 feature that allows to run a GNU/Linux terminal without the use of a VM. It can be installed following the official `documentation <https://docs.microsoft.com/en-us/windows/wsl/install-win10#requirements>`__. The Windows filesystem is mounted under ``/mnt``, e.g. ``cd mnt/c/Users/yourusername/Documents``.

When working on Linux machines, Docker produces output files that require root permissions to be handled (e.g. deleted), which is an issue for users who don't have root permissions. Running Docker on Windows does not have such problems because the user has more control on output files permissions. On macOS [TODO: it also shouldn't have root permissions]. On the contrary, Singularity has the advantage of always producing output files with non-root privileges.

The table below summarizes relevant differences between Singularity and Docker containers:

======================  ==============  ======
Feature                 Singularity     Docker
======================  ==============  ======
Runs on Linux           Y               Y
Runs on Windows         Y (WSL2 or VM)  Y
Runs on macOS           Y (beta or VM)  Y (VM)
Root privilege needed   N               Y
Mounting input files    N               Y
======================  ==============  ======


Singularity (suggested)
-----------------------
The machine where *ksrates* will be executed (either a local computer or a remote compute cluster) needs to have Singularity installed. More information can be found on the Singularity 3.7 installation `page <https://sylabs.io/guides/3.7/admin-guide/installation.html#installing-singularity>`__ or have an overview of the documentation `history <https://sylabs.io/docs/>`__ for up-to-date or version-specific instructions.

The Singularity container will be downloaded from ``vibpsb/ksrates`` repository on Docker Hub when launching the Nextflow pipeline. Successive runs will use the local copy.

When using the *ksrates* Nextflow pipeline, the only other dependency that must be installed is Nextflow (for more information see the official installation `page <https://www.nextflow.io/docs/latest/getstarted.html#requirements>`__).

*   First install Java 8 or later. ``default-djk`` works as well::

        sudo apt-get install default-jdk

*   Then install Nextflow in one of the following ways:

        *   Through ``wget``::
        
                wget -qO- https://get.nextflow.io | bash

        *   Through ``bioconda`` (for more info on how to setup ``bioconda`` see this `page <https://bioconda.github.io/user/install.html>`__)::

                conda install nextflow

*   Optionally make Nextflow accessible by your ``$PATH`` variable, for example::

        mv nextflow /usr/local/bin 


Docker
------

The machine where *ksrates* will be executed (either a local computer or a remote compute cluster) needs to have Docker installed. More information can be found on the Docker installation `page <https://docs.docker.com/get-docker/>`__.

The Docker container will be downloaded from ``vibpsb/ksrates`` repository on Docker Hub when launching the Nextflow pipeline. Successive runs will use the local copy.

When using the *ksrates* Nextflow pipeline, the only other dependency that must be installed is Nextflow (for more information see the official installation `page <https://www.nextflow.io/docs/latest/getstarted.html#requirements>`__).

*   First install Java 8 or later. ``default-djk`` works as well::

        sudo apt-get install default-jdk

*   Then install Nextflow in one of the following ways:

        *   Through ``wget``::
        
                wget -qO- https://get.nextflow.io | bash

        *   Through ``bioconda`` (for more info on how to setup ``bioconda`` see this `page <https://bioconda.github.io/user/install.html>`__)::

                conda install nextflow

*   Optionally make Nextflow accessible by your ``$PATH`` variable, for example::

        mv nextflow /usr/local/bin 


Local installation
==================

Without the use of a container the installation of *ksrates* and its dependencies has to be carried out manually.

1.  Clone the *ksrates* repository from `GitHub <https://github.com/VIB-PSB/ksrates>`__::

    	git clone https://github.com/VIB-PSB/ksrates

2.  Install the tool and its Python dependencies using ``pip3``::

    	cd ksrates
    	pip3 install .

3.  Non-Python dependencies can be installed in two possible ways.

		*   Through ``apt-get`` and ``wget``::

				sudo apt-get -yq install default-jdk build-essential ncbi-blast+ muscle mafft prank fasttree mcl phyml paml
				wget -qO- https://get.nextflow.io | bash

		*   Through ``apt-get`` and ``bioconda`` (more info on how to setup ``bioconda`` `here <https://bioconda.github.io/user/install.html>`__)::
				
				TODO: check these two commands
				sudo apt-get -yq install default-jdk build-essential
				conda install muscle blast mafft prank fasttree mcl phyml paml nextflow

   Optionally make Nextflow accessible by your ``$PATH`` variable, for example::

        mv nextflow /usr/local/bin 

4. Install I-ADHoRe 3.0 from its GitHub `page <https://github.com/VIB-PSB/i-ADHoRe>`__ (required only for collinearity analysis of genome data).


Testing your installation
=========================

1.  Clone the GitHub repository to get the use case dataset.

        git clone https://github.com/VIB-PSB/ksrates

2.  Access the ``test`` directory in a terminal::

        cd ksrates/test

2.  Launch *ksrates*:

    * Through a Singularity container::
     
        nextflow run VIB-PSB/ksrates --config ./config_elaeis.txt -with-singularity docker://vibpsb/ksrates:latest

    * Through a Docker container::
        
        nextflow run VIB-PSB/ksrates --config ./config_elaeis.txt -with-docker vibpsb/ksrates
    
    * With local installation::
    
        nextflow run VIB-PSB/ksrates --config ./config_elaeis.txt