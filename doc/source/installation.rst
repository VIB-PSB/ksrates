Installation
************

*ksrates* can be either installed locally on the machine that will run the pipeline (local computer or remote compute cluster) or can else be shipped in a Singularity and Docker container. Such containers are isolated portable environments in which the package and most dependencies are already installed, which comes in handy when local installation is not possible due for example to permission reasons.

Container availability
======================

Singularity runs natively only on Linux; on Windows it requires either WSL2 (suggested; see Note below) or a virtual machine (VM); on macOS it is available as a beta version or it requires a VM.
Singularity has the advantage over Docker of always producing output files with non-root permissions.

.. note::
   WSL2 (Windows Subsystem for Linux 2) is a native Windows 10 feature that allows to run a GNU/Linux terminal without the use of a VM. It can be installed following the official `documentation <https://docs.microsoft.com/en-us/windows/wsl/install-win10#requirements>`__.

Docker runs natively on Linux and Windows, while on macOS it can be installed as an app that makes use of a VM under the hood.
When working on Linux machines, Docker produces output files that require root permissions to be handled (e.g. deleted), which is an issue for users who don't have root permissions. Running Docker on Windows and macOS does not instead have such problems because the user has more control on output file permissions. 

The table below summarizes relevant differences between Singularity and Docker containers:

.. include:: <isopub.txt>
.. table:: Supported (|check|) and unsupported (|cross|) features.

    ======================  ====================  ============
    Feature                 Singularity           Docker
    ======================  ====================  ============
    Runs on Linux           |check|               |check|
    Runs on Windows         |check| (WSL2 or VM)  |check|
    Runs on macOS           |check| (beta or VM)  |check| (VM)
    Root privilege needed   |cross|               |check|
    ======================  ====================  ============


Singularity (suggested)
-----------------------
The machine where *ksrates* will be executed (either local computer or remote compute cluster) needs to have Singularity installed. More information can be found in the Singularity 3.7 installation `page <https://sylabs.io/guides/3.7/admin-guide/installation.html#installing-singularity>`__.
For Linux installation we suggest to follow the *Install from Source* `section <https://sylabs.io/guides/3.7/admin-guide/installation.html#before-you-begin>`__ (*Install Dependencies*, *Install Go*, *Download Singularity from a release* and *Compile Singularity*).
For up-to-date and version-specific instructions, please refer to this `page <https://sylabs.io/docs/>`__.

.. note::
   To allow users to run the pipeline from any directory in a cluster (i.e. not necessarily from their home directory), the `user bind control <https://sylabs.io/guides/3.7/admin-guide/configfiles.html?highlight=user%20bind%20control#bind-mount-management>`__ feature needs to be left active during Singularity installation [Default: "YES"].

When using the *ksrates* Nextflow pipeline, the only other dependency that must be installed is Nextflow (for more information see its official installation `page <https://www.nextflow.io/docs/latest/getstarted.html#requirements>`__).

*   First install Java 8 or later. ``default-djk`` works as well::

        sudo apt-get install default-jdk

*   Then install Nextflow::

        wget -qO- https://get.nextflow.io | bash

*   Optionally make Nextflow accessible by your ``$PATH`` variable, for example::

        sudo mv nextflow /usr/local/bin 

When launching the Nextflow pipeline with Singularity, the container will be downloaded from ``vibpsb/ksrates`` repository on Docker Hub and the local copy will be used for successive runs.

Docker
------

The machine where *ksrates* will be executed (either local computer or remote compute cluster) needs to have Docker installed. More information can be found on the Docker installation `page <https://docs.docker.com/get-docker/>`__.

When using the *ksrates* Nextflow pipeline, the only other dependency that must be installed is Nextflow (for more information see its official installation `page <https://www.nextflow.io/docs/latest/getstarted.html#requirements>`__).

*   First install Java 8 or later. ``default-djk`` works as well::

        sudo apt-get install default-jdk

*   Then install Nextflow::

        wget -qO- https://get.nextflow.io | bash

*   Optionally make Nextflow accessible by your ``$PATH`` variable, for example::

        sudo mv nextflow /usr/local/bin 

When launching the Nextflow pipeline with Docker, the container will be downloaded from ``vibpsb/ksrates`` repository on Docker Hub and the local copy will be used for successive runs.


Local installation
==================

Without the use of a container the installation of *ksrates* and its dependencies has to be carried out manually. The following commands guide through the installation on a Linux machine; Windows users can carry out the installation with the same commands by using either WSL2 (suggested; see Note below) or a virtual machine (VM) with Linux installed.

.. note::
   WSL2 (Windows Subsystem for Linux 2) is a native Windows 10 feature that allows to run a GNU/Linux terminal without the use of a VM. It can be installed following the official `documentation <https://docs.microsoft.com/en-us/windows/wsl/install-win10#requirements>`__.


1.  Most of non-Python dependencies can be installed with the following commands::

        sudo apt-get update && sudo apt-get -yq install python3-pip default-jdk build-essential ncbi-blast+ muscle mafft prank fasttree mcl phyml
	wget -qO- https://get.nextflow.io | bash

    Optionally make Nextflow accessible through your ``$PATH`` variable, for example::

        sudo mv nextflow /usr/local/bin
    
2.  Install PAML 4.9j from source (for more infromation see PAML installation `page <http://abacus.gene.ucl.ac.uk/software/#phylogenetic-analysis-by-maximum-likelihood-paml>`__) to avoid compatibility issues::

        wget http://abacus.gene.ucl.ac.uk/software/paml4.9j.tgz
        tar -xzf paml4.9j.tgz
        cd paml4.9j/src && make -f Makefile

    Then make the executable ``codeml`` available through the ``$PATH`` variable (the downloaded PAML directory can be deleted):
    
        *   Either move ``codeml`` to a directory already present in ``$PATH``, e.g. ``usr/local/bin``::

                sudo mv codeml usr/local/bin
        
        *   Or move ``codeml`` to another directory (here assumed to be ``~/bin``) and add this directory to ``$PATH``, for the Bash shell by copying the following line to the shell initialization file (e.g. ``.bashrc``)::

                export PATH=$PATH:~/bin

3. Install I-ADHoRe 3.0 from its GitHub `page <https://github.com/VIB-PSB/i-ADHoRe>`__ (required only for collinearity analysis of genome data).

4.  Clone the *ksrates* repository from `GitHub <https://github.com/VIB-PSB/ksrates>`__ and install the package and its Python dependencies::

        git clone https://github.com/VIB-PSB/ksrates
    	cd ksrates
    	pip3 install .


Testing your installation
=========================

.. note::
    WSL2 users can enter the Windows file system from the terminal through e.g. ``cd mnt/c/Users/your_username``.

1.  Clone the *ksrates* repository from `GitHub <https://github.com/VIB-PSB/ksrates>`__ to get the use case dataset::

        git clone https://github.com/VIB-PSB/ksrates

2.  Access the ``test`` directory in a terminal::

        cd ksrates/test

3.  Launch *ksrates* (the execution will take few minutes):

    * Through a Singularity container::
     
        nextflow run VIB-PSB/ksrates --config ./config_elaeis.txt -with-singularity vibpsb/ksrates:latest

    * Through a Docker container::
        
        nextflow run VIB-PSB/ksrates --config ./config_elaeis.txt -with-docker vibpsb/ksrates
    
    * With local installation::
    
        nextflow run VIB-PSB/ksrates --config ./config_elaeis.txt
