Installation
************

*ksrates* Python package and its dependencies can either be installed on the machine that will run the pipeline (a local computer or remote computer cluster) or can else be shipped in a Singularity or Docker container. Such containers are isolated portable environments in which the package and most dependencies are already installed, which comes in handy when local installation is not possible, for example due to permission issues. Moreover, beside being available through command line interface, *ksrates* can be launched through a Nextflow pipeline to ease the sequential concatenation of all required steps.

.. _`install_nextflow`:

Run as a Nextflow pipeline
==========================

It is possible to run the several *ksrates* commands as a single Nextflow pipeline. To install Nextflow and its dependencies, follow the commands below (for more information see its official installation `page <https://www.nextflow.io/docs/latest/getstarted.html#requirements>`__).

*   Make sure you have Bash 3.2 or later.

*   If you do not have `Java <https://www.oracle.com/java/>`__ installed, install `Java 8 (or later, up to 15) <https://www.oracle.com/java/technologies/javase-downloads.html>`__; on Linux you can for example use::

        sudo apt-get install default-jdk

*   Then install Nextflow::

        wget -qO- https://get.nextflow.io | bash

*   Optionally make Nextflow accessible by your ``$PATH`` variable, for example by moving the ``nextflow`` executable::

        sudo mv nextflow /usr/local/bin 



Container availability
======================

Singularity runs natively only on Linux; on Windows it requires either WSL2 (recommended; see Note below) or a virtual machine (VM); on macOS it is available as a beta version or it requires a VM.
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


Singularity (recommended)
-------------------------
The machine where *ksrates* will be executed (either a local computer or remote computer cluster) needs to have `Singularity <https://singularity.hpcng.org>`__ installed (*ksrates* has been tested with version 3.7). More information can be found on the `Singularity Quick Start page <https://singularity.hpcng.org/user-docs/master/quick_start.html>`__.
For Linux installation we suggest to follow the `Install from Source section <https://singularity.hpcng.org/admin-docs/master/installation.html#before-you-begin>`__ (*Install Dependencies*, *Install Go*, *Download Singularity from a release* and *Compile Singularity*).
For up-to-date and version-specific instructions, please refer to this `page <https://singularity.hpcng.org/docs/>`__.

.. note::
   To allow users to run the pipeline from any directory in a cluster (i.e. not necessarily from their home directory), the `user bind control <https://singularity.hpcng.org/admin-docs/master/configfiles.html?highlight=user%20bind%20control#bind-mount-management>`__ feature needs to be left active during Singularity installation [Default: "YES"].

When using the *ksrates* Nextflow pipeline, the only other dependency that must be installed is Nextflow (for more information see :ref:`install_nextflow` section).

When launching the Nextflow pipeline with Singularity, the container will be downloaded from the ``vibpsb/ksrates`` repository on Docker Hub and the local copy will be used for successive runs.

Docker
------

The machine where *ksrates* will be executed (either a local computer or remote computer cluster) needs to have Docker installed. More information can be found on the Docker installation `page <https://docs.docker.com/get-docker/>`__.

When using the *ksrates* Nextflow pipeline, the only other dependency that must be installed is Nextflow (for more information see :ref:`install_nextflow` section).

When launching the Nextflow pipeline with Docker, the container will be downloaded from the ``vibpsb/ksrates`` repository on Docker Hub and the local copy will be used for successive runs.


Manual installation
===================

Without the use of a container or to integrate the tool into existing environments/toolchains, the installation of *ksrates* Python package and its dependencies has to be carried out manually. The following commands guide through the installation on a Linux machine; Windows users can carry out the installation with the same commands by using either WSL2 (recommended; see Note below) or a virtual machine (VM) with Linux installed.

.. note::
   WSL2 (Windows Subsystem for Linux 2) is a native Windows 10 feature that allows to run a GNU/Linux terminal without the use of a VM. It can be installed following the official `documentation <https://docs.microsoft.com/en-us/windows/wsl/install-win10#requirements>`__.


1.  Most of non-Python dependencies can be installed with the following commands::

        sudo apt-get update && sudo apt-get -yq install python3-pip default-jdk build-essential ncbi-blast+ muscle mafft prank fasttree mcl phyml | bash

    When using the *ksrates* Nextflow pipeline, install Nextflow following the steps listed in :ref:`install_nextflow` section.
    
2.  Install PAML 4.9j from source (for more information see PAML installation `page <http://abacus.gene.ucl.ac.uk/software/#phylogenetic-analysis-by-maximum-likelihood-paml>`__) to avoid compatibility issues::

        wget http://abacus.gene.ucl.ac.uk/software/paml4.9j.tgz
        tar -xzf paml4.9j.tgz
        cd paml4.9j/src && make -f Makefile

    Then make the executable ``codeml`` available through the ``$PATH`` variable (the downloaded PAML directory can be deleted):
    
        *   Either move ``codeml`` to a directory already present in ``$PATH``, e.g. ``usr/local/bin``::

                sudo mv codeml usr/local/bin
        
        *   Or move ``codeml`` to another directory (here assumed to be ``~/bin``) and add this directory to ``$PATH``, for the Bash shell by copying the following line to the shell initialization file (e.g. ``.bashrc``)::

                export PATH=$PATH:~/bin

3. Install i-ADHoRe 3.0 from its GitHub `page <https://github.com/VIB-PSB/i-ADHoRe>`__ (required only for collinearity analysis of genome data).

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

3.  Launch *ksrates* (the execution will take few minutes)::
     
        nextflow run VIB-PSB/ksrates --config ./config_elaeis.txt

    Nextflow will download *ksrates* and will by default run the test pipeline on the ``local`` executor using the *ksrates* Singularity container, as configured in the included ``nextflow.config`` Nextflow configuration file (automatically detected). If needed, please adapt the configuration to the available resources (e.g. available CPUs/cores or switching to a Docker container or no container at all for a local installation) as described in the :ref:`nextflow_config_section` section.


Updating your installation
==========================

* To update the *ksrates* Nextflow pipeline to the latest release, run the following command::

        nextflow pull VIB-PSB/ksrates

* To update the Docker container image, run the following command to pull the new image from `Docker Hub <https://hub.docker.com/r/vibpsb/ksrates>`__::

        docker pull vibpsb/ksrates:latest

* To update the Singularity container image, first remove the old image (when using the *ksrates* Nextflow pipeline the image is stored in the ``cacheDir`` directory set in the ``nextflow.config`` or, if not set, by default in ``work/singularity`` in the project folder)::

        rm vibpsb-ksrates-latest.img

  The next time the pipeline is launched, Nextflow will automatically pull the new image from `Docker Hub <https://hub.docker.com/r/vibpsb/ksrates>`__.

 Alternatively, run the following command in the same directory of the old image to manually pull the new image from Docker Hub::
        
        singularity pull vibpsb-ksrates-latest.img docker://vibpsb/ksrates:latest

* To update your manual installation, uninstall the old version of *ksrates* package, clone the *ksrates* repository from `GitHub <https://github.com/VIB-PSB/ksrates>`__ and re-install the package::

    pip3 uninstall ksrates
    git clone https://github.com/VIB-PSB/ksrates
    cd ksrates
    pip3 install .