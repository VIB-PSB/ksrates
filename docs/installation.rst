Installation
************

*ksrates* is available as a `Singularity <https://singularity.hpcng.org>`__ or `Docker <https://www.docker.com>`__ container, which bundle *ksrates* and all required external software dependencies, and as a Python package, which requires manual installation of the package and its dependencies but allows for more flexibility in integrating it into existing bioinformatics environments and toolchains. In addition to a simple and easy to use command-line interface (CLI), we also provide a user-friendly `Nextflow <https://www.nextflow.io>`__ pipeline that allows to run a complete *ksrates* analysis fully automated.

*ksrates* runs on any Linux or macOS system or on Windows with Windows Subsystem for Linux 2 (WSL2) or a virtual machine installed. However, *ksrates* analyses are computationally demanding, and therefore we recommend the use of a computer cluster (or cloud platform) for any but the simplest data sets.

.. _`install_nextflow`:

Nextflow (recommended)
======================

The *ksrates* Nextflow pipeline makes it easy to execute the individual steps of a *ksrates* analysis as a single command. Nextflow also makes it easy to configure the execution of the pipeline on a variety of computer clusters (see the :ref:`nextflow_config_section` section).

To install `Nextflow <https://www.nextflow.io>`__ and its dependencies, follow the commands below (or the official `Nextflow installation instructions <https://www.nextflow.io/docs/latest/getstarted.html>`__).

*   Make sure you have Bash 3.2 (or later) installed.

*   If you do not have `Java <https://www.oracle.com/java/>`__ installed, install `Java 11 or later <https://www.oracle.com/java/technologies/javase-downloads.html>`__; on Linux you can for example use::

        sudo apt-get install default-jdk  
   
   .. note::
      Nextflow versions before ``22.01.x-edge`` `require <https://www.nextflow.io/blog/2022/evolution-of-nextflow-runtime.html>`__ Java version 8 up to 15 for their execution.

*   Then install Nextflow using either::

        wget -qO- https://get.nextflow.io | bash

    or::

        curl -fsSL https://get.nextflow.io | bash


    This creates the ``nextflow`` executable file in the current directory.

   .. note::
      As from `ksrates` ``v2.0.0``, the Nextflow pipeline has been ported to DSL2 syntax and requires at least Nextflow version ``22.03.0-edge``. Concerning older `ksrates` versions written with DSL1, from ``v1.1.3`` to ``v1.1.5`` they support from ``22.03.0-edge`` until before ``22.12.0-edge``, while ``v1.1.2`` and previous ones require earlier versions than ``22.03.0-edge``.

*   Optionally make the ``nextflow`` executable accessible by your ``$PATH`` variable, for example by moving it::

        sudo mv nextflow /usr/local/bin 


The *ksrates* Nextflow pipeline itself does not need to be installed and will be automatically downloaded and set up simply when you execute the launch of the *ksrates* Nextflow pipeline for the first time. The same applies to the *ksrates* package and its dependencies when using the *ksrates* Singularity or Docker container (see below). In other words, if you plan on only using the *ksrates* Nextflow pipeline with a container it is not necessary to manually download or install *ksrates* itself, the only other software you may need to install is the Singularity or Docker software itself (next section). (Note that the *ksrates* Nextflow pipeline code is however also included in the *ksrates* GitHub repository and can thus also be executed from a manually installed *ksrates* package (see the Manual installation section below).)


Singularity and Docker containers
=================================

Containers are standalone portable runtime environments that package everything needed to run a software, including application code, external software dependencies and operating system libraries and runtime, and are thus executable in any computing environment for which Singularity and Docker container engines are available. This comes in handy, for example, when local installation of *ksrates* and its software dependencies are not possible, for instance due to permission issues, or for deploying *ksrates* to a computer cluster or cloud.

Availability and dependencies
-----------------------------

`Singularity <https://singularity.hpcng.org>`__ runs natively only on Linux. On Windows it requires either WSL2 (recommended; see Note below) or a virtual machine (VM). On macOS it is available as a beta version or it also requires a VM.
Singularity has the advantage over Docker of always producing output files with non-root permissions.

.. note::
   WSL2 (Windows Subsystem for Linux 2) is a native Windows 10 feature that allows to run a GNU/Linux terminal without the use of a VM. It can be installed following the official `documentation <https://docs.microsoft.com/en-us/windows/wsl/install-win10#requirements>`__.

`Docker <https://www.docker.com>`__ runs natively on both Linux and Windows, while on macOS it can be installed as an application that makes use of a VM under the hood.
When working on Linux machines, Docker produces output files that require root permissions to be handled (e.g. to delete them), which is an issue for users who don't have root permissions. Running Docker on Windows and macOS does not have such problems because the user has more control on output file permissions.

The table below summarizes relevant differences between Singularity and Docker containers/engines:

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
When using the *ksrates* Singularity container, either to run the *ksrates* CLI or Nextflow pipeline, the machine (i.e. a local computer or a remote computer cluster or cloud node) needs to have `Singularity <https://singularity.hpcng.org>`__ installed (*ksrates* has been tested with version 3.7). More information can be found on the `Singularity Quick Start page <https://singularity.hpcng.org/user-docs/master/quick_start.html>`__.
For a Linux installation we suggest to follow the `Install from Source section <https://singularity.hpcng.org/admin-docs/master/installation.html#before-you-begin>`__ (*Install Dependencies*, *Install Go*, *Download Singularity from a release* and *Compile Singularity*).
For up-to-date and version-specific instructions, please refer to this `page <https://singularity.hpcng.org/docs/>`__.

.. note::
   To allow users to run the pipeline from any directory in a cluster (i.e. not necessarily from their home directory), the `user bind control <https://singularity.hpcng.org/admin-docs/master/configfiles.html?highlight=user%20bind%20control#bind-mount-management>`__ feature needs to be left active during Singularity installation [Default: "YES"].

When using the *ksrates* Nextflow pipeline with the *ksrates* Singularity container, the container will be automatically downloaded from the ``vibpsb/ksrates`` repository on Docker Hub on first launch (this may take awhile depending on your Internet connection speed since the container has a size of about 1 GB) and will then be stored and reused for successive runs.

Docker
------

When using the *ksrates* Docker container, either to run the *ksrates* CLI or Nextflow pipeline, the machine (i.e. a local computer or a remote computer cluster or cloud node) needs to have Docker installed. More information can be found on the Docker installation `page <https://docs.docker.com/get-docker/>`__.

When using the *ksrates* Nextflow pipeline with the *ksrates* Docker container, the container will be automatically downloaded from the ``vibpsb/ksrates`` repository on Docker Hub on first launch (this may take awhile depending on your Internet connection speed since the container has a size of about 1 GB) and will then be stored and reused for successive runs.

.. _`manual_installation`:

Manual installation
===================

When not using or not being able to use one of the *ksrates* containers, for example to integrate the tool into existing bioinformatics environments and toolchains, the installation of the *ksrates* Python package and its dependencies can or has to be carried out manually. The following commands guide you through the installation on a Linux machine. Windows users can carry out the installation with the same commands by using either WSL2 (recommended; see Note below) or a virtual machine (VM) with Linux installed. macOS users can for example use `Homebrew <https://brew.sh>`__ instead of ``apt-get``.

.. note::
   WSL2 (Windows Subsystem for Linux 2) is a native Windows 10 feature that allows to run a GNU/Linux terminal without the use of a VM. It can be installed following the official `documentation <https://docs.microsoft.com/en-us/windows/wsl/install-win10#requirements>`__.


1.  Most of the non-Python dependencies can be installed with the following commands::

            sudo apt-get update && sudo apt-get -yq install python3-pip default-jdk build-essential ncbi-blast+ muscle fasttree mcl phyml | bash

    .. note::
       The workflow was developed using BLAST 2.6.0+ and it might error when using a much older versions such as 2.12.0+ (see GitHub issue `#48 <https://github.com/VIB-PSB/ksrates/issues/48>`__ ). Moreover it is only compatible with MUSCLE v3 (`v3.8.31 <https://drive5.com/muscle/downloads_v3.htm>`__) and not with v5 (see GitHub issue `#41 <https://github.com/VIB-PSB/ksrates/issues/41>`__).


2.  Install PAML 4.9j from source (for more information see PAML installation `page <http://abacus.gene.ucl.ac.uk/software/#phylogenetic-analysis-by-maximum-likelihood-paml>`__) to avoid compatibility issues::

        wget http://abacus.gene.ucl.ac.uk/software/paml4.9j.tgz
        tar -xzf paml4.9j.tgz
        cd paml4.9j/src && make -f Makefile

    Then make the executable ``codeml`` available through the ``$PATH`` variable (the downloaded PAML directory can be deleted):
    
        *   Either move ``codeml`` to a directory already present in ``$PATH``, e.g. ``usr/local/bin``::

                sudo mv codeml usr/local/bin
        
        *   Or move ``codeml`` to another directory (here assumed to be ``~/bin``) and add this directory to ``$PATH``, for the Bash shell by copying the following line to the shell initialization file (e.g. ``.bashrc``)::

                export PATH=$PATH:~/bin

3. Install i-ADHoRe 3.0 from its GitHub `page <https://github.com/VIB-PSB/i-ADHoRe>`__ (required only for collinearity analysis of genome data for the focal species).

4.  Clone the *ksrates* repository from `GitHub <https://github.com/VIB-PSB/ksrates>`__ and install the package and its Python dependencies::

        git clone https://github.com/VIB-PSB/ksrates
        # Starting from v2.0.0, download also this compressed file:
        wget https://zenodo.org/records/15225340/files/original_angiosperm_sequences.tar.gz -P ksrates/ksrates/reciprocal_retention
    	# Move to the ksrates subdirectory and install the package with ``pip``
        cd ksrates
    	pip3 install .


Testing your installation
=========================

.. note::
    WSL2 users can enter the Windows file system from the terminal through e.g. ``cd mnt/c/Users/your_username``.

1.  Clone the *ksrates* repository from `GitHub <https://github.com/VIB-PSB/ksrates>`__ to get the use case dataset::

        git clone https://github.com/VIB-PSB/ksrates
        # Starting from v2.0.0, download also this compressed file:
        wget https://zenodo.org/records/15225340/files/original_angiosperm_sequences.tar.gz -P ksrates/ksrates/reciprocal_retention


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
    # Starting from v2.0.0, download also this compressed file:
    wget https://zenodo.org/records/15225340/files/original_angiosperm_sequences.tar.gz -P ksrates/ksrates/reciprocal_retention
    # Move to the ksrates subdirectory and install the package with ``pip``
    cd ksrates
    pip3 install .