.. _`ks_database`:

*K*:sub:`S` database
********************

.. _`ks_database_overview`:

Reusing previously computed *K*:sub:`S` data
============================================

Paralog and ortholog *K*:sub:`S` data can be stored in a self-hosted `sqld <https://github.com/tursodatabase/libsql>`__
database, instead of existing only in local output TSV files (see :ref:`ks_estimate_output_wgd` and the ortholog TSV databases in :ref:`other_output`).
This is particularly useful to centralize *K*:sub:`S` data across multiple independent *ksrates* analyses,
which can then also read and write it simultaneously.

An analysis layout may consist of several independent directories, each running *ksrates* on
its own dataset. Without a shared database, each dataset computes and stores its *K*:sub:`S` data independently -
even for species that recur across datasets, whose expensive paranome, anchor-pair, reciprocally
retained and i-ADHoRe computations then get needlessly repeated every time::

    /path/to/ks_analaysis_dir/
    ├── dataset1/    (ksrates configuration and output for this dataset)
    ├── dataset2/
    ├── dataset3/
    └── ...

Pointing every dataset at the same running server removes this redundancy: once *K*:sub:`S` data
for a species has been computed by any one dataset, it's available to every other dataset that
later needs it.
If a later request asks for an additional analysis type for a species already in the database,
the workflow attempts to compute only that new type and reuse the rest. 

.. note::
    ``colinearity`` is the one common exception: it requires paranome's output as local files
    (e.g. the MCL gene-family output in ``paralog_distributions``), not from the database, so a
    local copy is required regardless of what's already stored. If a directory lacks those files -
    e.g. because they were cleaned up locally, or because paranome was originally computed by
    another dataset - the paranome step reruns to regenerate them, even though paranome
    *K*:sub:`S` data is already in the database. The recomputed data then replaces the old
    database entry for consistency.

    To prevent this, request ``paranome`` and ``colinearity`` together the *first* time a species is
    processed, where possible.


.. _`ks_database_setup`:

Setting up the database server
==============================

Configure the database server
-----------------------------

Enable the use of the database via parameter ``use_ks_database_server`` in each expert configuration 
file of the *ksrates* datasets of interest; this is off by default, and without it *ksrates* relies as usual 
solely on the local TSV output files. The database will cover both paralog and ortholog data.

Choose a directory that will host the database server, i.e. where the address file (``ks_database_server_address.txt``), ``ks_database_sqld_data/`` and ``ks_database_sqld_keys/`` will be generated. Note: if you want several datasets to share one database, prefer a central directory above/alongside them (e.g. ``/path/to/ks_analaysis_dir/ks_database``).

Point then every dataset's configuration file at the address file, by setting parameter ``ks_database_server_address_path`` to the absolute path of ``<location>/<address_filename>``::

    ks_database_server_address_path = /path/to/ks_analaysis_dir/ks_database/ks_database_server_address.txt

Start the server with one of the methods below, providing the chosen directory via ``--location``. This step must be executed once, independently and before of any *ksrates* analysis, and left running for as long as needed (e.g. a self-standing job on the cluster).

Starting it via the Nextflow pipeline (recommended)
---------------------------------------------------

Launch the database server through a Nextflow pipeline (``setup_database_server.nf``)::

    nextflow run VIB-PSB/ksrates -main-script setup_database_server.nf \
        -c nextflow.config \
        -profile apptainer \
        -process.executor=local \
        --location /path/to/ks_analaysis_dir/ks_database

Provide via ``-c nextflow.config`` the Nextflow configuration of a *ksrates* analysis, to make it use the container;
configure ``-profile`` with the container engine in use (e.g. ``apptainer``);
with ``-process.executor=local``, the pipeline steps are not further submitted to the cluster and remain instead local to the ``nextflow run`` invocation itself;
``--location`` is the chosen database directory;
``--port`` (port ``sqld`` listens on, default 8080) and ``--address-filename`` (by default ``ks_database_server_address.txt``) can also be overridden if needed.


Starting it via the CLI command
-------------------------------

Alternatively, you can start the database server by directly invoking the *ksrates* command ``launch-ks-database-server``,
optionally leveraging a container as shown in :ref:`manual_pipeline`::

    ksrates launch-ks-database-server \
        --location /path/to/ks_analaysis_dir/ks_database

The same options explained above apply here (``--location``, ``--port`` and ``--address-filename``).
When running outside of a container, this command downloads the ``sqld`` binary on first launch to ``$XDG_CACHE_HOME/ksrates/sqld_bin``
(by default ``~/.cache/ksrates/sqld_bin``), and reuses it from there afterward.


Restarting the server
=====================

The server process and the data it serves are two separate things. Cancelling or losing the
server job (e.g. it hits a walltime limit) does not delete the database data. What goes stale is the address file, since it recorded a host and port that no longer has a server listening on it.

To pick the database back up, relaunch the server the same way (see :ref:`ks_database_setup`), with the same ``--location``. It reuses the existing ``ks_database_sqld_data/``
(so no data is lost) and the existing ``ks_database_sqld_keys/`` (so previously-issued tokens keep
working), and simply writes a fresh ``host:port`` to the address file. Every dataset's configuration file keeps pointing at the same address file path, so no reconfiguration is needed.


Inspecting and managing the database content
==============================================

The command ``inspect-ks-db`` dumps the database content to files, for inspection outside of *ksrates*:
a TSV file with the paralog data plus per-species i-ADHoRe files, and two TSV files with the ortholog peaks and *K*:sub:`S` lists
in the exact format of the ortholog TSV databases (``peak_database_path`` and ``ks_list_database_path``).
The address file path must be passed as argument. An optional ``NAME_FILTER`` argument can be
added to restrict the action to species matching that string: use the latin name to filter paralogs data, or the two latin names joined by an underscore for ortholog data. The match is case-insensitive and must be quoted if it includes spaces. Option ``--type`` restricts the action to ``paralog`` or ``ortholog`` data (default: ``all``)::


    ksrates inspect-ks-db \
                ks_database_server_address.txt  "Elaeis guineensis" --type paralog

Two more options change the command behavior:

* ``--list``: for a quick overview of what's in the
  database, instead exports a TSV with one row per species (paralogs; with a True/False column per analysis
  type: paranome, anchors, reciprocally retained, i-ADHoRe files) and one with one row per species pair (orthologs; peak and *K*:sub:`S` list)::

            ksrates inspect-ks-db \
                ks_database_server_address.txt --list

* ``--delete``: instead deletes every species or species pair matching ``NAME_FILTER`` (required), 
  after listing the matches and asking for confirmation; requires ``--type paralog`` or ``--type ortholog``::

            ksrates inspect-ks-db \
                ks_database_server_address.txt "Elaeis guineensis" --type paralog --delete


Backfilling from existing TSV output
======================================

Species processed without making use of the database already have their *K*:sub:`S` data available in local TSV files, but nothing stored in the
database yet. Command ``populate-ks-db`` backfills the database from that existing output,
without re-running any pipeline.

For paralog data::

    ksrates populate-ks-db \
        --paralogs configs/ paralog_distributions/ --database ks_database_server_address.txt

Here, ``configs/`` is a directory containing at least one *ksrates* configuration file (used to resolve each species' latin
name); ``paralog_distributions/`` is the directory containing the *K*:sub:`S` TSV output to copy from; ``--database`` is the path to the address file.

For ortholog data::

    ksrates populate-ks-db \
        --orthologs ortholog_peak_db.tsv ortholog_ks_list_db.tsv --database ks_database_server_address.txt

The peaks and *K*:sub:`S` lists already stored in the ortholog TSV databases are copied as they are.

Both options can be given in the same call. Since only the data a species
(pair) doesn't already have in the database are added by default, use ``--force`` to re-extract and overwrite everything instead.
