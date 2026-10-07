.. _`paralog_ks_database`:

Paralog *K*:sub:`S` database
****************************

.. _`paralog_ks_database_overview`:

Reusing previously computed *K*:sub:`S` data
============================================

Paralog *K*:sub:`S` data can be stored in a self-hosted `sqld <https://github.com/tursodatabase/libsql>`__
database, instead of existing only in local output TSV files (see :ref:`ks_estimate_output_wgd`).
This is particularly useful to centralize *K*:sub:`S` data across multiple independent *ksrates* analyses.

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


.. _`paralog_ks_database_setup`:

Setting up the database server
==============================

Configure the database server
-----------------------------

Enable the use of the database via parameter ``use_paralog_ks_database`` in each expert configuration 
file of the *ksrates* datasets of interest; this is off by default, and without it *ksrates* relies as usual 
solely on the local TSV output files.

Choose a directory that will host the database server, i.e. where the address file (``paralog_ks_server_address.txt``), ``paralog_ks_sqld_data/`` and ``paralog_ks_sqld_keys/`` will be generated. Note: if you want several datasets to share one database, prefer a central directory above/alongside them (e.g. ``/path/to/ks_analaysis_dir/paralog_ks_database``).

Point then every dataset's configuration file at the address file, by setting parameter ``ks_list_paralog_database_path`` to the absolute path of ``<location>/<address_filename>``::

    ks_list_paralog_database_path = /path/to/ks_analaysis_dir/paralog_ks_database/paralog_ks_server_address.txt

Start the server with one of the methods below, providing the chosen directory via ``--location``. This step must be executed once, independently and before of any *ksrates* analysis, and left running for as long as needed (e.g. a self-standing job on the cluster).

Starting it via the Nextflow pipeline (recommended)
---------------------------------------------------

A separate Nextflow pipeline (``setup_database_server.nf``) is available::

    nextflow run VIB-PSB/ksrates -main-script setup_database_server.nf \
        -c nextflow.config \
        -profile apptainer \
        -process.executor=local \
        --location /path/to/ks_analaysis_dir/paralog_ks_database

Provide via ``-c nextflow.config`` the Nextflow configuration of a *ksrates* analysis, to make it use the container;
configure ``-profile`` with the container engine in use (e.g. ``apptainer``);
with ``-process.executor=local``, the pipeline steps are not further submitted to the cluster and remain instead local to the ``nextflow run`` invocation itself;
``--location`` is the central directory that will hold the server's data/keys directories and
address file;
``--port`` (port ``sqld`` listens on, default 8080) and ``--address-filename`` (``paralog_ks_server_address.txt`` by default) can also be overridden if needed.


Starting it via the CLI command
-------------------------------

Alternatively, you can start the database server by directly invoking the *ksrates* command ``launch-paralog-ks-server``,
optionally leveraging a container as shown in :ref:`manual_pipeline`::

    ksrates launch-paralog-ks-server \
        --location /path/to/ks_analaysis_dir/paralog_ks_database

The same options explined above apply here (``--location``, ``--port`` and ``--address-filename``).
When running outside of a container, this command downloads the ``sqld`` binary on first launch to ``$XDG_CACHE_HOME/ksrates/sqld_bin``
(by default ``~/.cache/ksrates/sqld_bin``), and reuses it from there afterward.


Restarting the server
=====================

The server process and the data it serves are two separate things. Cancelling or losing the
server job (e.g. it hits a walltime limit) does not delete the database data. What goes stale is the address file, since it recorded a host and port that no longer has a server listening on it.

To pick the database back up, relaunch the server the same way (see :ref:`paralog_ks_database_setup`), with the same ``--location``. It reuses the existing ``paralog_ks_sqld_data/``
(so no data is lost) and the existing ``paralog_ks_sqld_keys/`` (so previously-issued tokens keep
working), and simply writes a fresh ``host:port`` to the address file. Every dataset's configuration file keeps pointing at the same address file path, so no reconfiguration is needed.


Inspecting and managing the database content
==============================================

The command ``inspect-paralog-ks-db`` dumps the full database content to a TSV
(gene pairs) plus per-species i-ADHoRe files, for inspection outside of *ksrates* (Excel, pandas...).
Two options change this behavior:

* ``--list``: for a quick overview of what's in the
  database, instead exports a TSV with one row per species and a True/False column per analysis
  type (paranome, anchors, reciprocally retained, i-ADHoRe files).
* ``--delete``: deletes every species matching the ``SPECIES_FILTER`` argument, after
  listing the matches and asking for confirmation; useful to discard incomplete or outdated data.

  .. note ::
        ``SPECIES_FILTER`` matches latin names as a case-insensitive substring; quote it if it includes
        spaces (e.g. ``"Elaeis guineensis"``)::

            ksrates inspect-paralog-ks-db paralog_ks_server_address.txt --list
            ksrates inspect-paralog-ks-db paralog_ks_server_address.txt "Elaeis guineensis" --delete


Backfilling from existing TSV output
======================================

Species processed without making use of the database already have their *K*:sub:`S` data stored in local TSV files, but nothing in the shared
database yet. Command ``populate-paralog-ks-db`` backfills the database from that existing output,
without re-running any pipeline::

    ksrates populate-paralog-ks-db configs/ paralog_distributions/ --database paralog_ks_server_address.txt

Here, ``configs/`` is a directory of *ksrates* configuration files (used to resolve each species' latin
name); ``paralog_distributions/`` is the directory containing the ``wgd_*`` subdirectories with the
TSV output to extract from.

Like the skip logic above, only analysis types a species doesn't already have in the database
are added by default (add ``--force`` to re-extract and overwrite everything instead),
so this command is safe to (re-)run any time.
