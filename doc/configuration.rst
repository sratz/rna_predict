.. _configuration:

Configuration
-------------

``rna_predict`` reads the config file ``~/.rna_predict/sysconfig`` if present.


Rosetta path and suffix
^^^^^^^^^^^^^^^^^^^^^^^

By default, ``rna_predict`` searches for the Rosetta binaries in the ``PATH``
environment variable and assumes *.linuxgccrelease* as binary suffix. This
can be overridden in the config file as follows::

    [rosetta]
    exe_path = /path/to/rosetta/directory/bin/
    exe_suffix =


Other options
^^^^^^^^^^^^^

Currently, there is only one other option to control buffering of the output
of called programs. To enable line-buffered output, set the following::

    [rna_predict]
    subprocess_buffsize = L

Note that this is done by prefixing all commands with ``stdbuf -<option>`` and
therefore requires the ``stdbuf`` program to be available.
