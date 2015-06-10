Welcome to rna_predict's documentation!
==========================================

The ``rna_predict`` Python package models a complete workflow for RNA tertiary
structure prediction using the secondary structure information as well as
optional, additional constraints.

The modelling is done using `Rosetta <https://www.rosettacommons.org/>`_
and its ``rna_denovo`` protocol.

Additional constraints, for example from co-evolutional anaylsis, can be
included to improve the prediction quality.

A utility to convert residue-residue contacts into atom-atom constraints
is included.

``rna_predict`` also comes with tools to cluster and visualize the prediction
results and to compare them to a native crystal strcture, if available.


Contents
========

.. toctree::
   :maxdepth: 4

   rosetta

   configuration
   usage
   atom-mapping
   filters

   rna_predict package and API documentation <rna_predict>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
