``rna_predict``
==================

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


Installation and Configuaration
-------------------------------

The ``rna_predict`` Python package has the following dependencies:

* ``biopython``
* ``matplotlib``
* ``numpy``


To install the package (and the dependencies, if necessary), run::

    python setup.py install


Additionally, ``rna_predict`` depends on *Rosetta*. However, a few minor
modifications to the *Rosetta* source are required. For a detailed Rosetta
compilation and installation guide see ``doc/rosetta.rst`` or the HTML
documentation.

A few configuration options can be overridden in ``~/.rna_predict/sysconfig``.
For details refer to the documentation.


Usage
-----

``rna_predict`` comes with a command line interface. To see all available
options run::

    rna_predict --help

For detailed usage information refer to the documentation.


Documentation
-------------

The HTML documentation can be found `here <https://sratz.github.io/rna_predict/docs/>`_.

A PDF version can be downloaded from `here <https://sratz.github.io/rna_predict/docs/rna_predict.pdf>`_.

If you downloaded a source package or checked out the source from Git,
you can also build the HTML documentation yourself. For that, ``rna_predict``
needs to be installed as well as the ``Sphinx`` Python package.

To build the docs run::

    python setup.py build_sphinx

After that, the HTML output can be found in ``build/sphinx/html``.
