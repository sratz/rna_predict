# `rna_prediction`

This python project models a complete workflow for RNA tertiary structure
prediction using the secondary structure information as well as additional
constraints.

The modelling is done using Rosetta and its `rna_denovo` protocol.

Additional constraints may come from a direct coupling analysis (DCA).
A utility to convert those results into all-atom constraints is included.

`rna_prediction` also comes with tools to cluster and compare the results with
a native crystal strcture .pdb as well as tools to visualize the results.


## Installation

### Dependencies

* Rosetta (Some minor modifications to the Rosetta source are required.
  For details check `docs/rosetta.txt`)
* Python packages:
  * `biopython`
  * `matplotlib`
  * `numpy`

### Installation of the python package

To install the `rna_prediction` package, run

    python setup.py install

This will also install the Python dependencies listed above, if necessary.

## Configuration

`rna_predict` reads the config file `~/.rna_predict/sysconfig` if present.

### Rosetta path and suffix

By default, `rna_predict` searches for the Rosetta binaries in the `PATH`
environment variable and assumes *.linuxgccrelease* as binary suffix. This
can be overridden in the config file as follows:

    [rosetta]
    exe_path = /path/to/rosetta/directory/bin/
    exe_suffix =

### Other options

Currently, there is only one other option to control buffering of the output
of called programs. To enable line-buffered output, set the following:

    [rna_predict]
    vsubprocess_buffsize = L

Note that this is done by prefixing all commands with `stdbuf -<option>` and
therefore requires the `stdbuf` program to be available.


## Usage

`rna_prediction` comes with a command line interface. To see all available
options run:

    rna_predict --help

# Workflow
RNA prediction is split into multiple steps.

TODO



## Important workflow information:

### Native pdb files:

When using the `--native` option while preparing, make sure the pdb file has
been adjusted to work with rosetta. This means that the residues have to be
reordered starting at 1, so that pdbslice.py takes the corrent ones. To
achieve this, prepare the file using the Rosetta tool:

    make_rna_rosetta_ready.py <native_pdb>

### Number of created structures in relation to the number of threads:
The `--nstruct` parameter behaves differently for the motif generation and the
assembly steps.

* Motif generation: The total number of structures is given and distributed
  among the threads.
* Assembly: Each thread is started with the total number of structures.

The reasoning behind this is that the number of structures for assembly is
not reached anyways. Instead, the current workflow uses a constant walltime.

