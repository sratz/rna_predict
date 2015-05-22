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

- Rosetta (Some minor modifications to the Rosetta source are required.
  For details check [doc/rosetta.md](doc/rosetta.md))
- Python packages:
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
    subprocess_buffsize = L

Note that this is done by prefixing all commands with `stdbuf -<option>` and
therefore requires the `stdbuf` program to be available.


## Usage

`rna_prediction` comes with a command line interface. To see all available
options run:

    rna_predict --help


## Peparation workflow

These steps need to be run only once for every structure.


### 1. Directory Preparation

The first step is to prepare a dedicated prediction directory. You need to have
at least a fasta file containing the sequence and another file containing the
secondary structure (in dot-bracket notation). Additionally, if available, you
can supply a PDB file containing the native crystal structure. To prepare the
directory, run the following:

    rna_predict prepare \
        --name <name> \
        --sequence <fasta_file> \
        --secstruct <secstruct_file> \
        --native <native_pdb_file>

All the options are optional and default values are used if not given.
See the output of `rna_predict prepare --help` for additional information.

This will parse out stems and motifs into the `preparation` subdirectory.

Note: When using the `--native` option, make sure the PDB file has neen adjusted
to work with rosetta. This means that the residues have to be reordered starting
at 1, so that Rosetta's `pdbslice.py` program uses the corrent ones. To achieve
this, prepare the file using the Rosetta tool:

    make_rna_rosetta_ready.py <native_pdb>


### 2. Create ideal helix models

For each helical part of the structure, an ideal helix model needs to be
generated:

    rna_predict create-helices


## Main workflow


### 1. Constraints creation

When tertiary constraints are to be used, `.cst` file should be created and put
int the `constraints` subdirectory. The syntax is explained in the Rosetta
documentation: <https://www.rosettacommons.org/docs/latest/constraint-file.html>

The file can be created by hand or using the following tools:


#### 1a. Generation from DCA files

To create a tertiary constraints from DCA predictions, use the following command:

    rna_predict make-constraints \
        --dca-file <dca.dat> \
        --mapping-mode <mapping_mode> \
        --dca-count <count> \
        --cst-function <function> \
        --cst-out-file <output_name> \
        --filter <filter>

The input dca file should contain at least two columns containing the residue
numbers of the contact. The file should be sorted by DCA score. Optionally, the
first line of the file may contain a comment specifying a PDB-mapping such as the
following:

    # pdb-mapping: 10-12,44,80-90

This defines how the residue numbers in the DCA file are to be mapped. Rosetta
always uses 1,2,3,... internally, so the mapping above would, for example, result
in the residue number 12 in the DCA file to be mapped to prediction residue 3.

The `--mapping-mode` parameter specifies the method to map residue-residue contacts
to atom-atom contacts. Options are:

- allAtomWesthof
- pOnly

For details about the mapping, see [doc/atom-mapping.md](doc/atom-mapping.md).

The `--dca-count` option limits the number of predictions in the DCA input file.

The `--cst-function` sets the Rosetta function to use. See
<https://www.rosettacommons.org/docs/latest/constraint-file.html#Function-Types>
for details. The default function for constraints creation (`FADE -100 26 20 -2 2`)
uses a spline smoothed square well potential (represented by the "FADE" function)
and a default parameter set. After the generation of the cst file, it can of
course be fine-tuned by further modification in any text editor.

The `--cst-out-file` option specifies an output filename.

The `--filter` option allows the DCA contacts to be passed through a chain of
filters first. For the filter documentation see [doc/filters.md](doc/filters.md).


#### 1b. Editing existing file

To simply replace the Rosetta function in an existing `.cst` file you can use

    rna_predict edit-constraints \
        --cst <input_cst> \
        --cst-function <function> \
        --cst-out-file <output_name>

For option explanation see above.

This is pretty much the same as using search-and-replace in any text editor.


### 2. Prepare constraints

To run a simulation with a specific set of constraints (or none), another
preparation step needs to be run:

    rna_predict prepare-cst \
        --cst <cst> \
        --override-motifs-cst <motif_cst>

The `--cst` option selects the constraints from the `constraints` directory
to be prepared. If not given, a prediction called 'none' for no tertiary
constraints is created.

Optionally, it is possible to use a different set of motifs for the assembly.
For example you can create a common set of motif models and use this in all
future assemblies. To do this, specify the `--override-motifs-cst` option.


### 3. Motif creation

For all non-helical parts (loop regions, etc.) multiple models need to be
created. To do this, run the following:

    rna_predict create-motifs \
        --cst <cst> \
        --cycles <cycles> \
        --nstruct <nstruct> \
        --seed <random_seed> \
        --use-native

As always, `--cst` selects the constraints.

The `--cycles` option sets the number of monte-carlo cycles to run for
generating each model.

The `--nstruct` option sets the number of models created for each motif.

To override the initial random seed, you can specify `--seed`.

And to have Rosetta automatically calculate RMSD values to a native structure
you can supply the `--use-native` option.


### 4. Assembly

To combine helix and motif models an assembly simulation is run:

    rna_predict assemble \
        --cst <cst> \
        --cycles <cycles> \
        --nstruct <nstruct> \
        --seed <random_seed> \
        --use-native

The options are the same as the ones for `create-motifs`, but their default
values vary.

Note: The assembly step does not check how many models have already been
created so far.


### 5. Evaluation


#### 5a. Evaluation using Rosetta clustering and scoring

When the assembly has finished, you can evaluate the simulation. This means:

- Cluster the models
- Calculate RMSD values to the native structure, if available, and to
  the model with the best score.

Usage:

    rna_predict evaluate \
        --cst <cst> \
        --cluster-cutoff <cutoff> \
        --cluster-limit <limit> \
        --full-eval

The `--cluster-cutoff` option specifies the RMSD radius in angstrom after
which to create a new cluster.

The `--cluster-limit` option limits the maximum number of clusters to be
created.

The `--full-eval` option forces the whole evaluation to be run again, and
ignore any previous results stored.


#### 5b. Custom scoring

Due to the fact that DCA predictions are not perfect, a custom scoring method
was created. For each DCA prediction neighboring residues are included and if
the distance between any of these residue paris are in contact the score
is increased.

    rna_predict evaluate-custom \
        --cst <cst> \
        --dca-file <dca_file> \
        --dca-count <count> \
        --radius <radius> \
        --threshold <threshold> \
        --full-eval

For the `--dca-file` and `--dca-count` options see `make-constraints`.

The `--radius` option sets the number of neighboring residues to take into
account.

The `--threshold` option sets the distance threshold under which a residue
pair is treated as in-contact.


## Utilities


### Comparison with native structure

To print a summary of all predictions and the smallest RMSD-to-native
values of the first 1, 5, and 10 clusters, run

    rna_predict compare


### Status information

To print a summary of all predictions and their current state, run

    rna_predict status

The output table contains the columns P, M, A, and E to indicate the following:

- P: preparation step
- M: motif generation:
- A: assembly:
- E: evaluation

If a step is completed, "X" is shown, "-" otherwise.
For motif generation a "*" may be shown to indicate that models from a different
set of constraints are used.


### Model information and extraction

To print model information or extract PDB files use the following subcommands:

    rna_predict print-models|extract-models \
        --cst <cst> \
        --mode <selection_mode> \
        model [model ...]

The `--mode` option selects the way to look up models:

- cluster: Cluster number to reference the cluster primary model
- cluster_ntop: Clusters sorted by the RMSDs of their representatives
- ntop: Models sorted by RMSD to native structure
- tag: Internal model name
- top: Models sorted by Rosetta score

The `model` options may be string (if `mode` is 'tag'), or numbers. For `mode=cluster_ntop`
it may also be in the form of `n/m`, meaning the `n`th best cluster out of the first `m`
clusters.

Examples:

    rna_predict extract-models --mode=tag S_000289 S_000100  # extract two models by tag
    rna_predict extract-models --mode=top 1 2 3 4 5  # extract the two best-scoring models
    rna_predict extract-models --mode=ntop 1  # extract the model with the lowest native RMSD
    rna_predict extract-models --mode=cluster 1 2 3 4 5  # extract the cluster primaries of the first 5 clusters
    rna_predict extract-models --mode=cluster_ntop 1/5  # extract the lowest native RMSD cluster out of the first 5 clusters


### Evaluation tools

Plot generation and other tools can be accessed using

    rna_predict tools <tool> ...


#### plot-clusters

Plot score over native RMSD.

    rna_predict tools plot-clusters \
        --score-weights <score:weight,score:weight,...> \
        --max-models <max> \
        cst

The `--score-weights` options allows to calculate a different total model score
using the individual Rosetta scores. The score name "default" can be given to
set a default value for all other, non-specified scores.

For example, to only visualize the score of additional constraints, use:

    --score-weights default:0,atom_pair_constraint:1

For a list of score names, refer to the Rosetta documentation or use the
`print-models` command.

The `--max-models` option limits the number of models either by specifying a
number (> 1) or a fraction (<= 1.0).


#### plot-constraint-quality

This visualizes the distances of constraints by comparing it to a reference
(native) PDB structure.

    rna_predict tools plot-constraint-quality \
        --dca-mode \
        reference-pdb \
        cst|dca|filter [cst|dca|filter ...]

When `--dca-mode` is given, residue-residue distances are plotted, atom-atom
contacts otherwise.

For the `filter`-syntax see [doc/filters.md](doc/filters.md).


#### plot-contact-atoms

Plots atoms involved in forming nucleotide contacts that satisfy the cutoff
condition in the contact database.

    rna_predict tools plot-contact-atoms \
        --mean-cutoff <cutoff> \
        --std-cutoff <cutoff>

The `--mean-cutoff` and `--std-cutoff` options select the limits for the
average contact distances standard deviations.


#### plot-contact-distances

Plots histogram for each nucleotide pair contact containing the distances of
the atoms involved.

    rna_predict tools plot-contact-distances


#### plot-dca-contacts-in-pdb

Visualizes how well DCA contacts are fullfiled in PDB files.

    rna_predict tools plot-dca-contacts-in-pdb \
        dca-file \
        pdb-file [pdb-file -..]


#### plot-pdb-comparison

Compare PDB files by plotting the distances of the residues.

    rna_predict tools plot-pdb-comparison \
        ref-pdb \
        sample-pdb [sample-pdb ...]


#### plot-gdt

Create a GDT (gloabl distance test) plot.

This plots a distance cutoff on the y-axis and the percent of residues which
are below the cutoff on the x-axis.

    rna_predict tools plot-gdt \
        ref-pdb \
        sample-pdb [sample-pdb ...]
