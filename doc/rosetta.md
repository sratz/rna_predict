# Rosetta compilation and installation

`rna_prediction` relies on a few modifications of the Rosetta source code to
work.

## 1. Get Rosetta

Download Rosetta from the official website. The current set of patches targets
build `2015wk02` so it's advised to get a fairly recent weekly build.

## 2. Apply patches

The needed patches are included in the `contrib/rosetta_patches` directory.

To apply them `cd` to the unpacked Rosetta source and run:

    patch -p1 < /path/to/patches/*.patch

If they don't apply cleanly, for example when targeting a newer version of
Rosetta, they need to be applied manually. The changes are only minor so it
should be possible to do without much effort.

## 3. Build and install

Build Rosetta normally according the official documentation.

## 4. `PATH` setup

`rna_predict` needs to find the Rosetta binaries as well as some of the helper
scripts from the `tools` package.

Also, `rna_predict` does not specify the `-database` parameter when invoking
Rosetta commands, so the database needs to reside in the default location
or the `ROSETTA3_DB` environment variable needs to be set correctly.

The following lines in `.bashrc` or similar should be sufficient to make
everything work:

    export ROSETTA3_DB=/path/to/rosetta_database
    export RNA_TOOLS=/path/to/rosetta_tools/rna_tools
    export PYTHONPATH=$PYTHONPATH:$RNA_TOOLS/bin/
    export PATH=$PATH:/path/to/rosetta_source/bin:$ROSETTA_TOOLS/bin
