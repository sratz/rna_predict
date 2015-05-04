# DCA Filtering Syntax

Sometimes it may be necessary to filter and modify DCA contacts.

To specify filters on the command line the following syntax is used:

    filter:<arg>:<arg>:<arg>:...

The number of arguments depends on the filter being used.

Multiple filters may be chained by separating them with `,`.

The following filters are implemented:

## `none`

The `none` filter works as a dummy filter and does nothing. It takes no
options.

## `threshold`

The `threshold` filter looks up DCA contacts in a target PDB file,
calculates their distances, and discards or includes them depending on
the threshold setting.

The syntax is the following:

    threshold:<n>:<cst>:<mode>:<moodel>

Arguments:

- `n`: threshold (< 0: keep below, > 0: keep above)
- `cst`: prediction / constraints selection to look up PDB file
- `mode`, `model`: model selection mode (see [README.md](../README.md) for details)

Example:

    threshold:8.0:100dca_FADE_-100_26_20_-2_2:cluster:1,threshold:-6.0:100dca_FADE_-100_26_20_-2_2:cluster:1

This uses the first cluster of the `100dca_FADE_-100_26_20_-2_2`
predicton as target PDB file. All DCA contacts whose distance is between
6 and 8 Ångström in the target PDB are kept.
