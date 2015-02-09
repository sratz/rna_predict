#!/usr/bin/env/python

# This script depends on the 'rnsd_b' script to be loaded.

# This script aligns a moving structure to a fixed one.
# The distances between the atom pairs are stored in the b-value fields of the moved structure.
# The b-values are then used to automatically colorize the structure using a spectrum.


from pymol import cmd


def align_and_color(moving, fixed, maximum=None):
    alignment_name = "alignment_" + moving
    cmd.align(moving, fixed, object=alignment_name, cycles=0)
    cmd.disable(alignment_name)
    rmsd_b(moving + " & " + alignment_name, fixed + " & " + alignment_name)
    cmd.color("white", fixed)
    cmd.spectrum("b", "green_red", moving, 0, maximum)
    cmd.show("cartoon")


cmd.extend("align_and_color", align_and_color)
