#!/usr/bin/env python

# helper script to quickly align structures and export them as a png file
# in the current directory

# requires align_and_color and make_export_ready to be loaded


from pymol import cmd


def mke(source, reference, maximum=None, export_filename=None):
    align_and_color(source, reference, maximum=maximum)
    make_export_ready()
    if export_filename:
        cmd.png(export_filename, height=1000, ray=1)


cmd.extend("mke", mke)
