#!/usr/bin/env python

# This script changes the representation to make it suitable for a png export.


from pymol import cmd


def make_export_ready():
    cmd.show_as("cartoon")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("cartoon_tube_radius", "0.8")


cmd.extend("make_export_ready", make_export_ready)


