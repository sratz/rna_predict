#!/usr/bin/env python2
# coding: utf8

"""
Created on Sep 22, 2014

@author: sra
"""

import itertools
import pprint
import matplotlib.pyplot as plt
import numpy as np
import os.path
from os.path import basename, splitext

from rna_prediction import dcatools
from rna_prediction import pdbtools
from rna_prediction.simulation import RNAPrediction
from rna_prediction.simulation import EvalData
from rna_prediction.sysconfig import SysConfig


def plot_contact_distances():
        """Plots histogram for each nucleotide pair contact containing the distances of the atoms involved."""
        distance_map = dcatools.get_contact_distance_map()

        f, axarr = plt.subplots(4, 4)

        f.suptitle("normed frequency / distance")
        i = 0
        for res_pair, distance_map_res_pair in sorted(distance_map.iteritems()):
            print res_pair
            dists = []
            for atom_pair, distances_atom_pair in distance_map_res_pair.iteritems():
                dists += distances_atom_pair

            axarr[i / 4, i % 4].hist(dists, normed=True, bins=50)
            axarr[i / 4, i % 4].set_title(res_pair)
            axarr[i / 4, i % 4].set_xlim([0, 25])

            i += 1

        plt.show()


def plot_contact_atoms(mean_cutoff, std_cutoff):
    """
    Plots atoms involved in forming nucleotide contacts that satisfy the cutoff condition in the contact database.

    :param mean_cutoff: limit for average distance
    :param std_cutoff: limit for the standard deviation
    """
    distance_map = dcatools.get_contact_distance_map()
    mean_distance_map = dcatools.get_contact_distance_map_mean(distance_map, mean_cutoff, std_cutoff)

    f, axarr = plt.subplots(4, 4)

    f.suptitle("mean distance (cutoff: %s), standard deviation (cutoff: %s)" % (mean_cutoff, std_cutoff))
    f.subplots_adjust(left=0.02, bottom=0.05, right=0.99, top=0.95, hspace=0.3, wspace=0.1)
    i = 0

    for res_pair, mean_distance_map_res_pair in sorted(mean_distance_map.items()):
        print res_pair
        data = sorted(mean_distance_map_res_pair.items(), key=lambda entry: entry[1])
        x_names = [d[0] for d in data]
        x = np.arange(len(x_names))
        y = [d[1][0] for d in data]
        y_err = [d[1][1] for d in data]

        axarr[i / 4, i % 4].errorbar(x, y, y_err)
        axarr[i / 4, i % 4].set_xticks(x)
        axarr[i / 4, i % 4].set_xticklabels(x_names, rotation=90, size="x-small")
        axarr[i / 4, i % 4].text(0.05, .95, res_pair, transform=axarr[i / 4, i % 4].transAxes, va="top")
        i += 1

    plt.show()


def plot_constraint_quality(comparison_pdb, sources, dca_mode=False):
    """
    Plot constraint quality.

    This visualizes the distances of constraints by comparing it to a reference (native) PDB structure.

    :param comparison_pdb: filename to a pdb file to compare to
    :param sources: list of dca files or constraints (depending on dca_mode). may also use 'filter:' to filter on-the-fly
    :param dca_mode: visualize residue-residue DCA instead of atom-atom constraints
    """
    sim = RNAPrediction(SysConfig())
    pdb = pdbtools.parse_pdb("foo", comparison_pdb)
    chain = pdb[0].child_list[0]

    print "Comparison PDB: %s" % comparison_pdb
    plt.figure(figsize=(14, 8))

    if dca_mode:
        print "Mode: residue-residue"
        bins = np.linspace(0, 60, 60)
    else:
        print "Mode: atom-atom"
        bins = np.linspace(0, 60, 180)

    for arg in sources:
        print "Processing: %s" % arg
        if arg.startswith("filter:"):
            title = arg
            _, dca_file, filtertext = arg.split(":", 2)
            print "  Applying filters: %s" % filtertext
            dca = dcatools.parse_dca_data(dca_file)
            dca_filter_chain = sim.parse_dca_filter_string(filtertext)
            dcatools.filter_dca_data(dca_data=dca, dca_filter_chain=dca_filter_chain, quiet=True)

            if not dca_mode:
                print "  Mapping to atom-atom."
                cst_info = dcatools.build_cst_info_from_dca_contacts(dca, sequence=sim.config["sequence"], mapping_mode="allAtomWesthof", cst_function="FADE -100 26 20 -2 2", number_dca_predictions=100, quiet=True)
        else:
            title = splitext(basename(arg))[0]
            print "  Reading file directly."

            if dca_mode:
                dca = dcatools.parse_dca_data(arg)
            else:
                cst_info = sim.parse_cst_file(sim.parse_cst_name_and_filename(arg)[1])

        dists = []
        if dca_mode:
            j = 0
            for d in dca:
                if not d.use_contact:
                    continue
                j += 1
                if j > 100:
                    break
                average_heavy, minimum_heavy, minimum_pair = dcatools.get_contact_information_in_pdb_chain(d, chain)
                dists.append(minimum_heavy)
        else:
            for cst in cst_info:
                atom1, res1, atom2, res2, fuction = cst
                try:
                    dists.append(chain[res1][atom1] - chain[res2][atom2])
                except KeyError:
                    print atom1, res1, atom2, res2, "   NOT FOUND"

        print "  Average distance:", np.average(dists)

        plt.hist(dists, bins, alpha=0.6, label=title)
    plt.legend()
    plt.title("%s Constraints Quality %s" % ("Residue-Residue" if dca_mode else "Atom-Atom", os.path.basename(os.getcwd())))
    plt.xlabel(u"Native distance / Å")
    plt.ylabel("Number of contraints")
    plt.savefig("/tmp/rna_tools_quality_%s.png" % os.path.basename(os.getcwd()), bbox_inches="tight")
    plt.show()


def plot_dca_contacts_in_pdb(dca_prediction_filename, pdb_files):
    """
    Visualize how well DCA contacts are fullfiled in PDB files.

    This plots the distances of dca contacts in one or more PDB files.

    :param dca_prediction_filename: input DCA filename
    :param pdb_files: list of PDB filenames
    """
    dca = dcatools.parse_dca_data(dca_prediction_filename=dca_prediction_filename)

    # loop over all pdb files
    chains = []
    for p in pdb_files:
        pdb = pdbtools.parse_pdb("foo", p)
        chains.append({"chain": pdb[0].child_list[0], "name": p})

    maximum = min(len(dca), 100)
    averages = np.zeros((len(chains), maximum))
    minimums = np.zeros((len(chains), maximum))

    # loop over all dca entries
    for j, d in enumerate(dca):
        if j >= maximum:
            break

        print d

        for i in xrange(0, len(chains)):

            average_heavy, minimum_heavy, minimum_pair = dcatools.get_contact_information_in_pdb_chain(d, chains[i]["chain"])
            averages[i][j] = average_heavy
            minimums[i][j] = minimum_heavy
            print "        %-70s  avg: %.05f   min: %-3s  %-3s %f" % (chains[i]["name"], average_heavy, minimum_pair[0].name, minimum_pair[1].name, minimum_heavy)

    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    x = np.array(range(0, maximum))
    width = 0.8
    for i in xrange(0, len(chains)):
        color = colors.next()
        plt.bar(x + 1 - width / 2, averages[i], width, label=chains[i]["name"], color=color, alpha=0.5)
    plt.xticks(x + 1, rotation="vertical")
    plt.xlim((0, maximum + 1))
    plt.legend()
    plt.title(dca_prediction_filename)
    plt.xlabel("DCA Prediction Number")
    plt.ylabel(u"Distance / Å")
    plt.savefig("/tmp/rna_tools_checkcontacts_%s.png" % os.path.basename(os.getcwd()), bbox_inches="tight")
    plt.show()


def plot_clusters(cst, max_models=0.99):
    """
    Plot score over native rmsd.

    :param cst: constraints
    :param max_models: limit to number of models if > 1, or relative percentage if <= 1
    """
    sim = RNAPrediction(SysConfig())
    cst_name, cst_file = sim.parse_cst_name_and_filename(cst)
    eval_data = EvalData.load_from_cst(cst_name)

    number = int(max_models)
    if number <= 1:
        number = int(eval_data.get_model_count() * max_models)

    models = eval_data.get_models(range(1, number + 1), "top")

    models_c = []
    plots = []
    descs = []

    fig = plt.figure(figsize=(14, 8))

    # no cluster
    models_c.append([m for m in models if "cluster" not in m])
    descs.append("no cluster (%d)" % (len(models_c[0])))

    # clusters
    i = 1
    while True:
        models_current_cluster = [m for m in models if "cluster" in m and m["cluster"] == i]
        if len(models_current_cluster) == 0:
            n_clusters = i - 1
            break
        models_c.append(models_current_cluster)
        descs.append("cluster %d (%d)" % (i, len(models_current_cluster)))
        i += 1

    if "rmsd_native" in models_c[1][0]:
        comparison = "rmsd_native"
        plt.xlabel("native rmsd / A")
    else:
        comparison = "rmsd_cluster_1"
        plt.xlabel("rmsd to best structure / A")

    # create plots
    colors = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#ff6600', '#006000', '#600060', '#F7C2CA']
    shapes = ['o', 'D', 's', 'p']

    def pick_event(event):
        cluster = plots.index(event.artist)
        for ind in event.ind:
            print "cluster %d: %s" % (cluster, models_c[cluster][ind]["tag"])
            pprint.pprint(models_c[cluster][ind])

    for i in range(0, n_clusters + 1):
        color = "w" if i == 0 else colors[(i - 1) % len(colors)]
        shape = "o" if i == 0 else shapes[((i - 1) / len(colors)) % len(shapes)]
        plots.append(plt.scatter([x[comparison] for x in models_c[i]], [x["score"] for x in models_c[i]], s=35, c=color, marker=shape, picker=True))

    fig.canvas.mpl_connect('pick_event', pick_event)
    plt.title("model clusters %s, constraints: %s" % (sim.config["name"], cst_name))
    plt.ylabel("rosetta score")
    plt.legend(plots, descs, loc="upper right", prop={'size': 12})
    plt.savefig("/tmp/rna_tools_rmsdscore_%s.png" % os.path.basename(os.getcwd()), bbox_inches="tight")
    plt.show()


def plot_pdb_comparison(pdb_ref_filename, pdbs_sample_filenames):
    """
    Compare PDB files by plotting the distance of the residues.

    :param pdb_ref_filename: reference PDB filename
    :param pdbs_sample_filenames: list of sample PDB filenames
    """
    pdb_ref = pdbtools.parse_pdb("foo", pdb_ref_filename)
    pdbs_sample = [pdbtools.parse_pdb(i, i) for i in pdbs_sample_filenames]

    sim = RNAPrediction(SysConfig())
    print sim.config["sequence"]
    print sim.config["secstruc"]

    plt.figure(figsize=(14, 8))
    ax = plt.gca()

    plt.title("Residue Distances %s" % (sim.config["name"]))
    plt.xlabel("Residue")
    plt.ylabel(u"Distance to Native Structure / Å")

    x_values = np.arange(1, len(sim.config["sequence"]) + 1)
    x_labels = ["%s\n%s" % ("\n".join(list("%.2d" % (i + 1))), x) for i, x in enumerate(sim.config["sequence"])]
    x_min = 0.5
    x_max = len(x_values) + 0.5
    plt.xlim([x_min, x_max])
    plt.xticks(x_values, x_labels)

    y_min = 0
    y_max = 0  # y limits will be set later when we know the final size

    for pdb_sample in pdbs_sample:

        dists_res, dists_atom, rmsd, rotran = pdbtools.align_structure(pdb_ref, pdb_sample, assign_b_factors=True)
        plt.plot(x_values, dists_res, "-o", label=pdb_sample.id)
        plt.plot((x_min, x_max), (rmsd, rmsd), "--", color=ax.lines[-1].get_color())
        y_max = max(y_max, max(dists_res))

    y_max *= 1.025

    plt.ylim([y_min, y_max])

    hatch = ""
    last_x = ""
    start_i = 0
    for i, x in enumerate(sim.config["secstruc"]):
        if x != last_x and last_x != "":
            if last_x != ".":
                if last_x == "(":
                    hatch = "//"
                elif last_x == ")":
                    hatch = "\\\\"
                plt.fill_between([start_i + 0.5, i + 0.5], y_min, y_max, hatch=hatch, color="black", alpha=0.15)
            start_i = i
        last_x = x

    plt.legend(prop={'size': 12})

    plt.savefig("/tmp/rna_tools_seqdist_%s.png" % os.path.basename(os.getcwd()), bbox_inches="tight")

    plt.show()
