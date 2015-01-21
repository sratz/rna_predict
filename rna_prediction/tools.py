#!/usr/bin/env python2
# coding: utf8

'''
Created on Sep 22, 2014

@author: sra
'''

import Bio.PDB
import itertools
import pprint
import sys
import pickle
import matplotlib.pyplot as plt
import numpy as np
import os.path
from os.path import basename, splitext

from rna_prediction import dcatools
from rna_prediction import pdbtools
from rna_prediction.simulation import parseCstFile
from rna_prediction.simulation import RNAPrediction
from rna_prediction.sysconfig import SysConfig

def tools():
    if sys.argv[1] == "hist":

        distanceMap = dcatools.getContactDistanceMap(westhofVector=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        #distanceMap = pickle.load(open("lol.dat"))

        f, axarr = plt.subplots(4, 4)

        f.suptitle("normed frequency / distance")
        i = 0
        for resPair, distanceMapResPair in sorted(distanceMap.iteritems()):
            print resPair
            dists = []
            for atomPair, distancesAtomPair in distanceMapResPair.iteritems():
                dists += distancesAtomPair

            axarr[i / 4, i % 4].hist(dists, normed=True, bins=50)
            axarr[i / 4, i % 4].set_title(resPair)
            axarr[i / 4, i % 4].set_xlim([0, 25])

            i += 1

        plt.show()

    if sys.argv[1] == "hist2":

        if len(sys.argv) >= 3:
            meanCutoff = float(sys.argv[2])
        else:
            meanCutoff = 6.0

        if len(sys.argv) >= 4:
            stdCutoff = float(sys.argv[3])
        else:
            stdCutoff = 3.0

        distanceMap = dcatools.getContactDistanceMap(westhofVector=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        #distanceMap = pickle.load(open("lol.dat"))
        meanDistanceMap = dcatools.getMeanDistanceMapMean(distanceMap, meanCutoff, stdCutoff)

        f, axarr = plt.subplots(4, 4)

        f.suptitle("mean distance (cutoff: %s), standard deviation (cutoff: %s)" % (meanCutoff, stdCutoff))
        f.subplots_adjust(left=0.02, bottom=0.05, right=0.99, top=0.95, hspace=0.3, wspace=0.1)
        i = 0

        for resPair, meanDistanceMapResPair in sorted(meanDistanceMap.items()):
            print resPair
            data = sorted(meanDistanceMapResPair.items(), key=lambda x: x[1])
            x_names = [d[0] for d in data]
            x = np.arange(len(x_names))
            y = [d[1][0] for d in data]
            y_err = [d[1][1] for d in data]


            axarr[i / 4, i % 4].errorbar(x, y, y_err)
            axarr[i / 4, i % 4].set_xticks(x)
            axarr[i / 4, i % 4].set_xticklabels(x_names, rotation=90, size="x-small")
            axarr[i / 4, i % 4].text(0.05, .95, resPair, transform=axarr[i / 4, i % 4].transAxes, va="top")
            i += 1

        plt.show()


    if sys.argv[1] == "comp2":
        sim = RNAPrediction(SysConfig(), ".")
        models = sim.getModels("100rna_r26_w2", [1, 10, 50, 100, 500], "top")
        pprint.pprint(models)
        print
        models = sim.getModels("100rna_r26_w2", [1, 2, 10], "cluster")
        pprint.pprint(models)
        print
        models = sim.getModels("100rna_r26_w2", ["S_000123_5", "S_000100"], "tag")
        pprint.pprint(models)


    if sys.argv[1] == "comp":
        if sys.argv[2] == "--dca":
            dca_mode = True
            sys.argv.pop(2)
        else:
            dca_mode = False

        comparison_pdb = sys.argv[2]
        sources = sys.argv[3:]

        sim = RNAPrediction(SysConfig(), ".")
        pdb = pdbtools.parsePdb("foo", comparison_pdb)
        chain = pdb[0].child_list[0]

        print "Comparison PDB: %s" % (comparison_pdb)
        plt.figure(figsize=(14,8))

        if dca_mode:
            print "Mode: residue-residue"
            bins = np.linspace(0, 60, 60)
        else:
            print "Mode: atom-atom"
            bins = np.linspace(0, 60, 180)

        for arg in sources:
            print "Processing: %s" % (arg)
            if arg.startswith("filter:"):
                title = arg
                _, dca_file, filtertext = arg.split(":", 2)
                print "  Applying filters: %s" % (filtertext)
                dca = dcatools.parseDcaData(dca_file)
                dcaFilterChain = dcatools.parseFilterLine(filtertext, sim)
                dcatools.filterDcaData(dcaData=dca, dcaFilterChain=dcaFilterChain, quiet=True)

                if not dca_mode:
                    print "  Mapping to atom-atom."
                    distanceMap = dcatools.getContactDistanceMap()
                    distanceMapMean = dcatools.getMeanDistanceMapMean(distanceMap, meanCutoff=6.0, stdCutoff=3.0)
                    cst_info = dcatools.buildCstInfoFromDcaContacts(dca, sequence=sim.config["sequence"], distanceMapMean=distanceMapMean, cstFunction="FADE -100 26 20 -2 2", numberDcaPredictions=100, quiet=True)
            else:
                title = splitext(basename(arg))[0]
                print "  Reading file directly."


                if dca_mode:
                    dca = dcatools.parseDcaData(arg)
                else:
                    cst_info = parseCstFile(sim._parseCstNameAndFilename(arg)[1])

            dists = []
            if dca_mode:
                j = 0
                for d in dca:
                    if not d.useContact:
                        continue
                    j += 1
                    if j > 100:
                        break
                    average_heavy, minimum_heavy, minimum_pair = dcatools.getContactInformationInPdbChain(d, chain)
                    dists.append(minimum_heavy)
            else:
                for cst in cst_info:
                    atom1, res1, atom2, res2, fuction = cst
                    try:
                        dist = np.linalg.norm(chain[res1][atom1].coord - chain[res2][atom2].coord)
                        dists.append(dist)
                        #print atom1, res1, atom2, res2, "   distance:", dist
                    except:
                        print atom1, res1, atom2, res2, "   NOT FOUND"

            print "  Average distance:", np.average(dists)

            plt.hist(dists, bins, alpha=0.6, label=title)
        plt.legend()
        plt.title("%s Constraints Quality %s" % ("Residue-Residue" if dca_mode else "Atom-Atom", os.path.basename(os.getcwd())))
        plt.xlabel(u"Native distance / Å")
        plt.ylabel("Number of contraints")
        plt.savefig("/tmp/rna_tools_quality_%s.png" % os.path.basename(os.getcwd()), bbox_inches="tight")
        plt.show()


    if sys.argv[1] == "checkcontacts":
        dca = dcatools.parseDcaData(dcaPredictionFileName=sys.argv[2])

        # loop over all pdb files
        chains = []
        for i in xrange(3, len(sys.argv)):
            pdb = pdbtools.parsePdb("foo", sys.argv[i])
            chains.append({"chain": pdb[0].child_list[0], "name":sys.argv[i]})


        maximum = min(len(dca), 100)
        averages = np.zeros((len(chains), maximum))
        minimums = np.zeros((len(chains), maximum))

        # loop over all dca entries
        for j, d in enumerate(dca):
            if j >= maximum:
                break

            print d

            for i in xrange(0, len(chains)):

                average_heavy, minimum_heavy, minimum_pair = dcatools.getContactInformationInPdbChain(d, chains[i]["chain"])
                averages[i][j] = average_heavy
                minimums[i][j] = minimum_heavy
                print "        %-70s  avg: %.05f   min: %-3s  %-3s %f" % (chains[i]["name"], average_heavy, minimum_pair[0].name, minimum_pair[1].name, minimum_heavy)


        colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
        x = np.array(range(0, maximum))
        width = 0.8
        width_single = width / len(chains)
        for i in xrange(0, len(chains)):
            color = colors.next()
            plt.bar(x + 1 - width / 2, averages[i], width, label=chains[i]["name"], color=color, alpha=0.5)
        plt.xticks(x + 1, rotation="vertical")
        plt.xlim((0, maximum + 1))
        plt.legend()
        plt.title(sys.argv[2])
        plt.xlabel("DCA Prediction Number")
        plt.ylabel(u"Distance / Å")
        plt.savefig("/tmp/rna_tools_checkcontacts_%s.png" % os.path.basename(os.getcwd()), bbox_inches="tight")
        plt.show()


    # plot native rmsd over rosetta score
    # usage: rna_tools rmsdscore <cst> <max_models> [output_image]
    if sys.argv[1] == "rmsdscore":
        cst = sys.argv[2]
        number = int(sys.argv[3])


        sim = RNAPrediction(SysConfig(), ".")
        cst_name, cst_file = sim._parseCstNameAndFilename(cst)
        models = sim.getModels(cst_name, range(1, number), "top")

        # TODO: this only works for max 10 clusters so far
        colors = ['w', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#ff6600', '#006000', '#600060', '#F7C2CA']

        models_c = []
        plots = []
        descs = []

        plt.figure(figsize=(14,8))

        # cluster0
        models_c.append([m for m in models if "cluster" not in m])
        descs.append("no cluster")

        # cluster1-10
        for i in range(1,11):
            models_c.append([m for m in models if "cluster" in m and m["cluster"] == i])
            descs.append("cluster %d" % i)

        if "native_rmsd" in models_c[1][0]:
            comparison = "native_rmsd"
            plt.xlabel("native rmsd / A")
        else:
            comparison = "rmsd_cluster_1"
            plt.xlabel("rmsd to best structure / A")

        # create plots
        for i in range(0,11):
            plots.append(plt.scatter([x[comparison] for x in models_c[i]], [x["score"] for x in models_c[i]], s=50, c=colors[i]))

        plt.title("model clusters %s, constraints: %s" % (sim.config["name"], cst_name))
        plt.ylabel("rosetta score")
        plt.legend(plots, descs, loc="upper right", prop={'size':12})
        plt.savefig("/tmp/rna_tools_rmsdscore_%s.png" % os.path.basename(os.getcwd()), bbox_inches="tight")
        plt.show()

    # align pdb files and plot distance for each atom/residue in the sequence
    if sys.argv[1] == "seqdist":
        pdb_ref_filename = sys.argv[2]
        pdbs_sample_filenames = sys.argv[3:]
        pdb_ref = pdbtools.parsePdb("foo", pdb_ref_filename)
        pdbs_sample = [pdbtools.parsePdb(i, i) for i in pdbs_sample_filenames]

        sim = RNAPrediction(SysConfig(), ".")
        print sim.config["sequence"]
        print sim.config["secstruc"]

        plt.figure(figsize=(14,8))
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
        y_max = 0 # y limits will be set later when we know the final size

        for pdb_sample in pdbs_sample:

            dists_res, dists_atom, rmsd, rotran = pdbtools.alignStructure(pdb_ref, pdb_sample, assignBFactors=True)
            plt.plot(x_values, dists_res, "-o", label=pdb_sample.id)
            plt.plot((x_min, x_max), (rmsd, rmsd), "--", color=ax.lines[-1].get_color())
            y_max = max(y_max, max(dists_res))

        y_max = y_max * 1.025

        plt.ylim([y_min, y_max])

        hatch = ""
        last_x = ""
        start_i = 0
        for i, x in enumerate(sim.config["secstruc"]):
            if x != last_x and last_x != "" :
                if last_x != ".":
                    if last_x == "(":
                        hatch = "//"
                    elif last_x == ")":
                        hatch = "\\\\"
                    plt.fill_between([start_i + 0.5, i + 0.5], y_min, y_max, hatch=hatch, color="black", alpha=0.15)
                start_i = i
            last_x = x

        plt.legend(prop={'size':12})

        plt.savefig("/tmp/rna_tools_seqdist_%s.png" % os.path.basename(os.getcwd()), bbox_inches="tight")

        plt.show()