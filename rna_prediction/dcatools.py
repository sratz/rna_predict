import os
import numpy
import math
import Bio.PDB
import warnings
from eSBMTools import PdbFile
from Bio.PDB.PDBExceptions import PDBConstructionWarning

'''
Created on Sep 10, 2014

@author: sebastian, blutz
'''
def _getAtomsBackbone(termPhosphate=False):
    atoms = ["P", "OP1", "OP2"] if termPhosphate else []
    return atoms + ["O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]

def getAtomsForRes(res, termPhosphate=False):
    atoms = _getAtomsBackbone(termPhosphate)
    if res == "A":
        return atoms + ["N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"]
    elif res == "U":
        return atoms + ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]
    elif res == "G":
        return atoms + ["N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"]
    elif res == "C":
        return atoms + ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]


# the westhofVector can be used to apply different weights to the bonding family classes
def buildContactDistanceMap(pdbDirectory, structureDirectory, westhofVector=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]):
    # make sure directory names have a trailing slash
    pdbDirectory = os.path.normpath(pdbDirectory) + os.sep
    structureDirectory = os.path.normpath(structureDirectory) + os.sep

    try:
        os.makedirs(pdbDirectory)
    except:
        pass

    # TODO: add exception handling in case structureDirectory does not exist or is missing needed files
    pdbStructureDict = {}
    distanceMap = {}
    for nt1 in ["A", "U", "G", "C"]:
        for nt2 in ["A", "U", "G", "C"]:

            distanceMapResPair = {}

            pdbCodes = []
            residues = []

            for line in open(structureDirectory + nt1 + '-' + nt2 + '.txt'):
                if not line.startswith('-'):
                    pdbCodes.append(line.rstrip('\n').split(' ')[0].upper())
                    residues.append((int(line.rstrip('\n').split(' ')[1]),int(line.rstrip('\n').split(' ')[2])))
                else:
                    pdbCodes.append('-')
                    residues.append('-')

            for pdbCode in pdbCodes:
                if not os.path.exists(pdbDirectory + pdbCode + '.pdb'):
                    if not pdbCode == '-':
                        PdbFile.downloadPdbFile(pdbDirectory, pdbCode)

            for index, pdbCode in enumerate(pdbCodes):
                if pdbCode == '-':
                    continue

                if pdbCode not in pdbStructureDict:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", PDBConstructionWarning)
                        pdbStructureDict[pdbCode] = Bio.PDB.PDBParser().get_structure(pdbCode, pdbDirectory + pdbCode + '.pdb')

                model = pdbStructureDict[pdbCode][0]

                found1 = False
                found2 = False
                for chain in model:
                    try:
                        res1 = chain[residues[index][0]]
                        assert res1.get_resname().strip() == nt1
                        found1 = True
                        if found2:
                            break
                    except:
                        pass

                    try:
                        res2 = chain[residues[index][1]]
                        assert res2.get_resname().strip() == nt2
                        found2 = True
                        if found1:
                            break
                    except:
                        pass

                    if found1 and found2:
                        break

                if not found1 or not found2:
                    raise Exception("Could not find residue contact in pdb file: %s-%s %s %s %s" % (nt1, nt2, pdbCode, residues[index][0], residues[index][1]))

                print "%s-%s %s %s %s" % (nt1, nt2, pdbCode, residues[index][0], residues[index][1])
                for atom1 in res1:
                    for atom2 in res2:
                        if not (atom1.name.startswith('H') or atom2.name.startswith('H')):
                            contactKey = str(atom1.name) + '-' + str(atom2.name)
                            x1, y1, z1 = atom1.coord
                            x2, y2, z2 = atom2.coord
                            distance = westhofVector[index] * math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
                            if not contactKey in distanceMapResPair:
                                distanceMapResPair[contactKey] = [distance]
                            else:
                                distanceMapResPair[contactKey].append(distance)

            distanceMap[nt1 + nt2] = distanceMapResPair
    return distanceMap


def buildMeanDistanceMapMean(distanceMap, meanCutoff=3, stdCutoff=6):
    meanDistanceMap = {}
    for resPair, distanceMapResPair in distanceMap.iteritems():
        meanDistanceMapRes = {}
        for atomPair, distancesAtomPair in distanceMapResPair.iteritems():
            mean = numpy.asarray(distancesAtomPair).mean()
            std = numpy.asarray(distancesAtomPair).std()
            if mean < meanCutoff and std < stdCutoff:
                meanDistanceMapRes[atomPair] = [mean, std]
        meanDistanceMap[resPair] = meanDistanceMapRes
    return meanDistanceMap


def createPdbMapping(sequence, mapping):
    # parse a range in the form of 1-7,80,100-120,8-
    pdbMapping = {}
    i = 1
    ranges = mapping.split(",")
    for r in ranges:
        r = r.split("-")
        if len(r) == 1:
            # single number
            pdbMapping[int(r[0])] = i
            i += 1
        elif r[1] == "":
            # open end range, fill till the end of the sequence
            x = int(r[0])
            while i <= len(sequence):
                pdbMapping[x] = i
                x += 1
                i += 1
        else:
            # regular start-end
            for x in range(int(r[0]), int(r[1]) + 1):
                pdbMapping[x] = i
                i += 1
    return pdbMapping
