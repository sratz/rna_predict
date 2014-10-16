import glob
import os
import pickle
import numpy as np
import math
import re

from . import pdbtools
from sysconfig import SysConfig


'''
Created on Sep 10, 2014

@author: sebastian, blutz
'''

INFO_DIRECTORY = SysConfig.SYSCONFIG_LOCATION + os.sep + "structure_info"
CACHE_DIRECTORY = SysConfig.SYSCONFIG_LOCATION + os.sep + "cache"
CACHE_DISTANCEMAP = CACHE_DIRECTORY + os.sep + "distancemap.dat"


class DcaException(Exception):
    pass


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
def getContactDistanceMap(structureDirectory=INFO_DIRECTORY, westhofVector=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], forceRebuild=False):
    # try to use a cached version of the distance map if found and recent and forceRebuild is False
    if not forceRebuild:
        try:
            cacheOk = True
            if os.path.isfile(CACHE_DISTANCEMAP):
                cacheTimestamp = os.path.getmtime(CACHE_DISTANCEMAP)
                for d in glob.glob(INFO_DIRECTORY + os.sep + "*.txt"):
                    if (os.path.getmtime(d) > cacheTimestamp):
                        cacheOk = False
                        print "Contact map cache out of date. Rebuilding..."
                        break
                if cacheOk:
                    with open(CACHE_DISTANCEMAP, "r") as f:
                        return pickle.load(f)
        except Exception:
            print "Contact map cache broken. Rebuilding..."

    print "Building contact distance map:"

    # make sure directory names have a trailing slash
    structureDirectory = os.path.normpath(structureDirectory) + os.sep

    # TODO: add exception handling in case structureDirectory does not exist or is missing needed files
    pdbStructureDict = {}
    distanceMap = {}
    for nt1 in ["A", "U", "G", "C"]:
        for nt2 in ["A", "U", "G", "C"]:

            distanceMapResPair = {}

            pdbCodes = []
            residues = []

            # read the structures for the 12 edge-to-edge bonding families
            for line in open(structureDirectory + nt1 + '-' + nt2 + '.txt'):
                fields = line.strip().split(" ")
                if fields[0] != "-":
                    pdbCodes.append(fields[0].upper())
                    residues.append((int(fields[1]), int(fields[2])))
                else:
                    pdbCodes.append(None)
                    residues.append(None)

            for index, pdbCode in enumerate(pdbCodes):
                if pdbCode is None:
                    continue

                if pdbCode not in pdbStructureDict:
                    pdbStructureDict[pdbCode] = pdbtools.getPdbByCode(pdbCode)

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

    # save distance map in cache
    try:
        os.makedirs(CACHE_DIRECTORY)
    except:
        pass
    with open(CACHE_DISTANCEMAP, "w") as f:
        pickle.dump(distanceMap, f)

    return distanceMap


def getMeanDistanceMapMean(distanceMap, meanCutoff=None, stdCutoff=None):
    meanDistanceMap = {}
    for resPair, distanceMapResPair in distanceMap.iteritems():
        meanDistanceMapRes = {}
        for atomPair, distancesAtomPair in distanceMapResPair.iteritems():
            mean = np.asarray(distancesAtomPair).mean()
            std = np.asarray(distancesAtomPair).std()
            if (meanCutoff is None or mean < meanCutoff) and (stdCutoff is None or std < stdCutoff):
                meanDistanceMapRes[atomPair] = [mean, std]
        meanDistanceMap[resPair] = meanDistanceMapRes
    return meanDistanceMap


def createPdbMapping(mapping):
    # parse a range in the form of 1-7,80,100-120,8-9
    try:
        pdbMapping = {}
        i = 1
        ranges = mapping.split(",")
        for r in ranges:
            r = r.split("-")
            if len(r) == 1:
                # single number
                pdbMapping[int(r[0])] = i
                i += 1
            else:
                # regular start-end
                for x in range(int(r[0]), int(r[1]) + 1):
                    pdbMapping[x] = i
                    i += 1
        return pdbMapping
    except:
        raise DcaException("Invalid pdb mapping string: %s" % (mapping))


def parseDcaData(dcaPredictionFileName, pdbMappingOverride=None):
    pdbMapping = None
    if pdbMappingOverride is not None:
        print "  %s pdbMapping (user): %s" % (" " if pdbMappingOverride is None else "*", pdbMappingOverride)
        pdbMapping = createPdbMapping(pdbMappingOverride)
    pattern_parameter = re.compile(r"^#\s(\S+)\s+(.*)$")
    dca = []
    with open(dcaPredictionFileName) as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == "#":
                # comment / parameter line
                m = pattern_parameter.match(line)
                if m:
                    if m.group(1) == "pdb-mapping":
                        # if pdbMapping is already set, it was overridden on invocation, skip this line!
                        if pdbMapping is not None:
                            continue
                        print "  %s pdbMapping (file): %s" % (" " if pdbMapping is not None else "*", m.group(2))
                        pdbMapping = createPdbMapping(m.group(2))
                continue
            # data line
            parts = line.split(" ")
            if pdbMapping is not None:
                try:
                    dca.append([pdbMapping[int(parts[0])], pdbMapping[int(parts[1])]])
                except:
                    raise DcaException("Invalid PDB mapping. Could not access residue: %s" % (parts[:2]))
            else:
                dca.append([int(parts[0]), int(parts[1])])
    return dca


# return True if contact is not realized
# TODO: implement more stuff parameters?
def filterOutDcaContact(contact, pdbChain, threshold):
    average_heavy, minimum_heavy, minimum_pair = getContactInformationInPdbChain(contact, pdbChain)
    print "Threshold filter: min_distance: %f, threshold: %f" % (minimum_heavy, threshold)
    # filter out contact if realized minimum distance is larger than the threshold
    if minimum_heavy > threshold:
        return True
    return False


# returns distance information about dca contact in a realized pdb chain
def getContactInformationInPdbChain(dcaContact, pdbChain):
    res1, res2 = (pdbChain[dcaContact[0]], pdbChain[dcaContact[1]])
    # calculate average distance
    average_heavy = np.linalg.norm(pdbtools.getCenterOfRes(res1) - pdbtools.getCenterOfRes(res2))

    minimum_heavy = 9999
    minimum_pair = []
    for atom1 in pdbtools.filterAtoms(res1, heavyOnly=True):
        for atom2 in pdbtools.filterAtoms(res2, heavyOnly=True):
            dist = np.linalg.norm(atom1.coord - atom2.coord)
            if dist < minimum_heavy:
                minimum_heavy = dist
                minimum_pair = [atom1, atom2]
    return (average_heavy, minimum_heavy, minimum_pair)
