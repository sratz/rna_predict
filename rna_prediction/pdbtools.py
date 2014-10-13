'''
Created on Oct 13, 2014

@author: sebastian
'''

import numpy as np
import os
import warnings
import Bio.PDB
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from eSBMTools import PdbFile

from sysconfig import SysConfig


PDB_DIRECTORY = SysConfig.SYSCONFIG_LOCATION + os.sep + "pdbs"


def writePdb(filename, data, model=1, remark=None, append=False):
    with open(filename, "a" if append else "w") as output_pdb_fd:
        output_pdb_fd.write("MODEL %d\n" % (model))
        output_pdb_fd.write("REMARK %s\n" % (remark))
        output_pdb_fd.write(data)
        output_pdb_fd.write("ENDMDL\n")


def extractPOnly(filename):
    def isValidAtom(atom):
        if atom == "P":
            return True
        if atom == "OP1":
            return True
        if atom == "OP2":
            return True
        return False

    p_only = ""
    with open(filename, "r") as fd:
        # only extract first chain
        chain_id = None
        for line in fd:
            fields = line.split()
            if fields[0] == "ATOM" and isValidAtom(fields[2]):
                if chain_id is None:
                    chain_id = fields[4]
                elif chain_id != fields[4]:
                    break
                p_only += line
        p_only += "TER\n"
        return p_only


def getPdbByCode(pdbCode, pdbDirectory=PDB_DIRECTORY):
    # make sure directory names have a trailing slash
    pdbDirectory = os.path.normpath(pdbDirectory) + os.sep

    try:
        os.makedirs(pdbDirectory)
    except:
        pass
    pdbFile = pdbDirectory + pdbCode + '.pdb'
    if not os.path.exists(pdbFile):
        PdbFile.downloadPdbFile(pdbDirectory, pdbCode)
    return parsePdb(pdbCode, pdbFile)


def parsePdb(pdbCode, pdbFile):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        return Bio.PDB.PDBParser().get_structure(pdbCode, pdbFile)


def filterAtoms(atoms, heavyOnly=False, pOnly=False):
    ret = []
    for atom in atoms:
        if pOnly and atom.name != "P":
            continue
        elif heavyOnly and atom.name.startswith("H"):
            continue
        ret.append(atom)
    return ret


def getCenterOfRes(res):
    coords = []
    # loop over all atoms
    for atom in filterAtoms(res, heavyOnly=True):
        coords.append(atom.coord)
    return np.mean(coords, axis=0)
