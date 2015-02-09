"""
Created on Oct 13, 2014

@author: sebastian
"""

import itertools
import numpy as np
import os
import warnings
import Bio.PDB
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from sysconfig import SysConfig
from . import utils


PDB_DIRECTORY = SysConfig.SYSCONFIG_LOCATION + os.sep + "pdbs"


def writePdb(filename, data, model=1, remark=None, append=False):
    with open(filename, "a" if append else "w") as output_pdb_fd:
        output_pdb_fd.write("MODEL %d\n" % model)
        output_pdb_fd.write("REMARK %s\n" % remark)
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
    # only extract first chain
    chain_id = None
    for line in utils.readFileLineByLine(filename):
        fields = line.split()
        if fields[0] == "ATOM" and isValidAtom(fields[2]):
            if chain_id is None:
                chain_id = fields[4]
            elif chain_id != fields[4]:
                break
            p_only += line + "\n"
    p_only += "TER\n"
    return p_only


def downloadPdbFile(pdbCode, pdbDirectory=PDB_DIRECTORY):
    import urllib2

    # make sure directory names have a trailing slash
    pdbDirectory = os.path.normpath(pdbDirectory) + os.sep

    response = urllib2.urlopen("http://www.rcsb.org/pdb/files/%s.pdb" % pdbCode.upper())
    with open(pdbDirectory + pdbCode + ".pdb", "w") as f:
        f.write(response.read())


def getPdbByCode(pdbCode, pdbDirectory=PDB_DIRECTORY):
    # make sure directory names have a trailing slash
    pdbDirectory = os.path.normpath(pdbDirectory) + os.sep

    utils.mkdir_p(pdbDirectory)
    pdbFile = pdbDirectory + pdbCode + '.pdb'
    if not os.path.exists(pdbFile):
        downloadPdbFile(pdbCode, pdbDirectory)
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


def alignStructure(refPdb, movingPdb, assignBFactors=True):
    """Aligns movingPdb to refPdb. Returns (resDists, atomDists, rmsd, transformationMatrix)"""
    chain_ref = refPdb[0].child_list[0]
    chain_sample = movingPdb[0].child_list[0]

    ref_atoms = []
    ref_res = []
    sample_atoms = []
    sample_res = []

    for res1 in chain_ref:
        if res1.id not in chain_sample:
            print "skipping %s" % res1
            continue

        res2 = chain_sample[res1.id]

        ref_res.append(res1)
        sample_res.append(res2)

        for atom1 in res1:
            if atom1.id not in res2:
                print "  skipping %s" % atom1
                continue

            atom2 = res2[atom1.id]

            ref_atoms.append(atom1)
            sample_atoms.append(atom2)

    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(movingPdb)

    # in case we want to assign b-factors to the moving structure, set them all to 0 before.
    if assignBFactors:
        for res in chain_sample:
            for atom in res:
                atom.set_bfactor(0)

    dists_atom = []
    for atom1, atom2 in itertools.izip(ref_atoms, sample_atoms):
        dist = np.linalg.norm(atom1.coord - atom2.coord)
        dists_atom.append(dist)
        if assignBFactors:
            atom2.set_bfactor(dist)

    print super_imposer.rotran
    print super_imposer.rms

    dists_res = []
    for res1, res2 in itertools.izip(ref_res, sample_res):
        dist = np.linalg.norm(getCenterOfRes(res1) - getCenterOfRes(res2))
        dists_res.append(dist)

    return dists_res, dists_atom, super_imposer.rms, super_imposer.rotran
