# coding: utf-8

import itertools
import numpy as np
import os
import urllib2
import warnings

import Bio.PDB
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from . import utils
from .sysconfig import SysConfig


PDB_DIRECTORY = SysConfig.SYSCONFIG_LOCATION + os.sep + "pdbs"


def get_pdb_by_code(pdb_code, pdb_directory=PDB_DIRECTORY):
    """
    Get PDB file by code. Download if necessary.

    :param pdb_code: PDB code to download
    :param pdb_directory: directory lookup and store the PDB file, defaults to the rna_prediction PDB cache directory
    :return: PDB filename
    """

    # make sure directory names have a trailing slash
    pdb_directory = os.path.normpath(pdb_directory) + os.sep

    utils.mkdir_p(pdb_directory)
    pdb_file = pdb_directory + pdb_code + '.pdb'
    if not os.path.exists(pdb_file):
        response = urllib2.urlopen("http://www.rcsb.org/pdb/files/%s.pdb" % pdb_code.upper())
        with open(pdb_file, "w") as f:
            f.write(response.read())
    return pdb_file


def parse_pdb(pdb_code, pdb_file):
    """
    Parse PDB file using Biopyhon.
    :param pdb_code: internal id
    :param pdb_file: PDB filename
    :return: PDB structure object
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        return Bio.PDB.PDBParser().get_structure(pdb_code, pdb_file)


def filter_atoms(atoms, heavy_only=False, p_only=False):
    """
    Filter list of atoms.
    :param atoms: list of atoms
    :param heavy_only: only keep heavy atoms
    :param p_only: only keep P atoms
    :return: filtered list of atoms
    """
    ret = []
    for atom in atoms:
        if p_only and atom.name != "P":
            continue
        elif heavy_only and atom.name.startswith("H"):
            continue
        ret.append(atom)
    return ret


def get_center_of_res(res):
    """
    Calculate the center of a residue.
    :param res: residue object
    :return: center coordinates
    """
    coords = []
    # loop over all atoms
    for atom in filter_atoms(res, heavy_only=True):
        coords.append(atom.coord)
    return np.mean(coords, axis=0)


def align_structure(ref_pdb, moving_pdb, assign_b_factors=True):
    """
    Align one PDB structure to another.

    :param ref_pdb: reference PDB structure object
    :param moving_pdb: moving PDB structure object
    :param assign_b_factors: place distance values in the b-factor field
    :return: tuple (res_dists, atom_dists, rmsd, transformation_matrix)
    """
    chain_ref = ref_pdb[0].child_list[0]
    chain_sample = moving_pdb[0].child_list[0]

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
    super_imposer.apply(moving_pdb)

    # in case we want to assign b-factors to the moving structure, set them all to 0 before.
    if assign_b_factors:
        for res in chain_sample:
            for atom in res:
                atom.set_bfactor(0)

    dists_atom = []
    for atom1, atom2 in itertools.izip(ref_atoms, sample_atoms):
        dist = np.linalg.norm(atom1.coord - atom2.coord)
        dists_atom.append(dist)
        if assign_b_factors:
            atom2.set_bfactor(dist)

    print super_imposer.rotran
    print super_imposer.rms

    dists_res = []
    for res1, res2 in itertools.izip(ref_res, sample_res):
        dist = np.linalg.norm(get_center_of_res(res1) - get_center_of_res(res2))
        dists_res.append(dist)

    return dists_res, dists_atom, super_imposer.rms, super_imposer.rotran
