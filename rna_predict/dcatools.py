# coding: utf-8

import os
from pkg_resources import ResourceManager
try:
    import cPickle as pickle
except ImportError:
    import pickle
import numpy as np
import re

from . import pdbtools
from . import utils
from .sysconfig import SysConfig


INFO_DIRECTORY = SysConfig.SYSCONFIG_LOCATION + os.sep + "structure_info"
CACHE_DIRECTORY = SysConfig.SYSCONFIG_LOCATION + os.sep + "cache"
CACHE_DISTANCEMAP = CACHE_DIRECTORY + os.sep + "distancemap.dat"


class DcaException(Exception):
    """
    Custom exception class used for foreseeable, DCA related errors.
    """
    pass


def _get_atoms_backbone(term_phosphate=False):
    """
    Get list of backbone atoms.

    :param term_phosphate: add P atoms
    :return: list of atoms
    """
    atoms = ["P", "OP1", "OP2"] if term_phosphate else []
    return atoms + ["O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]


def get_atoms_for_res(res, term_phosphate=False):
    """
    Get list of atoms for residue.

    :param res: nucleotide (A,U,G,C)
    :param term_phosphate: add P atoms
    :return: list of atoms
    """
    atoms = _get_atoms_backbone(term_phosphate)
    if res == "A":
        return atoms + ["N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"]
    elif res == "U":
        return atoms + ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]
    elif res == "G":
        return atoms + ["N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"]
    elif res == "C":
        return atoms + ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]


def get_atoms_for_res_sequence(sequence):
    """
    Get list of atoms for a sequence of nucleotides

    :param sequence: sequence as text
    :return: list of atoms
    """
    atoms = []
    first = True
    for res in sequence:
        res = res.upper()
        atoms.append((res, get_atoms_for_res(res, term_phosphate=(not first))))
        first = False
    return atoms


def get_contact_distance_map(structure_directory=INFO_DIRECTORY, westhof_vector=None, force_rebuild=False):
    """
    Returns contact distance map

    The contact distance map is cached it in the user directory and updated when newer files are found.

    :param structure_directory: directory to look up structure information text files
    :param westhof_vector: list of factors to apply different weights to the bonding family classes (defaults to ``[1, 1, ... ]``)
    :param force_rebuild: force rebuilding the distance map
    """
    # default: same weight for all families
    if not westhof_vector:
        westhof_vector = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    nucleotides = ["A", "U", "G", "C"]

    # build a dict of filenames
    # if a local file is present in the user directory it will take precedence over the system wide shipped version
    structure_filenames = {}
    resource_manager = ResourceManager()
    for nt1 in nucleotides:
        for nt2 in nucleotides:
            ntpair = "%s-%s" % (nt1, nt2)
            local_file = structure_directory + os.sep + ntpair + ".txt"
            if os.path.isfile(local_file):
                structure_filenames[ntpair] = local_file
            else:
                structure_filenames[ntpair] = resource_manager.resource_filename(__name__, "structure_info/%s.txt" % ntpair)

    # try to use a cached version of the distance map if found and recent and force_rebuild is False
    if not force_rebuild:
        try:
            cache_ok = True
            if os.path.isfile(CACHE_DISTANCEMAP):
                cache_timestamp = os.path.getmtime(CACHE_DISTANCEMAP)
                for d in structure_filenames.itervalues():
                    if os.path.getmtime(d) > cache_timestamp:
                        cache_ok = False
                        print "Contact map cache out of date. Rebuilding..."
                        break
                if cache_ok:
                    with open(CACHE_DISTANCEMAP, "r") as f:
                        return pickle.load(f)
        except (IOError, pickle.PickleError, AttributeError, EOFError, IndexError):
            print "Contact map cache broken. Rebuilding..."

    print "Building contact distance map:"

    pdb_structure_dict = {}
    distance_map = {}
    for nt1 in nucleotides:
        for nt2 in nucleotides:

            distance_map_res_pair = {}

            pdb_codes = []
            residues = []

            # read the structures for the 12 edge-to-edge bonding families
            for line in utils.read_file_line_by_line(structure_filenames[nt1 + '-' + nt2]):
                fields = line.split(" ")
                if fields[0] != "-":
                    pdb_codes.append(fields[0].upper())
                    residues.append((int(fields[1]), int(fields[2])))
                else:
                    pdb_codes.append(None)
                    residues.append(None)

            # loop over all pdbcodes and their index in the list (0-11)
            for index, pdb_code in enumerate(pdb_codes):
                # skip if we don't have any entry for this family
                if pdb_code is None:
                    continue

                # download pdb if necessary
                if pdb_code not in pdb_structure_dict:
                    pdb_structure_dict[pdb_code] = pdbtools.parse_pdb(pdb_code, pdbtools.get_pdb_by_code(pdb_code))

                # extract model from pdb
                model = pdb_structure_dict[pdb_code][0]

                # try to find the residue contact specified. this is done by looping over all chains in the model,
                # and checking if the residue is in there and is the correct nucleotide
                def find_res(res, resname):
                    for chain in model:
                        try:
                            if chain[res].get_resname().strip() == resname:
                                return chain[res]
                        except KeyError:
                            pass
                    return None

                res1 = find_res(residues[index][0], nt1)
                res2 = find_res(residues[index][1], nt2)

                if not res1 or not res2:
                    raise Exception("Could not find residue contact in pdb file: %s-%s %s %s %s" % (nt1, nt2, pdb_code, residues[index][0], residues[index][1]))

                print "%s-%s %s %s %s" % (nt1, nt2, pdb_code, residues[index][0], residues[index][1])

                # add all atom-atom contacts to the distance map for the current residue pair
                for atom1 in res1:
                    for atom2 in res2:
                        if not (atom1.name.startswith('H') or atom2.name.startswith('H')):
                            contact_key = str(atom1.name) + '-' + str(atom2.name)
                            distance = westhof_vector[index] * (atom1 - atom2)
                            if contact_key not in distance_map_res_pair:
                                distance_map_res_pair[contact_key] = [distance]
                            else:
                                distance_map_res_pair[contact_key].append(distance)

            distance_map[nt1 + nt2] = distance_map_res_pair

    # save distance map in cache
    utils.mkdir_p(CACHE_DIRECTORY)
    with open(CACHE_DISTANCEMAP, "w") as f:
        pickle.dump(distance_map, f)

    return distance_map


def get_contact_distance_map_mean(distance_map, mean_cutoff=None, std_cutoff=None):
    """
    Return an average distance map containing only those contacts with average distance and standard deviation satisfiying a cutoff.

    :param distance_map: full distance map
    :param mean_cutoff: limit for average
    :param std_cutoff: limit for standard deviation
    :return: average distance map
    """
    mean_distance_map = {}
    for res_pair, distance_map_res_pair in distance_map.iteritems():
        mean_distance_map_res = {}
        for atom_pair, distances_atom_pair in distance_map_res_pair.iteritems():
            mean = np.asarray(distances_atom_pair).mean()
            std = np.asarray(distances_atom_pair).std()
            if (mean_cutoff is None or mean < mean_cutoff) and (std_cutoff is None or std < std_cutoff):
                mean_distance_map_res[atom_pair] = [mean, std]
        mean_distance_map[res_pair] = mean_distance_map_res
    return mean_distance_map


def read_pdb_mapping_from_file(dca_prediction_filename):
    """
    Read a PDB mapping from DCA file if present and return it as text

    :param dca_prediction_filename: DCA input filename
    :return: PDB mapping text
    """
    pattern_parameter = re.compile(r"^#\s(\S+)\s+(.*)$")
    for line in utils.read_file_line_by_line(dca_prediction_filename):
        if line[0] == "#":
            # comment / parameter line
            m = pattern_parameter.match(line)
            if m:
                if m.group(1) == "pdb-mapping":
                    return m.group(2)
            continue
        # we only allow comments on top, so no need to check the rest of the lines:
        break
    return None


def parse_dca_data(dca_prediction_filename):
    """
    Read a DCA file, adjust the sequence numbers to match the alignment of the PDB file, and create a list of DcaContacts

    :param dca_prediction_filename: DCA input filename
    :return: list of DcaContact objects
    """
    print "Parsing dca file %s..." % dca_prediction_filename

    # read pdb mapping from file
    pdb_mapping_text = read_pdb_mapping_from_file(dca_prediction_filename)
    if pdb_mapping_text is None:
        print "pdb-mapping: 1-N (no pdb-mapping found in header of dca file)"
    else:
        print "pdb-mapping: %s" % pdb_mapping_text

    # parse pdb mapping to dict
    # for example 1-7,80,100-120,8-9 --> {1: 1, 2: 2, ..., 80: 8, 81: 9, ...}
    try:
        pdb_mapping = None if pdb_mapping_text is None else {x: i + 1 for i, x in enumerate(utils.comma_separated_ranges_to_list(pdb_mapping_text))}
    except ValueError as e:
        raise DcaException("Invalid pdb mapping string: %s" % e.message)

    dca = []
    for line in utils.read_file_line_by_line(dca_prediction_filename):
        if line[0] == "#":
            continue
        # data line
        parts = line.split(" ")
        if pdb_mapping is not None:
            try:
                res1 = pdb_mapping[int(parts[0])]
                res2 = pdb_mapping[int(parts[1])]
            except KeyError:
                # our dca file might contain more than what we are interested in for the current prediction.
                # simply ignore all dca pairs that don't fit the mapping range, but print a warning.
                print "Warning: At least one residue of contact %s not found in PDB-Mapping. Ignoring..." % parts[:2]
                continue
        else:
            res1 = int(parts[0])
            res2 = int(parts[1])

        dca.append(DcaContact(res1, res2))
    return dca


def get_contact_information_in_pdb_chain(dca_contact, pdb_chain, heavy_only=True):
    """
    Returns distance information about a DCA contact in a realized PDB chain

    Return value is a tuple of:

    - Average distance: Mean distance of all atoms in the contacts.
    - Minimum distance: Minimum distance between two atoms in the contact.
    - Minimum pair: List of ``[atom1, atom2]`` forming the minimal contact

    :param dca_contact: DcaContact object
    :param pdb_chain: PDB chain structure object
    :param heavy_only: Only use heavy atoms
    :return: tuple ``(average_dist, minimum_dist, minimum_pair)``. In case the contact cannot be found in the PDB file ``(0, 0, None)`` is returned.
    """
    try:
        res1, res2 = (pdb_chain[dca_contact.res1], pdb_chain[dca_contact.res2])
    except KeyError:
        return 0, 0, None
    # calculate average distance
    average_heavy = np.linalg.norm(pdbtools.get_center_of_res(res1) - pdbtools.get_center_of_res(res2))

    minimum_heavy = 9999
    minimum_pair = []
    for atom1 in pdbtools.filter_atoms(res1, heavy_only):
        for atom2 in pdbtools.filter_atoms(res2, heavy_only):
            dist = np.linalg.norm(atom1.coord - atom2.coord)
            if dist < minimum_heavy:
                minimum_heavy = dist
                minimum_pair = [atom1, atom2]
    return average_heavy, minimum_heavy, minimum_pair


def build_cst_info_from_dca_contacts(dca_data, sequence, mapping_mode, cst_function, number_dca_predictions, quiet=False):
    """
    Maps DCA residue contacts to atom-atom constraints.

    :param dca_data: list od DcaContacts
    :param sequence: sequence as text
    :param mapping_mode: atom-to-atom mapping mode to use, supported values: "minAtom" or "pOnly"
    :param cst_function: rosetta function and parameters as text string
    :param number_dca_predictions: maximum number of DCA predictions to use
    :param quiet: reduce output verbosity
    :return: list of constraint information
    """
    mapping_mode = mapping_mode.lower()
    if mapping_mode not in ["minAtom", "ponly"]:
        raise DcaException("build_cst_info: Invalid mapping mode given: %s" % mapping_mode)

    if mapping_mode == "minAtom":
        # load contact map for atom-atom contacts
        distance_map = get_contact_distance_map()
        distance_map_mean = get_contact_distance_map_mean(distance_map, mean_cutoff=6.0, std_cutoff=3.0)

    atoms = get_atoms_for_res_sequence(sequence)

    cst_info = []
    predictions_used = 0
    for i, d in enumerate(dca_data):
        if predictions_used >= number_dca_predictions:
            if not quiet:
                print "Limit of %d used predictions reached. Stopping..." % number_dca_predictions
            break

        # print some information about the dca contact
        if not quiet:
            print "Contact %d: %s" % (i + 1, d)

        # skip contact completely?
        if not d.use_contact:
            if not quiet:
                print "  Dca contact skipped."
            continue

        # lookup contact
        try:
            res1 = atoms[d.res1 - 1]
            res2 = atoms[d.res2 - 1]
            contact_key = res1[0] + res2[0]
        except IndexError:
            # This happens only when there was no pdb-mapping specified in the dca file header.
            # TODO: Maybe this whole mapping parsing thing should be worked over so that these
            # two cases can be unified.
            print "  Warning: Residues of contact cannot be found in current sequence. Ignoring..."
            continue

        predictions_used += 1
        if not quiet:
            print "  Dca contact used (%d)." % predictions_used

        # build atom-atom constraints
        if mapping_mode == "minAtom":
            for atom1 in res1[1]:
                for atom2 in res2[1]:
                    atom_contact_key = atom1 + '-' + atom2
                    if atom_contact_key in distance_map_mean[contact_key]:
                        distance = distance_map_mean[contact_key][atom_contact_key][0] / 10.0
                        if not quiet:
                            print "[%s, %s] %s %s %s" % (d.res1, d.res2, contact_key, atom_contact_key, distance)
                        cst_info.append([atom1, d.res1, atom2, d.res2, d.get_rosetta_function(cst_function)])
        elif mapping_mode == "ponly":
            # TODO: rosetta does not add P atoms for the first residue, so just skip those?
            if d.res1 == 1 or d.res2 == 1:
                continue
            atom_contact_key = "P-P"
            if not quiet:
                print "[%s, %s] %s %s" % (d.res1, d.res2, contact_key, atom_contact_key)
            cst_info.append(["P", d.res1, "P", d.res2, d.get_rosetta_function(cst_function)])
    return cst_info


class DcaContact(object):
    """Class representing a DCA contact."""
    def __init__(self, res1, res2, use_contact=True, weight=1):
        """
        Create new DCA contact.

        :param res1: number of first residue (from 1)
        :param res2: number of second residue (from 1)
        :param use_contact: True or False if contact is to be used
        :param weight: assign a different weight to the contact (default 1)
        """
        self.res1 = res1
        self.res2 = res2
        self.use_contact = use_contact
        self.weight = weight

    def __str__(self):
        return "[%s, %s], use_contact=%s, weight=%f" % (self.res1, self.res2, self.use_contact, self.weight)

    def get_rosetta_function(self, function="FADE -100 26 20 -2 2"):
        """
        Return Rosetta function of the contact as a list while applying weight.

        :param function: rosetta function as text
        :return: rosetta function as a list with applied weight
        """
        function = function.split()

        # parse numeric arguments
        # noinspection PyShadowingNames
        def parse_function(function, types):
            for arg, argtype in types:
                try:
                    function[arg] = argtype(function[arg])
                except ValueError:
                    raise DcaException("Invalid Rosetta function: Could not parse argument '%s' as '%s'" % (function[arg], argtype.__name__))

        try:
            if function[0] == "FADE":
                parse_function(function, [(1, int), (2, int), (3, int), (4, float), (5, float)])  # TODO: make all of these floats?
                return [function[0], function[1], function[2], function[3], function[4] * self.weight, function[5] * self.weight]
            else:
                raise DcaException("Invalid Rosetta function: '%s' is not implemented! Only 'FADE' function is recognized at this time." % function[0])
        except IndexError:
            raise DcaException("Invalid Rosetta function: Not enough arguments")


# DCA FILTERING

def filter_dca_data(dca_data, dca_filter_chain, quiet=False):
    """
    Run list of DCA contacts through a chain of filters.

    :param dca_data: list of DcaContact objects
    :param dca_filter_chain: list of DcaFilter objects
    :param quiet: reduce output verbosity
    """
    if dca_filter_chain is None:
        return dca_data
    for dca_filter in dca_filter_chain:
        if dca_filter is not None:
            for contact in dca_data:
                dca_filter.apply(contact, quiet)


class DcaFilter(object):
    """Filter base class."""

    def apply(self, contact, quiet=False):
        """Apply filter to a DCA contact.

        :param contact: DCA contact
        :param quiet: reduce verbosity
        """
        raise NotImplementedError


class DcaFilterThreshold(DcaFilter):
    """Filter to skips DCA contact if below or above a threshold."""

    def __init__(self, pdb_chain, threshold, keep_below=True, mode="minimum_heavy"):
        """
        Create a new threshold filter.

        :param pdb_chain: PDB chain
        :param threshold: threshold below or above to keep a contact
        :param keep_below: True to keep below threshold, False to keep above
        :param mode: What distance to compare to (average_heavy, minimum_heavy)
        """
        self.pdb_chain = pdb_chain
        self.threshold = threshold
        self.keep_below = keep_below
        self.mode = mode

    def __repr__(self):
        return "Threshold filter: keep %s %f (%s)" % ("below" if self.keep_below else "above", self.threshold, self.mode)

    def apply(self, contact, quiet=False):
        # do not touch contacts that are already disabled
        if not contact.use_contact:
            return

        # get contact information
        average_heavy, minimum_heavy, minimum_pair = get_contact_information_in_pdb_chain(contact, self.pdb_chain)

        # skip if not found
        if minimum_pair is None:
            contact.use_contact = False
            print "%s: not found --> skip"
            return

        # which distance should be used?
        if self.mode == "average_heavy":
            distance = average_heavy
        else:
            distance = minimum_heavy

        skip = False
        # set contact to disabled if filter failed
        if (self.keep_below and distance >= self.threshold) or (not self.keep_below and distance <= self.threshold):
            skip = True
            contact.use_contact = False

        if not quiet:
            print "%s: %f --> %s" % (self, distance, "skip" if skip else "keep")
