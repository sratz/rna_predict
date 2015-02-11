import os
import pickle
import numpy as np
import re
from pkg_resources import ResourceManager

from . import pdbtools
from . import utils
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


def _get_atoms_backbone(term_phosphate=False):
    atoms = ["P", "OP1", "OP2"] if term_phosphate else []
    return atoms + ["O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]


def get_atoms_for_res(res, term_phosphate=False):
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
    atoms = []
    first = True
    for res in sequence:
        res = res.upper()
        atoms.append((res, get_atoms_for_res(res, term_phosphate=(not first))))
        first = False
    return atoms


# the westhofVector can be used to apply different weights to the bonding family classes
def get_contact_distance_map(structure_directory=INFO_DIRECTORY, westhof_vector=None, force_rebuild=False):
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
            for index, pdbCode in enumerate(pdb_codes):
                # skip if we don't have any entry for this family
                if pdbCode is None:
                    continue

                # download pdb if necessary
                if pdbCode not in pdb_structure_dict:
                    pdb_structure_dict[pdbCode] = pdbtools.get_pdb_by_code(pdbCode)

                # extract model from pdb
                model = pdb_structure_dict[pdbCode][0]

                # try to find the residue contact specified. this is done by looping over all chains in the model,
                # and checking if the residue is in there and is the correct nucleotide
                res1 = None
                for chain in model:
                    try:
                        res1 = chain[residues[index][0]]
                        assert res1.get_resname().strip() == nt1
                        break
                    except (KeyError, AssertionError):
                        res1 = None

                res2 = None
                for chain in model:
                    try:
                        res2 = chain[residues[index][1]]
                        assert res2.get_resname().strip() == nt2
                        break
                    except (KeyError, AssertionError):
                        res2 = None

                if not res1 or not res2:
                    raise Exception("Could not find residue contact in pdb file: %s-%s %s %s %s" % (nt1, nt2, pdbCode, residues[index][0], residues[index][1]))

                print "%s-%s %s %s %s" % (nt1, nt2, pdbCode, residues[index][0], residues[index][1])

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
    mean_distance_map = {}
    for resPair, distanceMapResPair in distance_map.iteritems():
        mean_distance_map_res = {}
        for atomPair, distancesAtomPair in distanceMapResPair.iteritems():
            mean = np.asarray(distancesAtomPair).mean()
            std = np.asarray(distancesAtomPair).std()
            if (mean_cutoff is None or mean < mean_cutoff) and (std_cutoff is None or std < std_cutoff):
                mean_distance_map_res[atomPair] = [mean, std]
        mean_distance_map[resPair] = mean_distance_map_res
    return mean_distance_map


def create_pdb_mapping_from_string(mapping):
    # parse a range in the form of 1-7,80,100-120,8-9
    if mapping is None or mapping == "":
        return None
    try:
        pdb_mapping = {}
        i = 1
        ranges = mapping.split(",")
        for r in ranges:
            r = r.split("-")
            if len(r) == 1:
                # single number
                pdb_mapping[int(r[0])] = i
                i += 1
            else:
                # regular start-end
                for x in range(int(r[0]), int(r[1]) + 1):
                    pdb_mapping[x] = i
                    i += 1
        return pdb_mapping
    except:
        raise DcaException("Invalid pdb mapping string: %s" % mapping)


# read a pdb mapping from dca file if present and return it as text
def read_pdb_mapping_from_file(dca_prediction_filename):
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


# reads a dca file and adjusts the sequence numbers to match the alignment of the pdb file
def parse_dca_data(dca_prediction_filename):
    print "Parsing dca file %s..." % dca_prediction_filename
    # read pdb mapping from file
    pdb_apping_text = read_pdb_mapping_from_file(dca_prediction_filename)
    pdb_apping = create_pdb_mapping_from_string(pdb_apping_text)
    print "pdb-mapping: %s" % pdb_apping_text

    if pdb_apping_text is None:
        print "Warning: no pdb-mapping found in header of dca file %s, assuming '1-N'" % dca_prediction_filename

    dca = []
    for line in utils.read_file_line_by_line(dca_prediction_filename):
        if line[0] == "#":
            continue
        # data line
        parts = line.split(" ")
        if pdb_apping is not None:
            try:
                res1 = pdb_apping[int(parts[0])]
                res2 = pdb_apping[int(parts[1])]
            except:
                raise DcaException("Invalid PDB mapping. Could not access residue: %s" % (parts[:2]))
        else:
            res1 = int(parts[0])
            res2 = int(parts[1])

        dca.append(DcaContact(res1, res2))
    return dca


# returns distance information about dca contact in a realized pdb chain
def get_contact_information_in_pdb_chain(dca_contact, pdb_chain):
    res1, res2 = (pdb_chain[dca_contact.res1], pdb_chain[dca_contact.res2])
    # calculate average distance
    average_heavy = np.linalg.norm(pdbtools.get_center_of_res(res1) - pdbtools.get_center_of_res(res2))

    minimum_heavy = 9999
    minimum_pair = []
    for atom1 in pdbtools.filter_atoms(res1, heavy_only=True):
        for atom2 in pdbtools.filter_atoms(res2, heavy_only=True):
            dist = np.linalg.norm(atom1.coord - atom2.coord)
            if dist < minimum_heavy:
                minimum_heavy = dist
                minimum_pair = [atom1, atom2]
    return average_heavy, minimum_heavy, minimum_pair


# maps dca residue contacts to atom-atom constraints
# mode: mapping mode, can be: allAtomWesthof,pOnly
def build_cst_info_from_dca_contacts(dca_data, sequence, mapping_mode, cst_function, number_dca_predictions, quiet=False):
    mapping_mode = mapping_mode.lower()
    if mapping_mode not in ["allatomwesthof", "ponly"]:
        raise DcaException("buildCstInfo: Invalid mapping mode given: %s" % mapping_mode)

    if mapping_mode == "allatomwesthof":
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
        if not d.useContact:
            if not quiet:
                print "  Dca contact skipped."
            continue
        if not quiet:
            print "  Dca contact used (%d)." % (predictions_used + 1)
        predictions_used += 1

        # build atom-atom constraints
        res1 = atoms[d.res1 - 1]
        res2 = atoms[d.res2 - 1]
        contact_key = res1[0] + res2[0]

        if mapping_mode == "allatomwesthof":
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
    def __init__(self, res1, res2, use_contact=True, weight=1):
        self.res1 = res1
        self.res2 = res2
        self.useContact = use_contact
        self.weight = weight

    def __str__(self):
        return "[%s, %s], useContact=%s, weight=%f" % (self.res1, self.res2, self.useContact, self.weight)

    def get_rosetta_function(self, function="FADE -100 26 20 -2 2"):
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
                parse_function(function, [(1, int), (2, int), (3, int), (4, int), (5, int)])  # TODO: make this float?
                return [function[0], function[1], function[2], function[3], function[4] * self.weight, function[5] * self.weight]
            else:
                raise DcaException("Invalid Rosetta function: '%s' is not implemented! Only 'FADE' function is recognized at this time." % function[0])
        except IndexError:
            raise DcaException("Invalid Rosetta function: Not enough arguments")


# DCA FILTERING

# run dca data through a chain of filters
def filter_dca_data(dca_data, dca_filter_chain, quiet=False):
    if dca_filter_chain is None:
        return dca_data
    for dcaFilter in dca_filter_chain:
        if dcaFilter is not None:
            for d in dca_data:
                dcaFilter(d, quiet)


# parses a text string and turns it into a filter chain
# multiple filters are separated b "," and their fields are separated using ":"
# example: threshold:8.0:100rnaDCA_FADE_-100_26_20_-2_2:cluster:1,threshold:-6.0:100rnaDCA_FADE_-100_26_20_-2_2:cluster:1
# TODO: document filter format
def parse_filter_line(line, simulation):
    filter_chain = []

    # filters are split by ","
    for filtertext in line.split(","):
        # initialize to None to suppress warning
        dca_filter = None
        # filter fields are split by ":"
        fields = filtertext.split(":")
        if fields[0] == "none":
            continue
        elif fields[0] == "threshold":
            threshold = float(fields[1])
            cst_name = fields[2]
            model_kind = fields[3]
            model_i = fields[4]
            model = simulation.get_models(cst_name, [model_i], model_kind)[0]
            pdbfile = model["pdbfile"]
            if not os.path.isfile(pdbfile):
                simulation.extract_pdb(cst_name, model)
            filter_pdb_chain = pdbtools.parse_pdb("", pdbfile)[0].child_list[0]
            if threshold > 0:
                dca_filter = dca_filter_threshold_minimum_keep_below(threshold, filter_pdb_chain)
            else:
                dca_filter = dca_filter_threshold_minimum_keep_above(abs(threshold), filter_pdb_chain)
        filter_chain.append(dca_filter)

    return filter_chain


# use contact if realized minimum distance in a single pdb is smaller than a threshold
def dca_filter_threshold_minimum_keep_below(threshold, pdb_chain):
    return _dca_filter_threshold_minimum_keep(threshold, pdb_chain, below=True)


# use contact if realized minimum distance is a single pdb is larger than a threshold
def dca_filter_threshold_minimum_keep_above(threshold, pdb_chain):
    return _dca_filter_threshold_minimum_keep(threshold, pdb_chain, below=False)


def _dca_filter_threshold_minimum_keep(threshold, pdb_chain, below=True):
    def f(contact, quiet):
        # do not touch contacts that are already disabled
        if not contact.useContact:
            return

        # get contact information
        average_heavy, minimum_heavy, minimum_pair = get_contact_information_in_pdb_chain(contact, pdb_chain)
        if not quiet:
            print "Threshold filter: min_distance: %f, threshold: %f, keep %s" % (minimum_heavy, threshold, "below" if below else "above")

        # set contact to disabled if filter failed
        if (below and minimum_heavy >= threshold) or (not below and minimum_heavy <= threshold):
            contact.useContact = False

    # return filter function
    return f
