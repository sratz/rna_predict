# coding: utf-8

import glob
import numpy as np
import os
from os.path import splitext, basename, abspath
import pickle
import pprint
import re
import shutil
import string
import struct
import subprocess
import sys
import time

from . import dcatools
from . import pdbtools
from . import utils


def check_file_existence(path, alternative_error_text=None):
    """
    Make sure a file exists and raise an exception otherwise.

    :param path: filename to check
    :param alternative_error_text: alternative error text passed to exception
    """
    if not os.path.isfile(path):
        raise SimulationException(alternative_error_text if alternative_error_text is not None else "Cannot find file: %s" % path)


def check_dir_existence(path, alternative_error_text=None):
    """
    Make sure a directory exists and raise an exception otherwise.

    :param path: directory to check
    :param alternative_error_text: alternative error text passed to exception
    """
    if not os.path.isdir(path):
        raise SimulationException(alternative_error_text if alternative_error_text is not None else "Cannot find directory: %s" % path)


def delete_glob(pattern, print_notice=True):
    """
    Helper function to delete files while expanding shell globs.

    :param pattern: pattern to delete
    :param print_notice: if True print verbose notice
    """
    for f in glob.glob(pattern):
        if print_notice:
            print "deleting %s..." % f
        try:
            if os.path.isdir(f):
                shutil.rmtree(f, ignore_errors=True)
            else:
                os.remove(f)
        except IOError:
            pass


def merge_silent_files(target, sources):
    """
    Merges rosetta silent files (.out).

    All source files are appended to target.
    Model numbers are incremented uniquely.
    Source files are deleted after a successful merge.

    :param target: target filename
    :param sources: source filenames
    """
    n = 0
    pattern_header = re.compile("^(?:SEQUENCE:|SCORE:\s+score).*")
    pattern_normal = re.compile("^(.*)S_\d+$")
    # if target already exists, see how many structures we have already
    if os.path.isfile(target):
        with open(target, "r") as t:
            line = None
            for line in t:
                pass
            if line is not None and 'line' in locals():
                m = re.match(".*S_(\d+)$",  line)
                if m:
                    n = int(m.group(1))

    # loop through all source files and append all new structures to target
    # while adjusting the numbering
    for source in glob.glob(sources):
        print "merging %s into %s" % (source, target)
        with open(target, "a+") as t:
            with open(source, "r") as s:
                for line in s:
                    m = pattern_header.match(line)
                    if m:
                        # only write header if this is he first one
                        if n == 0:
                            t.write(line)
                        continue
                    m = pattern_normal.match(line)
                    if m:
                        if line.startswith("SCORE:"):
                            n += 1
                        t.write("%sS_%06d\n" % (m.group(1), n))
                        continue
                    t.write(line)

    # delete source files
    delete_glob(sources, print_notice=True)

    # return total number of structures
    return n


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    """Helper function to be used as sort key in sorted() to naturally sort numeric parts."""
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]


# TODO: Correcting atom names IS MOST LIKELY WRONG! It is commented out for now and can probably be removed completely.
def fix_atom_names_in_cst(cst_info, sequence):
    """
    Switches N1 and N3 atom names in a residue.
    :param cst_info: list of constraints
    :param sequence: residue sequence in lower case
    :return: modified list of constraints
    """
    pyrimidines = ['c', 'u']
    purines = ['a', 'g']
    cst_info_new = []
    for cst in cst_info:
        atom_name1, res1, atom_name2, res2, function = cst
        if sequence[res1 - 1] in pyrimidines and atom_name1 == 'N1':
            atom_name1 = 'N3'
            print 'correcting atom name for ', res1
        if sequence[res2 - 1] in pyrimidines and atom_name2 == 'N1':
            atom_name2 = 'N3'
            print 'correcting atom name for ', res2
        if sequence[res1 - 1] in purines and atom_name1 == 'N3':
            atom_name1 = 'N1'
            print 'correcting atom name for ', res1
        if sequence[res2 - 1] in purines and atom_name2 == 'N3':
            atom_name2 = 'N1'
            print 'correcting atom name for ', res2
        cst_info_new.append([atom_name1, res1, atom_name2, res2, function])
    return cst_info_new


class SimulationException(Exception):
    """
    Custom exception class for foreseeable prediction errors.
    """
    pass


class Command(object):
    """
    Helper class to store external program calls plus additional flags.
    """

    def __init__(self, command, add_suffix=None, dry_run=False, print_commands=True, stdin=None, quiet=False):
        """
        Create a new command wrapper.

        :param command: list of external command and parameters
        :param add_suffix: type of suffix to first entry in command. currently only None or "rosetta" are supported
        :param dry_run: don't actually execute anything
        :param print_commands: print commandline when running
        :param stdin: optional text to use as standard input
        :param quiet: hide output of external program
        """
        self.command = command
        self.add_suffix = add_suffix
        self.dry_run = dry_run
        self.print_commands = print_commands
        self.stdin = stdin
        self.quiet = quiet

    def get_full_command(self, sysconfig):
        """
        If necessary, add the suffix to the command. Depending on user / system configuration.

        :param sysconfig: instance of SysConfig
        """
        com = self.command[:]
        if self.add_suffix == "rosetta":
            com[0] = sysconfig.rosetta_exe_path + com[0] + sysconfig.rosetta_exe_suffix
        if sysconfig.subprocess_buffsize is not None:
            com = ["stdbuf", "-o", sysconfig.subprocess_buffsize] + com
        return com


class EvalData(object):
    """
    Evaluation data storing model and cluster information.
    """
    def __init__(self):
        self.models = {}
        self.clusters = {}
        self.cluster_cutoff = 0.0

    @staticmethod
    def load_from_cst(constraints):
        """
        Load evaluation data from a prediction.

        :param constraints: name of a prediction
        """
        cst_name, cst_file = RNAPrediction.parse_cst_name_and_filename(constraints)
        return EvalData.load_from_file("predictions/%s/output/evaldata.dat" % cst_name)

    @staticmethod
    def load_from_file(filename):
        """
        Load evaluation data from a file.

        :param filename: path to a file containing the evaluation data
        """
        try:
            with open(filename, "r") as f:
                eval_data = pickle.load(f)
                if "rmsd_cluster_1" not in next(eval_data.models.itervalues()):
                    raise SimulationException("evaluation data file '%s' does not contain needed information" % filename)
                return eval_data
        except (IOError, pickle.PickleError, AttributeError, EOFError, IndexError, KeyError):
            # file broken or missing, force full evaluation
            raise SimulationException("evaluation data file '%s' missing or broken" % filename)

    def save_to_file(self, filename):
        """
        Dump evaldata to a file.

        :param filename: path to a file to store the data
        """
        with open(filename, "w") as f:
            pickle.dump(self, f)

    def get_model_count(self):
        """
        Returns the number of models in the vaulation data

        :return: number of models
        """
        return len(self.models)

    def get_models(self, model_list, kind="tag"):
        """Returns a selection of models with specific properties

        :param model_list: list of models to get. might a list of numbers, or a list of model names
        :param kind: model selection mode: "tag": string: internal name such as S_000123_5
                                           "top": number: models ordered by score
                                           "ntop": number: models ordered by rmsd_native
                                           "cluster": number: nth cluster decoy
                                           "cluster_ntop": n[/m]: nth cluster ordered by native rmsd [of first m]
        :return: list of selected models
        """
        # initialize list
        models_sorted = None

        results = []

        for model_i in model_list:
            # decide which model we want
            if kind == "cluster":
                model_i = int(model_i)
                tag = self.clusters[model_i]["primary_model"]
            elif kind == "cluster_ntop":
                if "/" in model_i:
                    i, limit = (int(x) for x in model_i.split("/"))
                    if i > limit:
                        raise SimulationException("get_models: Cluster number (%d) can't be larger than the limit (%d)" % (i, limit))
                else:
                    i = int(model_i)
                    limit = None
                model_selection = [self.models[m["primary_model"]] for m in self.clusters.values()[:limit]]
                models_sorted = sorted(model_selection, key=lambda i: i["rmsd_native"])
                tag = models_sorted[i - 1]["tag"]
            elif kind == "tag":
                tag = model_i
            elif kind == "top":
                model_i = int(model_i)
                # do we need to create a sorted list?
                if models_sorted is None:
                    models_sorted = sorted(self.models.items(), key=lambda x: x[1]["score"])
                tag = models_sorted[model_i - 1][0]
            elif kind == "ntop":
                model_i = int(model_i)
                if models_sorted is None:
                    try:
                        models_sorted = sorted(self.models.items(), key=lambda x: x[1]["rmsd_native"])
                    except KeyError:
                        raise SimulationException("get_models: No native rmsd information stored")
                tag = models_sorted[model_i - 1][0]
            else:
                raise SimulationException("get_models: Invalid 'kind' parameter: '%s'" % kind)

            if tag not in self.models:
                raise SimulationException("get_models: Invalid tag: '%s'" % tag)
            results.append(self.models[tag])

        return results

    def print_models(self, model_list, kind="tag"):
        """
        Print a list of models

        :param model_list: see get_models
        :param kind: see get_models
        """
        models = self.get_models(model_list, kind)
        for model in models:
            print "Model: %s" % model["tag"]
            pprint.pprint(model)


class RNAPrediction(object):
    """
    Base class used for prediction simulation
    """

    CONFIG_FILE = ".config"

    def load_config(self):
        """Load simulation configuration of current directory"""
        if os.path.exists(RNAPrediction.CONFIG_FILE):
            c = open(RNAPrediction.CONFIG_FILE)
            self.config = pickle.load(c)
            c.close()

    def save_config(self):
        """Save simulation configuration of current directorx"""
        c = open(RNAPrediction.CONFIG_FILE, "w")
        pickle.dump(self.config, c)
        c.close()

    def print_config(self):
        """Print simulation configuration."""
        print "Simulation configuration:"
        print "    path: %s" % (abspath(os.getcwd()))
        if self.config:
            for key, value in sorted(self.config.iteritems()):
                if key == "motif_res_maps" or key == "motif_stem_sets" or key == "motifs" or key == "stems" or key == "evaluate":
                    print "    %s: ..." % key
                else:
                    print "    %s: %s" % (key, value)
        else:
            print "    No configuration found."

    # TODO: support different data types?
    def modify_config(self, key, value):
        """
        Modify configuration entry.

        :param key: setting to modify
        :param value: new value ("-" to store None)
        """
        if key in self.config:
            self.config[key] = None if value == "-" else value
        else:
            raise SimulationException("No such config entry: %s" % key)

    def check_config(self):
        """Check if current directory contains a valid configuraion."""
        if not self.config:
            raise SimulationException("No config file found. Please run 'prepare' first!")

    def __init__(self, sysconfig):
        """Create a prediction simulation for the current directory and try to load an existing configuration."""
        self.config = {}
        self.sysconfig = sysconfig
        self.load_config()

    def execute_command(self, command, add_suffix=None, dry_run=False, print_commands=True, stdin=None, quiet=False):
        """
        Execute an external command.

        :param command: list of external command and parameters
        :param add_suffix: type of suffix to first entry in command. currently only None or "rosetta" are supported
        :param dry_run: don't actually execute anything
        :param print_commands: print commandline when running
        :param stdin: optional text to use as standard input
        :param quiet: hide output of external program
        """
        self.execute_commands([Command(command, add_suffix, dry_run, print_commands, stdin, quiet)])

    def execute_command_and_capture(self, command, add_suffix=None, dry_run=False, print_commands=True, stdin=None, quiet=False):
        """
        Execute an external command while capturing output.

        :param command: list of external command and parameters
        :param add_suffix: type of suffix to first entry in command. currently only None or "rosetta" are supported
        :param dry_run: don't actually execute anything
        :param print_commands: print commandline when running
        :param stdin: optional text to use as standard input
        :param quiet: hide output of external program
        :return: generator over lines of standard output
        """
        c = Command(command, add_suffix, dry_run, print_commands, stdin, quiet)
        command = c.get_full_command(self.sysconfig)
        if print_commands:
            print " ".join(command)
        try:
            p = subprocess.Popen(command, stdin=(subprocess.PIPE if c.stdin is not None else None), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            if c.stdin is not None:
                p.stdin.write(c.stdin)
        except OSError, e:
            raise SimulationException("Failed to execute command: %s, Reason: %s" % (" ".join(c.command), e))
        for line in p.stdout:
            if not quiet:
                sys.stdout.write(line)
            yield line
        p.wait()
        if p.returncode != 0:
            raise SimulationException("Non-zero return code from executed command: %s" % " ".join(c.command))

    def execute_commands(self, commands, threads=1):
        """
        Execute a list of commands parallelly.

        :param commands: list of Commands
        :param threads: number of parallel invocations
        """
        processes = []
        stdout = None
        try:
            while True:
                while commands and len(processes) < threads:
                    c = commands.pop(0)
                    command = c.get_full_command(self.sysconfig)
                    if c.print_commands:
                        print " ".join(command)
                    if not c.dry_run:
                        stdout = None
                        if c.quiet and stdout is None:
                            stdout = open(os.devnull, "w")
                        stderr = stdout
                        try:
                            p = subprocess.Popen(command, stdin=(subprocess.PIPE if c.stdin is not None else None), stdout=stdout, stderr=stderr)
                            # store original command in process
                            p.command = command
                            processes.append(p)
                            if c.stdin is not None:
                                p.stdin.write(c.stdin)
                        except OSError, e:
                            raise SimulationException("Failed to execute command: %s, Reason: %s" % (" ".join(c.command), e))

                for p in processes:
                    if p.poll() is not None:
                        processes.remove(p)
                        if p.returncode != 0:
                            raise SimulationException("Non-zero return code from executed command: %s" % " ".join(p.command))

                if not processes and not commands:
                    break
                else:
                    time.sleep(0.1)

        except:
            # kill all other processes
            for p in processes:
                p.kill()
            raise
        finally:
            if stdout is not None:
                stdout.close()

    @staticmethod
    def _make_tag_with_dashes(int_vector):
        """
        Transform a list of values into a list of strings while using "n-m" for adjacent values.

        Example: [1, 2, 3, 5, 7, 8] -> ["1-3", "5", "7-8"]

        :param int_vector: list of values
        :return: list of strings
        """
        tag = []

        start_res = int_vector[0]
        for i in range(1, len(int_vector) + 1):
            if i == len(int_vector) or int_vector[i] != int_vector[i - 1] + 1:

                stop_res = int_vector[i - 1]
                if stop_res > start_res:
                    tag += ["%d-%d" % (start_res, stop_res)]
                else:
                    tag += ["%d" % stop_res]

                if i < len(int_vector):
                    start_res = int_vector[i]

        return tag

    def prepare(self, fasta_file="sequence.fasta", params_file="secstruct.txt", native_pdb_file=None, data_file=None, torsions_file=None, name=None):
        """
        Preparation step. Parse out stems and motifs from sequence and secondary structure information and create necessary base files.

        :param fasta_file: fasta file containing the sequence
        :param params_file: text file containing the secondary structure
        :param native_pdb_file: native pdb file if available
        :param data_file: additional data file if available
        :param torsions_file: additional torsions file if available
        :param name: optional name for this set of predictions
        """
        if name is None:
            name = os.path.basename(os.path.abspath(os.getcwd()))
        for f in [fasta_file, params_file, native_pdb_file, data_file, torsions_file]:
            if f is not None:
                check_file_existence(f)

        utils.mkdir_p("constraints")
        utils.mkdir_p("preparation")
        utils.mkdir_p("dca")
        utils.mkdir_p("predictions")
        self.config["fasta_file"] = fasta_file
        self.config["params_file"] = params_file
        self.config["native_pdb_file"] = native_pdb_file
        self.config["data_file"] = data_file
        self.config["torsions_file"] = torsions_file
        self.config["name"] = name

        print
        print "Preparation:"

        # Read in files
        for line in utils.read_file_line_by_line(fasta_file):
            if line[0] == ">":
                continue
            self.config["sequence"] = line.lower()
            break
        print self.config["sequence"]
        numres = len(self.config["sequence"])

        # Read in data information
        self.config["data_info"] = []
        if data_file is not None:
            self.config["backbone_burial_info"] = []
            for line in utils.read_file_line_by_line(data_file):
                if len(line) > 6 and line[:6] == 'EXPOSE':
                    cols = string.split(line)
                    for i in range(len(cols) / 3):
                        self.config["data_info"].append([int(cols[3 * i + 1]), cols[3 * i + 2], cols[3 * i + 3]])
                if len(line) > 15 and line[:15] == 'BACKBONE_BURIAL':
                    cols = string.split(line)
                    for i in range(1, len(cols)):
                        self.config["backbone_burial_info"].append(int(cols[i]))

        pair_map = {}
        all_pairs = []

        complement = {'a': ['u'], 'u': ['a', 'g'], 'c': ['g'], 'g': ['c', 'u']}

        # Parse out stems
        basepair_mismatch = []
        cutpoints_original = []
        for line in utils.read_file_line_by_line(params_file):
            if line[:4] == 'STEM':
                cols = string.split(line)
                for i in range(len(cols)):
                    if cols[i] == 'PAIR':
                        # Offset to get to python numbering (starts with zero)
                        res1 = int(cols[i + 1]) - 1
                        res2 = int(cols[i + 2]) - 1
                        pair_map[res1] = res2
                        pair_map[res2] = res1
                        all_pairs.append([res1, res2])
                        # check if the secondary structure watson-crick pair is allowed
                        if not self.config["sequence"][res1] in complement[self.config["sequence"][res2]]:
                            basepair_mismatch += [[res1, res2]]
            elif line.count('(') > 0:  # maybe dot/bracket notation (((...)))
                print line
                self.config["secstruc"] = line
                left_brackets = []
                for i in range(len(line)):
                    if line[i] == '(':
                        left_brackets.append(i)
                    if line[i] == ')':
                        res1 = left_brackets[-1]
                        res2 = i
                        del (left_brackets[-1])
                        pair_map[res1] = res2
                        pair_map[res2] = res1
                        all_pairs.append([res1, res2])
                        # check if the secondary structure watson-crick pair is allowed
                        if not self.config["sequence"][res1] in complement[self.config["sequence"][res2]]:
                            basepair_mismatch += [[res1, res2]]
                if left_brackets:
                    raise SimulationException("unbalanced paranthesis in secondary structure")
            else:
                cols = string.split(line)
                res1 = int(cols[0]) - 1
                res2 = int(cols[1]) - 1
                pair_map[res1] = res2
                pair_map[res2] = res1
                all_pairs.append([res1, res2])
                # check if the secondary structure watson-crick pair is allowed
                if not self.config["sequence"][res1] in complement[self.config["sequence"][res2]]:
                    basepair_mismatch += [[res1, res2]]

            if line[:13] == 'CUTPOINT_OPEN':
                cutpoints_original = map(lambda x: int(x), string.split(line[14:]))

        # print all helix pairs that don't fit the allowed complements
        if basepair_mismatch:
            basepair_mismatch_str = [" "] * len(self.config["sequence"])
            for r in basepair_mismatch:
                basepair_mismatch_str[r[0]] = "^"
                basepair_mismatch_str[r[1]] = "^"
            print "".join(basepair_mismatch_str)

            # more verbose:
            for r in basepair_mismatch:
                res1 = r[0]
                res2 = r[1]
                resstr1 = str(res1 + 1) + self.config["sequence"][res1]
                resstr2 = str(res2 + 1) + self.config["sequence"][res2]
                print " " * res1 + "^" + " " * (res2 - res1 - 1) + "^"
                print " " * res1 + resstr1 + " " * (res2 - res1 - len(resstr1)) + resstr2

            raise SimulationException("some helix pairs in secondary structure cannot be realized (see above)")

        # print pair_map

        # Parse out stems
        already_in_stem = {}
        for i in range(numres):
            already_in_stem[i] = 0

        self.config["stems"] = []
        for i in range(numres):
            if i in pair_map and not already_in_stem[i]:  # In a base pair
                k = i
                stem_res = [[k, pair_map[k]]]

                already_in_stem[k] = 1
                already_in_stem[pair_map[k]] = 1

                # Can we extend in one direction?
                while (k + 1) in pair_map and pair_map[k + 1] == pair_map[k] - 1 and not already_in_stem[k + 1]:
                    k += 1
                    stem_res.append([k, pair_map[k]])
                    already_in_stem[k] = 1
                    already_in_stem[pair_map[k]] = 1

                # Do not allow single WC base pairs.
                if len(stem_res) < 2:
                    print 'All stems must have length > 1 bp '
                    print stem_res
                    exit()
                self.config["stems"].append(stem_res)

        # Parse out motifs
        already_in_motif = {}
        for i in range(numres):
            already_in_motif[i] = 0

        motif_cutpoints = []
        self.config["motif_res_maps"] = []
        self.config["motif_stem_sets"] = []
        self.config["motifs"] = []
        for i in range(numres):

            if not already_in_motif[i] and (not already_in_stem[i] or (i > 0
                                                                       and already_in_stem[i - 1]
                                                                       and pair_map[i] + 1 != pair_map[i - 1])):

                motif_res = []
                motif_stem_set = []
                cutpoints = []

                if i > 1:

                    # back up to beginning of stem.
                    k = i - 1

                    # first base pair.
                    motif_stem = [[k, pair_map[k]]]
                    motif_res.append(k)
                    motif_res.append(pair_map[k])

                    k -= 1
                    while k >= 0 and already_in_stem[k] and \
                            (pair_map[k] - 1 == pair_map[k + 1]):
                        motif_stem.append([k, pair_map[k]])
                        motif_res.append(k)
                        motif_res.append(pair_map[k])
                        k -= 1
                    motif_stem_set.append(motif_stem)
                    k += 1
                    cutpoints.append(pair_map[k])

                # print 'AFTER FIRST HELIX: ', motif_res

                k = i
                while k not in motif_res and k < numres:
                    # Move forward to next stem:
                    while k < numres and not already_in_stem[k]:
                        if already_in_motif[k]:
                            print 'Hey cant deal with pseudoknots!'
                            exit()
                        motif_res.append(k)
                        k += 1

                    stem_start = k

                    if k >= numres:
                        cutpoints.append(k - 1)
                        break

                    if k in motif_res:
                        break

                    motif_stem = [[k, pair_map[k]]]
                    motif_res.append(k)
                    motif_res.append(pair_map[k])
                    k += 1

                    while k < numres and already_in_stem[k] and pair_map[k - 1] == pair_map[k] + 1 and k not in motif_res:
                        motif_stem.append([k, pair_map[k]])
                        motif_res.append(k)
                        motif_res.append(pair_map[k])
                        k += 1
                    motif_stem_set.append(motif_stem)
                    cutpoints.append(k - 1)

                    # Next non-helical part..
                    k = pair_map[stem_start] + 1

                    # print 'AFTER NEXT HELIX: ', motif_res

                motif_res.sort()

                motif_res_map = {}
                for k in range(len(motif_res)):
                    motif_res_map[motif_res[k]] = k
                    already_in_motif[motif_res[k]] = 1

                self.config["motifs"].append(motif_res)
                self.config["motif_stem_sets"].append(motif_stem_set)
                motif_cutpoints.append(cutpoints)
                self.config["motif_res_maps"].append(motif_res_map)
                # print 'CUTPOINTS ', cutpoints

        for i in range(len(self.config["stems"])):

            # Fasta
            tag = 'preparation/stem%d.fasta' % (i + 1)
            fid = open(tag, 'w')
            fid.write('>' + tag + '\n')

            stem_res = self.config["stems"][i]
            stem_length = len(stem_res)

            for k in range(stem_length):
                fid.write(self.config["sequence"][stem_res[k][0]])
                # print stem_res[k][0]+1,
            for k in range(stem_length):
                fid.write(self.config["sequence"][stem_res[stem_length - k - 1][1]])
                # print stem_res[stem_length-k-1][1]+1,

            # print
            fid.write('\n')
            fid.close()
            print 'Created: ', tag

            # pdb_file
            if native_pdb_file is not None:
                command = ["pdbslice.py",
                           native_pdb_file,
                           "-segment",
                           "%d" % (stem_res[0][0] + 1),
                           "%d" % (stem_res[-1][0] + 1),
                           "%d" % (stem_res[-1][-1] + 1),
                           "%d" % (stem_res[0][-1] + 1),
                           "preparation/stem%d_" % (i + 1)]
                self.execute_command(command)
                native_pdb_file_subset = 'preparation/stem%d_%s' % (i + 1, native_pdb_file)
                print 'Created: ', native_pdb_file_subset

        # Output motif jobs
        for i in range(len(self.config["motifs"])):

            # Fasta
            motif_fasta_file = 'preparation/motif%d.fasta' % (i + 1)
            fid = open(motif_fasta_file, 'w')
            fid.write('>' + motif_fasta_file + '\n')

            motif_res = self.config["motifs"][i]
            motif_length = len(motif_res)

            for k in range(motif_length):
                fid.write(self.config["sequence"][motif_res[k]])
            fid.write('\n')
            fid.close()
            print 'Created: ', motif_fasta_file

            # params file
            motif_stem_set = self.config["motif_stem_sets"][i]
            motif_res_map = self.config["motif_res_maps"][i]
            motif_cutpoint = motif_cutpoints[i]

            motif_params_file = 'preparation/motif%d.params' % (i + 1)
            fid = open(motif_params_file, 'w')

            for k in range(len(motif_stem_set)):
                motif_stem = motif_stem_set[k]
                fid.write('STEM   ')
                for n in range(len(motif_stem)):
                    fid.write('  PAIR %d %d W W A' %
                              (motif_res_map[motif_stem[n][0]] + 1,
                               motif_res_map[motif_stem[n][1]] + 1))
                fid.write('\n')

                fid.write('OBLIGATE PAIR %d %d W W A \n\n' %
                          (motif_res_map[motif_stem[-1][0]] + 1,
                           motif_res_map[motif_stem[-1][1]] + 1))

            motif_cutpoint.sort()

            # print motif_res
            # print motif_cutpoint

            if len(motif_cutpoint) > 1:
                fid.write('CUTPOINT_OPEN ')
                for k in range(len(motif_cutpoint)):
                    if motif_res_map[motif_cutpoint[k]] < (len(motif_res) - 1):
                        fid.write(' %d' % (motif_res_map[motif_cutpoint[k]] + 1))
            fid.write('\n')
            fid.close()
            print 'Created: ', motif_params_file

            # pdb_file
            if native_pdb_file is not None:
                command = ["pdbslice.py",
                           native_pdb_file,
                           "-subset"]
                for k in range(motif_length):
                    command += ["%d" % (motif_res[k] + 1)]
                command += ["preparation/motif%d_" % (i + 1)]
                self.execute_command(command)
                native_pdb_file_subset = 'preparation/motif%d_%s' % (i + 1, native_pdb_file)
                print 'Created: ', native_pdb_file_subset

            if data_file is not None:
                motif_data_file = 'preparation/motif%d.data' % (i + 1)
                fid_data = open(motif_data_file, 'w')
                fid_data.write('EXPOSE')
                for data in self.config["data_info"]:
                    if data[0] - 1 in motif_res_map.keys():
                        fid_data.write('   %d %s %s ' % (motif_res_map[data[0] - 1] + 1, data[1], data[2]))
                fid_data.write('\n')

                if len(self.config["backbone_burial_info"]) > 0:
                    fid_data.write('BACKBONE_BURIAL ')
                    for k in self.config["backbone_burial_info"]:
                        if k - 1 in motif_res_map.keys():
                            fid_data.write(' %d' % (motif_res_map[k - 1] + 1))
                    fid_data.write('\n')
                fid_data.close()
                print 'Created: ', motif_data_file

        if len(cutpoints_original) > 0:
            # noinspection PyUnboundLocalVariable
            fid.write('CUTPOINT_OPEN ')
            for cutpoint in cutpoints_original:
                fid.write(' %d' % (cutpoint + 1))
            fid.write('\n')

        cutpoints = []
        for i in range(len(self.config["motifs"])):

            motif_stem_set = self.config["motif_stem_sets"][i]

            motif_stem = motif_stem_set[0]

            possible_cutpoints = [motif_stem[0][0], motif_stem[1][1]]
            possible_cutpoints.sort()
            # print possible_cutpoints
            if possible_cutpoints[0] not in cutpoints:
                cutpoints.append(possible_cutpoints[0])

        params_file = "preparation/sequence.params"
        fid = open(params_file, 'w')

        if len(cutpoints) > 0:
            fid.write('CUTPOINT_CLOSED ')
            # cutpoints.sort()
            for cutpoint in cutpoints:
                fid.write(' %d' % (cutpoint + 1))
            fid.write('\n')

        # for cutpoint in cutpoints:
        # fid.write( 'OBLIGATE   PAIR %d %d W W A\n' % (cutpoint+1, pair_map[cutpoint]+1) )

        for i in range(len(self.config["stems"])):
            stem_res = self.config["stems"][i]
            fid.write('STEM ')
            for k in range(len(stem_res)):
                fid.write(' PAIR %d %d W W A ' %
                          (stem_res[k][0] + 1, stem_res[k][1] + 1))
            fid.write('\n')
            fid.write('OBLIGATE PAIR %d %d W W A \n\n' %
                      (stem_res[-1][0] + 1,
                       stem_res[-1][1] + 1))

        fid.close()
        print 'Created: ', params_file

        # default cutpoint constraints
        assemble_cst_file = "preparation/cutpoints.cst"
        fid = open(assemble_cst_file, 'w')
        fid.write('[ atompairs ]\n')
        for cutpoint in cutpoints:
            fid.write("O3'  %d  P     %d  HARMONIC  1.619  2.0\n" % (cutpoint + 1, cutpoint + 2))
        fid.close()
        print 'Created: ', assemble_cst_file

    def create_helices(self, dry_run=False, threads=1):
        """
        Helix creation step. Create one ideal helix model for each helix.

        :param dry_run: don't actually run any external command
        :param threads: number of threads to use
        """
        self.check_config()
        commands = list()
        delete_glob("preparation/stem*.out")
        for i in range(len(self.config["stems"])):
            command = ["rna_helix",
                       "-fasta", "preparation/stem%d.fasta" % (i + 1),
                       "-out:file:silent", "preparation/stem%d.out" % (i + 1)]
            commands.append(Command(command, add_suffix="rosetta", dry_run=dry_run))
        self.execute_commands(commands, threads=threads)

        # rna_helix dumps <sequence>.pdb files in the working directory.
        # These can be useful for us, so move them to the preparation directory with a proper stemX.pdb filename.
        # TODO: What happens if there are multiple helices with the same sequence? In that case we need to move them
        # out of the way before running the next command.
        # Just disable multithreading here (not really needed with stem generation anyways)?
        for i in range(len(self.config["stems"])):
            sequence = [self.config["sequence"][j[0]] for j in self.config["stems"][i]]
            sequence += [self.config["sequence"][j[1]] for j in reversed(self.config["stems"][i])]
            shutil.move("%s.pdb" % ("".join(sequence)), "preparation/stem%d.pdb" % (i + 1))

    @staticmethod
    def parse_cst_name_and_filename(constraints):
        """
        Find and clean up constraints by name or filename

        :param constraints: constraints name or filename
        :return: tuple of constraints name and filename
        """
        if constraints is None or constraints == "none":
            cst_name = "none"
            cst_file = None
        elif os.path.isfile(constraints):
            cst_file = constraints
            cst_name = splitext(basename(cst_file))[0]
        else:
            constraints = basename(constraints)
            cst_name, cst_ext = splitext(constraints)
            if cst_ext != ".cst":
                cst_name = constraints
            cst_file = "constraints/%s.cst" % cst_name
            check_file_existence(cst_file)
        return cst_name, cst_file

    @staticmethod
    def parse_cst_file(constraints_file):
        """
        Parse .cst file as a list of constraints

        :param constraints_file: path to a .cst file
        :return: list of constraints
        """
        cst_info = []
        for line in utils.read_file_line_by_line(constraints_file):
            if len(line) > 6 and line[0] != '[':
                cols = string.split(line)
                atom_name1 = cols[0]
                res1 = int(cols[1])
                atom_name2 = cols[2]
                res2 = int(cols[3])
                cst_info.append([atom_name1, res1, atom_name2, res2, cols[4:]])
        return cst_info

    def prepare_cst(self, constraints=None, motifs_override=None):
        """
        Constraints preparation step. Prepare constraints files for motif generation and assembly.

        :param constraints: constraints selection
        :param motifs_override: optional name of a different set of constraints to use as motifs.
        """
        self.check_config()
        cst_name, cst_file = self.parse_cst_name_and_filename(constraints)
        print "Constraints preparation:"
        print "    constraints: %s" % cst_name

        dir_prediction = "predictions/%s" % cst_name
        dir_assembly = "predictions/%s/assembly" % cst_name
        dir_motifs = "predictions/%s/motifs" % cst_name

        if motifs_override:
            motifcst_name, motifcst_file = self.parse_cst_name_and_filename(motifs_override)
            print "    motif constraints: %s" % motifcst_name
            if motifcst_name == cst_name:
                raise SimulationException("Motif override cst name can't be the same as the cst name!")

        utils.mkdir_p("predictions/%s" % cst_name)
        utils.mkdir_p(dir_assembly)

        if motifs_override:
            # TODO: just make all that stuff a cst config file
            # create MOTIF_OVERRIDE file and store the cst name of the motifs
            file_motif_override = "%s/MOTIF_OVERRIDE" % dir_prediction
            with open(file_motif_override, "w") as f:
                f.write(motifcst_name)

        # assembly constraints: start with default harmonic cutpoints
        assembly_cst = "%s/assembly.cst" % dir_assembly
        shutil.copy("preparation/cutpoints.cst", assembly_cst)
        print "Created: %s" % assembly_cst

        # if cst is "none" we are done now
        if cst_file is None:
            return

        # load constraints information
        cst_info = self.parse_cst_file(cst_file)

        # probably wrong, commented out for now!
        # cst_info = fix_atom_names_in_cst(cst_info, self.config["sequence"])

        # add tertiary constraints
        with open(assembly_cst, "a") as ft:
            for cst in cst_info:
                ft.write('%s %d %s %d %s \n' % (cst[0], cst[1], cst[2], cst[3], " ".join(map(str, cst[4]))))
        print "Added tertiary constraints: %s" % assembly_cst

        # extract relevant constraints for motif generation from global cst file if necessary
        if not motifs_override:
            utils.mkdir_p(dir_motifs)
            for i in range(len(self.config["motifs"])):
                motif_res_map = self.config["motif_res_maps"][i]
                motif_cst_file = '%s/motif%d.cst' % (dir_motifs, i + 1)
                fid_cst = open(motif_cst_file, 'w')
                fid_cst.write('[ atompairs ]\n')
                for cst in cst_info:
                    if cst[1] - 1 in motif_res_map.keys() and cst[3] - 1 in motif_res_map.keys():
                        fid_cst.write('%s %d %s %d %s\n' % (cst[0], motif_res_map[cst[1] - 1] + 1, cst[2], motif_res_map[cst[3] - 1] + 1, " ".join(map(str, cst[4]))))
                fid_cst.close()
                print 'Created: ', motif_cst_file

    def create_motifs(self, nstruct=50000, cycles=20000, dry_run=False, seed=None, use_native_information=False, threads=1, constraints=None):
        """
        Motif generation step. Generate models for each motif.

        :param nstruct: number of models to create for each motif
        :param cycles: number of monte-carlo cycles per model
        :param dry_run: don't actually run any external command
        :param seed: optionally override random seed
        :param use_native_information: use native pdb file to calculate rmsd for each model
        :param threads: number of threads to use
        :param constraints: constraints selection
        """
        self.check_config()
        cst_name, cst_file = self.parse_cst_name_and_filename(constraints)
        print "Motif creation configuration:"
        print "    constraints: %s" % cst_name
        print "    cycles: %s" % cycles
        print "    nstruct: %s" % nstruct
        print "    dry_run: %s" % dry_run
        print "    random_seed: %s" % seed
        print "    threads: %s" % threads

        dir_motifs = "predictions/%s/motifs" % cst_name
        utils.mkdir_p("predictions/%s" % cst_name)
        utils.mkdir_p(dir_motifs)

        n_motifs = len(self.config["motifs"])

        # check if motif constraints were created correctly
        if cst_file is not None:
            for i in range(n_motifs):
                check_file_existence("%s/motif%d.cst" % (dir_motifs, i + 1), "Motif cst files not found. Please run the 'prepare-cst' step!")

        # merge all motifs abd check what we have so far
        completed = {}
        for i in range(n_motifs):
            completed[i] = merge_silent_files("%s/motif%d.out" % (dir_motifs, i + 1), "%s/motif%d_*.out" % (dir_motifs, i + 1))

        commands = list()
        for i in range(n_motifs):
            # split the rest of the work in multiple jobs
            structs_missing = nstruct - completed[i]
            print "motif %d:  completed: %d  missing: %d" % (i + 1, completed[i], structs_missing)
            if structs_missing <= 0:
                continue

            structs_threads = [0 for _ in range(threads)]
            for j in range(threads):
                structs_threads[j] = structs_missing / threads
            for j in range(structs_missing % threads):
                structs_threads[j] += 1

            command = ["rna_denovo",
                       "-fasta", "preparation/motif%d.fasta" % (i + 1),
                       "-params_file", "preparation/motif%d.params" % (i + 1),
                       "-cycles", "%d" % cycles,
                       "-mute", "all",
                       "-close_loops",
                       "-close_loops_after_each_move",
                       "-minimize_rna"]

            if cst_file is not None:
                command += ["-cst_file", '%s/motif%d.cst' % (dir_motifs, i + 1)]
            if self.config["native_pdb_file"] is not None and use_native_information:
                command += ["-native", "preparation/motif%d_%s" % (i + 1, self.config["native_pdb_file"])]
            if self.config["data_file"] is not None:
                command += ["-data_file", "preparation/motif%d.data" % (i + 1)]
            if self.config["torsions_file"] is not None:
                command += ["-vall_torsions", "preparation/motif%d.torsions" % (i + 1)]

            which_stems = []
            stem_chunk_res = []
            motif_stem_set = self.config["motif_stem_sets"][i]
            motif_res_map = self.config["motif_res_maps"][i]
            for k in range(len(motif_stem_set)):
                motif_stem = motif_stem_set[k]
                # need to find in stems
                for n in range(len(self.config["stems"])):
                    stem = self.config["stems"][n]
                    found_match = 0
                    for q in range(len(stem)):
                        if motif_stem[0][0] in stem[q]:
                            found_match = 1
                            break
                    if found_match:
                        which_stems.append(n)
                        break

                for q in range(len(stem)):
                    stem_chunk_res.append(motif_res_map[stem[q][0]] + 1)
                for q in range(len(stem)):
                    stem_chunk_res.append(motif_res_map[stem[len(stem) - q - 1][1]] + 1)

            command += ["-in:file:silent_struct_type", "rna",
                        "-in:file:silent"]
            for n in which_stems:
                command += ["preparation/stem%d.out" % (n + 1)]

            command += ["-in:file:input_res"]
            command += self._make_tag_with_dashes(stem_chunk_res)

            for j in range(threads):
                if structs_threads[j] == 0:
                    continue
                command_full = command + ["-out:file:silent", "%s/motif%d_%d.out" % (dir_motifs, i + 1, j + 1),
                                          "-nstruct", "%d" % (structs_threads[j])]
                if seed is not None:
                    command_full += ["-constant_seed", "-jran", "%d" % seed]
                    seed += 1

                commands.append(Command(command_full, add_suffix="rosetta", dry_run=dry_run))
        self.execute_commands(commands, threads=threads)

        # merge motifs
        for i in range(n_motifs):
            merge_silent_files("%s/motif%d.out" % (dir_motifs, i + 1), "%s/motif%d_*.out" % (dir_motifs, i + 1))

    def assemble(self, nstruct=50000, cycles=20000, constraints=None, dry_run=False, seed=None, use_native_information=False, threads=1):
        """
        Assembly step. Assemble helices and motifs into complete models.

        :param nstruct: number of models to create
        :param cycles: number of monte-carlo cycles per model
        :param constraints: constraints selection
        :param dry_run: don't actually run any external command
        :param seed: optionally override random seed
        :param use_native_information: use native pdb file to calculate rmsd for each model
        :param threads: number of threads to use
        """
        self.check_config()
        cst_name, cst_file = self.parse_cst_name_and_filename(constraints)
        print "Assembly configuration:"
        print "    constraints: %s" % cst_name
        print "    cycles: %s" % cycles
        print "    nstruct: %s" % nstruct
        print "    dry_run: %s" % dry_run
        print "    random_seed: %s" % seed

        # Check if a different set of motifs should be used.
        file_motif_override = "predictions/%s/MOTIF_OVERRIDE" % cst_name
        if os.path.isfile(file_motif_override):
            with open(file_motif_override, "r") as f:
                motifcst_name = f.read().strip()
                print "    using motifs: %s" % motifcst_name
        else:
            motifcst_name = cst_name

        dir_assembly = "predictions/%s/assembly" % cst_name
        dir_motifs = "predictions/%s/motifs" % motifcst_name
        file_assembly_cst = "%s/assembly.cst" % dir_assembly

        check_file_existence(file_assembly_cst, "Preparation incomplete. Please run 'prepare-cst'.")
        check_file_existence("%s/motif1.out" % dir_motifs, "Motifs not found. Please run 'create-motifs'.")

        commands = list()

        command = ["rna_denovo",
                   "-minimize_rna",
                   "-fasta", self.config["fasta_file"],
                   "-in:file:silent_struct_type", "binary_rna",
                   "-cycles", "%d" % cycles,
                   "-params_file", "preparation/sequence.params",
                   "-close_loops",
                   "-in:file:silent"]

        for i in range(len(self.config["stems"])):
            command += ["preparation/stem%d.out" % (i + 1)]

        for i in range(len(self.config["motifs"])):
            command += ["%s/motif%d.out" % (dir_motifs, i + 1)]

        if cst_file is not None:
            command += ["-cst_file", file_assembly_cst]

        chunk_res = []
        for n in range(len(self.config["stems"])):
            stem = self.config["stems"][n]
            for q in range(len(stem)):
                chunk_res.append(stem[q][0] + 1)
            for q in range(len(stem)):
                chunk_res.append(stem[len(stem) - q - 1][1] + 1)

        for n in range(len(self.config["motifs"])):
            motif_res = self.config["motifs"][n]
            for m in motif_res:
                chunk_res.append(m + 1)

        command += ["-in:file:input_res"]
        command += self._make_tag_with_dashes(chunk_res)

        if self.config["native_pdb_file"] is not None and use_native_information:
            command += ["-native", self.config["native_pdb_file"]]
        if self.config["torsions_file"] is not None:
            command += ["-vall_torsions", self.config["torsions_file"]]
        if self.config["data_file"] is not None:
            command += ["-data_file", self.config["data_file"]]

        # distribute nstruct among threads
        structs_threads = [0 for _ in range(threads)]
        for j in range(threads):
            structs_threads[j] = nstruct / threads
        for j in range(nstruct % threads):
            structs_threads[j] += 1

        for j in range(threads):
            if structs_threads[j] == 0:
                continue

            # In case no seed is specified, we do what rosetta does to get a random number, and seed it with that.
            # This way we know the number beforehand and can choose an appropriate output filename.
            # In case a seed was specified, we increment it here with the thread number.
            if seed is None:
                seed_used = struct.unpack("=i", os.urandom(4))[0]
            else:
                seed_used = seed + j

            command_full = command + ["-out:file:silent", "%s/assembly_%d.out" % (dir_assembly, seed_used),
                                      "-constant_seed", "-jran", "%d" % seed_used, "--nstruct", "%d" % structs_threads[j]]

            commands.append(Command(command_full, add_suffix="rosetta", dry_run=dry_run))

        self.execute_commands(commands, threads=threads)

    # TODO: set the cutoff value back to 4.0? It was set to 4.1 because the old bash script used integer comparison and even 4.09 was treated as 4.0
    def evaluate(self, constraints=None, cluster_limit=10, cluster_cutoff=4.1, full_evaluation=False):
        """
        Evaluation step. Extract model information, cluster models and calculate rmsd values.

        :param constraints: constraints selection
        :param cluster_limit: maximum number of clusters to produce
        :param cluster_cutoff: rmsd distance in A at which a new cluster is created
        :param full_evaluation: discard existing evaluation data, re-extract model information and re-calculate rmsd values.
        :return:
        """
        self.check_config()
        cst_name, cst_file = self.parse_cst_name_and_filename(constraints)
        print "Evaluation configuration:"
        print "    constraints: %s" % cst_name
        print "    cluster_limit: %s" % cluster_limit
        print "    cluster_cutoff: %s" % cluster_cutoff

        dir_assembly = "predictions/%s/assembly" % cst_name
        dir_output = "predictions/%s/output" % cst_name
        dir_tmp = "predictions/%s/temp" % cst_name
        file_evaldata = "%s/evaldata.dat" % dir_output

        if not os.path.isdir("predictions/%s" % cst_name):
            raise SimulationException("No prediction directory for constraint '%s' found! Maybe you should assemble first?" % cst_name)

        files_assembly = sorted(glob.glob("%s/assembly_*.out" % dir_assembly), key=natural_sort_key)
        if len(files_assembly) == 0:
            raise SimulationException("No assembly files for constraint '%s' found! Maybe you should assemble first?" % cst_name)

        if not full_evaluation:
            if not os.path.isdir(dir_output) or not os.path.isfile(file_evaldata):
                # first time, silently do a full eval
                full_evaluation = True
            else:
                try:
                    eval_data = EvalData.load_from_file(file_evaldata)
                except SimulationException:
                    # file broken or missing, force full evaluation
                    print "Warning: %s missing or broken, running full evaluation." % file_evaldata
                    full_evaluation = True

        print "    full_evaluation: %s" % full_evaluation

        # cleanup
        if full_evaluation:
            delete_glob(dir_output)
            delete_glob(dir_tmp)
            # create empty evaluation data
            eval_data = EvalData()
        else:
            # clean out old cluster information from models
            for model in eval_data.models.itervalues():
                if "cluster" in model:
                    del model["cluster"]

        eval_data.clusters = {}

        delete_glob(dir_output + os.sep + "*.pdb")
        delete_glob(dir_output + os.sep + "*.out")
        utils.mkdir_p(dir_output)
        utils.mkdir_p(dir_tmp)

        # extract info from .out files
        if full_evaluation:
            # store all the rosetta scores in a dict
            scores_headers = None

            # loop over all out files matching the constraint
            for i, f in enumerate(files_assembly):
                print "processing rosetta silent file: %s..." % f

                # read out file and store a dict of the scores and the full score line
                regex_score = re.compile("^SCORE:\s+(.*)$")
                for line in utils.read_file_line_by_line(f):
                    m = regex_score.match(line)
                    if m:
                        # split rosetta score entries into a list
                        scores = m.group(1).split()
                        # store the header names if not already known
                        if scores[-1] == "description":
                            if not scores_headers:
                                scores_headers = scores
                            continue
                        # create a score dict
                        scores_dict = {}
                        for s_index, s_name in enumerate(scores_headers):
                            # store float values as float
                            try:
                                scores_dict[scores_headers[s_index]] = float(scores[s_index])
                            except ValueError:
                                scores_dict[scores_headers[s_index]] = scores[s_index]

                        # name models exactly like rosetta would (that is, append _<num> when already present)
                        j = 1
                        tag = scores_dict["description"]
                        while tag in eval_data.models:
                            tag = "%s_%d" % (scores_dict["description"], j)
                            j += 1
                        eval_data.models[tag] = {"source_file": basename(f),
                                                 "score": scores_dict["score"],
                                                 "tag": tag,
                                                 "tag_source": scores_dict["description"],
                                                 "rosetta_scores": scores_dict}

        # clustering
        eval_data.cluster_cutoff = cluster_cutoff
        print "clustering models..."
        filename_clusters = "%s/clusters.out" % dir_output
        sys.stdout.write("  ")

        regex_clusters = re.compile("^.*RNA_PREDICT (old|new) cluster ([0-9]+) score ([-0-9.]+) model (S_[0-9_]+)$")
        for line in self.execute_command_and_capture(["rna_cluster", "-in:file:silent"] + files_assembly + ["-cluster:radius", "%s" % cluster_cutoff, "-nstruct", str(cluster_limit), "-out:file:silent", filename_clusters], add_suffix="rosetta", quiet=True):
            m = regex_clusters.match(line)
            if m:
                eval_data.models[m.group(4)]["cluster"] = int(m.group(2))
                if m.group(1) == "new":
                    eval_data.clusters[int(m.group(2))] = {"primary_model": m.group(4)}
                print "  %s cluster: %02s, model: %-10s, score: %.3f" % ("new" if m.group(1) == "new" else "   ", m.group(2), m.group(4), float(m.group(3)))

        # extract cluster pdbs
        print "extracting cluster decoy pdbs..."
        owd = os.getcwd()
        try:
            os.chdir(dir_output)
            self.execute_command(["rna_extract", "-in:file:silent", basename(filename_clusters), "-in:file:silent_struct_type", "rna"], add_suffix="rosetta")
        finally:
            os.chdir(owd)

        # loop over all extracted pdb files
        for i, fe in enumerate(sorted(glob.glob("%s/S_*.pdb" % dir_output), key=natural_sort_key)):
            shutil.move(fe, "%s/cluster_%d.pdb" % (dir_output, i + 1))

        # rmsdcalc helper function
        # noinspection PyShadowingNames
        def calculate_rmsd(name, filename_comparison):
            # calculate native rmsd values for all models if native pdb available
            filename_tmp = "%s/rmsd.out" % dir_tmp
            print "caluculating rmsd values to %s for all models..." % name
            sys.stdout.write("  ")

            regex_rmsd = re.compile("^All atom rmsd over moving residues: (S_[0-9_]+) ([-0-9.]+)$")
            for line in self.execute_command_and_capture(["rna_score", "-in:file:silent"] + files_assembly + ["-in:file:native", filename_comparison, "-score:just_calc_rmsd", "-out:file:silent", filename_tmp], add_suffix="rosetta", quiet=True):
                m = regex_rmsd.match(line)
                if m:
                    print m.group(1), m.group(2)
                    eval_data.models[m.group(1)]["rmsd_%s" % name] = float(m.group(2))

            delete_glob(filename_tmp, print_notice=False)

        # only calculate rmsds if full_evalation is needed or wanted
        if full_evaluation:
            # calculate rmsd to native structure if native pdb available
            if self.config["native_pdb_file"] is not None:
                calculate_rmsd("native", self.config["native_pdb_file"])

            # calculate rmsd to best model
            calculate_rmsd("cluster_1", "%s/cluster_1.pdb" % dir_output)

        # save evaluation data
        eval_data.save_to_file(file_evaldata)

    def evaluate_custom(self, constraints=None, dca_prediction_filename="dca/dca.txt", full_evaluation=False, threshold=7.5, radius=2, number_dca_predictions=100, threads=4):
        """
        Custom scoring algorithmy by inspecting neighboring residues of dca contact pairs

        :param constraints: constraints selection
        :param dca_prediction_filename: dca filename
        :param full_evaluation: discard existing evaluation data and extracted models, re-extract and re-calculate distance information
        :param threshold: threshold in A to count a contact
        :param radius: number of adjacent residues to inspect
        :param number_dca_predictions: maximum number of DCA predictions to use
        :param threads: number of threads to use for extraction
        """
        self.check_config()
        cst_name, cst_file = self.parse_cst_name_and_filename(constraints)
        print "Evaluation configuration:"
        print "    constraints: %s" % cst_name
        print "    dca_file: %s" % dca_prediction_filename
        print "    dca_count: %s" % number_dca_predictions
        print "    threshold: %f" % threshold
        print "    radius: %s" % radius

        dir_assembly = "predictions/%s/assembly" % cst_name
        dir_output = "predictions/%s/output" % cst_name
        dir_tmp = "predictions/%s/temp" % cst_name
        file_evaldata = "%s/evaldata.dat" % dir_output
        file_matrices = "%s/matrices.dat" % dir_tmp

        files_assembly = sorted(glob.glob("%s/assembly_*.out" % dir_assembly), key=natural_sort_key)

        try:
            eval_data = EvalData.load_from_file(file_evaldata)
        except SimulationException:
            raise SimulationException("Normal --evaluate step needs to be run before --evaluate-custom.")

        slen = len(self.config["sequence"])

        utils.mkdir_p(dir_output)
        utils.mkdir_p(dir_tmp)

        extract = full_evaluation
        build_matrices = full_evaluation

        # check if we can use existing data
        if not build_matrices:
            # try to load existing distance matrix
            print "loading existing distance matrices..."
            try:
                with open(file_matrices, "r") as f:
                    distance_matrices = pickle.load(f)
            except (IOError, pickle.PickleError, AttributeError, EOFError, IndexError, KeyError):
                print "error loading existing matrices. forcing rebuild..."
                build_matrices = True

                # check if we need to extract pdb files
                if not os.path.isfile("%s/%s.pdb" % (dir_tmp, next(eval_data.models.iterkeys()))):
                    print "no extraced pdb files found. forcing extraction..."
                    extract = True

        # extract all pdb files to temp dir
        if extract:
            delete_glob(dir_tmp)
            print "extracing pdb files..."
            commands = []
            for f in files_assembly:
                commands.append(Command(["rna_extract", "-in:file:silent", f, "-in:file:silent_struct_type", "rna", "-out:prefix", dir_tmp + os.sep + splitext(basename(f))[0] + "_"], add_suffix="rosetta"))
            self.execute_commands(commands, threads=threads)

            for model in eval_data.models.itervalues():
                shutil.move(dir_tmp + os.sep + splitext(basename(model["source_file"]))[0] + "_" + model["tag_source"] + ".pdb", dir_tmp + os.sep + model["tag"] + ".pdb")

        # build residue distance matrix for all models
        if build_matrices:
            delete_glob(file_matrices)
            print "building matrices for %d models..." % len(eval_data.models)
            distance_matrices = {}
            x = 1
            for model in eval_data.models.itervalues():
                matrix = np.zeros([len(self.config["sequence"]), len(self.config["sequence"])])
                chain = pdbtools.parse_pdb("", dir_tmp + os.sep + model["tag"] + ".pdb")[0].child_list[0]

                for i in range(2, slen + 1):
                    for j in range(2, slen + 1):
                        if j == i:
                            continue
                        matrix[i - 1][j - 1] = chain[i]["P"] - chain[j]["P"]
                distance_matrices[model["tag"]] = matrix
                if x % 10 == 0:
                    sys.stdout.write("%d..." % x)
                    sys.stdout.flush()
                x += 1
            if x % 10 != 0:
                print "%x..." % x
            print "done."
            with open(file_matrices, "w") as f:
                pickle.dump(distance_matrices, f)

        # calculate scores
        print "calculating scores..."
        dca = dcatools.parse_dca_data(dca_prediction_filename)
        for model in eval_data.models.itervalues():
            score_custom = 0
            matrix = distance_matrices[model["tag"]]
            c_count = 0
            for contact in dca:
                c_count += 1
                if c_count > number_dca_predictions:
                    break
                found = False
                for i in range(-radius, radius + 1):
                    if contact.res1 + i < 1 or contact.res1 + i > slen:
                        continue
                    for j in range(-radius, radius + 1):
                        if contact.res2 + j < 1 or contact.res2 + j > slen:
                            continue
                        if abs(contact.res1 + i - contact.res2 - j) <= 3:
                            continue
                        if matrix[contact.res1 + i - 1][contact.res2 + j - 1] < threshold:
                            found = True
                            score_custom -= 1
                            break
                    if found:
                        break
            print "model: %-11s, score: %8.3f, score_custom: %3d" % (model["tag"], model["score"], score_custom)
            model["score_custom"] = score_custom

        eval_data.save_to_file(file_evaldata)

    def compare(self):
        """Print table of native rmsd values for all predictions"""
        self.check_config()
        if self.config["native_pdb_file"] is None:
            raise SimulationException("Cannot compare without native information.")
        print self.config["name"]

        # noinspection PyShadowingNames
        def print_comparison_line(cst_name, comparisons):
            print "  %-035s %05s %05s %05s" % (cst_name, comparisons[0], comparisons[1], comparisons[2])

        # loop over all different constraint sets
        # that is either files in the "constraints" directory, or directories under "predictions", and always "none"
        cst_names = ["none"]
        cst_names += [splitext(basename(cst_file))[0] for cst_file in glob.glob("constraints/*.cst")]
        cst_names += [basename(cst_dir) for cst_dir in glob.glob("predictions/*")]
        for cst_name in sorted(set(cst_names), key=natural_sort_key):
            skip = False
            try:
                eval_data = EvalData.load_from_cst(cst_name)
                if "rmsd_native" not in next(eval_data.models.itervalues()):
                    skip = True
            except SimulationException:
                skip = True
            if skip:
                print_comparison_line(cst_name, ["-", "-", "-"])
                continue
            comparisons = []
            for c in (1, 5, 10):
                if c > len(eval_data.clusters):
                    comparisons.append("-")
                    continue
                min_rmsd = 999
                for c2 in range(1, c + 1):
                    model = eval_data.clusters[c2]["primary_model"]
                    rmsd = eval_data.models[model]["rmsd_native"]
                    if rmsd < min_rmsd:
                        min_rmsd = rmsd
                comparisons.append("%.2f" % min_rmsd)
            print_comparison_line(cst_name, comparisons)

    # returns the path to the extracted pdb file
    def extract_pdb(self, constraints, model):
        """
        Extract PDB file of a model to the tmp directory.
        :param constraints: constraints selection
        :param model: tag of the model
        :return: path to the extracted PDB file
        """
        cst_name, cst_file = self.parse_cst_name_and_filename(constraints)
        dir_assembly = "predictions/%s/assembly" % cst_name
        dir_tmp = "predictions/%s/temp" % cst_name
        prefix = dir_tmp + os.path.sep + "tmp_"
        pdbfile = "%s/%s.pdb" % (dir_tmp, model["tag"])
        self.execute_command(["rna_extract", "-in:file:silent", "%s/%s" % (dir_assembly, model["source_file"]), "-out:prefix", prefix, "-tags", model["tag_source"]], add_suffix="rosetta", quiet=True)
        shutil.move("%s/tmp_%s.pdb" % (dir_tmp, model["tag_source"]), pdbfile)
        return pdbfile

    @staticmethod
    def get_models(constraints, model_list, kind="tag"):
        """
        Returns a selection of models with specific properties

        :param constraints: constraints selection
        :param model_list: see EvalData.print_models
        :param kind: see EvalData.get_models
        """
        return EvalData.load_from_cst(constraints).get_models(model_list=model_list, kind=kind)

    @staticmethod
    def print_models(constraints, model_list, kind="tag"):
        """
        Print a set of models

        :param constraints: constraints selection
        :param model_list: see EvalData.print_models
        :param kind: see EvalData.print_models
        """
        EvalData.load_from_cst(constraints).print_models(model_list=model_list, kind=kind)

    def extract_models(self, constraints, model_list, kind="tag"):
        """
        Extract PDB files of a set of models

        :param constraints: constraints selection
        :param model_list: see EvalData.print_models
        :param kind: see EvalData.print_models
        """
        for model in EvalData.load_from_cst(constraints).get_models(model_list=model_list, kind=kind):
            pdbfile = self.extract_pdb(constraints, model)
            print "Extracted: %s" % pdbfile

    # TODO: include mapping_mode?
    @staticmethod
    def _create_constraints_output_filename(input_filename, output_filename, cst_function, number_dca_predictions=None, output_format="%n_%f"):
        """
        Create a reasonable constraints output filename based on input filename and a formatting string.

        The formatting string can include placeholdes

        :param input_filename: input cst filename
        :param output_filename: fixed output filename or basename or None to automatically format using output_format
        :param cst_function: rosetta function
        :param number_dca_predictions: maximum number of DCA predictions to use
        :param output_format: formatting string containing placeholders: %f: rosetta function
                                                                         %n: basename of source file
                                                                         %d: number of dca predictions
        :return: output filename
        """
        source_basename = splitext(basename(input_filename))[0]
        cst_function_underscore = cst_function.replace(" ", "_")

        # chose a default filename, is none given
        if output_filename is None:
            output_filename = "constraints/%s.cst" % output_format
        # check if a path was given, if not put it in the constraints dir
        elif basename(output_filename) == output_filename:
            output_filename = "constraints/%s" % output_filename

        # make sure the filename always ends in .cst
        if splitext(output_filename)[1] != ".cst":
            output_filename += ".cst"

        # replace placeholders with actual values
        output_filename = output_filename.replace("%f", cst_function_underscore).replace("%n", source_basename).replace("%d", str(number_dca_predictions))
        return output_filename

    def parse_dca_filter_string(self, line):
        """
        Parse a text string and turn it tinto a list of DCA filters

        Multiple filters are separated by command and have the folling format:
        Format: <filter>:<arg>:<arg>:...

        Threshold filter: Lookup dca contacts in a PDB file, discard or keep contact depending on whether the contact is realized.
        Format: threshold:<n>:<cst>:<mode>:<moodel>
          n: threshold (< 0: keep below, > 0: keep above)
          cst: constraints selection to look up PDB file
          mode, model: model selection mode (see EvalData.get_models for details)
        Example: threshold:8.0:100rnaDCA_FADE_-100_26_20_-2_2:cluster:1,threshold:-6.0:100rnaDCA_FADE_-100_26_20_-2_2:cluster:1

        None filter: Empty filter
        Format: none

        :param line: string to parse
        :return: list of DcaFilter objects
        """
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
                eval_data = EvalData.load_from_cst(cst_name)
                model = eval_data.get_models([model_i], model_kind)[0]
                pdbfile = self.extract_pdb(cst_name, model)
                filter_pdb_chain = pdbtools.parse_pdb("", pdbfile)[0].child_list[0]
                dca_filter = dcatools.DcaFilterThreshold(filter_pdb_chain, abs(threshold), keep_below=threshold > 0, mode="minimum_heavy")
            filter_chain.append(dca_filter)

        return filter_chain

    def make_constraints(self, dca_prediction_filename="dca/dca.txt", output_filename=None, number_dca_predictions=100, cst_function="FADE -100 26 20 -2 2", filter_text=None, mapping_mode="allAtomWesthof"):
        """
        Create a set of constraints from a DCA prediction.

        :param dca_prediction_filename: DCA input file
        :param output_filename: fixed output filename or None to automatically create
        :param number_dca_predictions: maximum number of DCA predictions to use
        :param cst_function: rosetta function and parameters as text string
        :param filter_text: optional: List of DCA filters (see parse_dca_filter_string for details)
        :param mapping_mode: atom-to-atom mapping mode to use, supported values: "allAtomWesthof" or "pOnly"
        """
        self.check_config()
        output_filename = self._create_constraints_output_filename(dca_prediction_filename, output_filename, cst_function, number_dca_predictions, "%d%n_%f")
        print "Constraints creation:"
        print "    dca_prediction_filename: %s" % dca_prediction_filename
        print "    output_filename: %s" % output_filename
        print "    number_dca_predictions: %d" % number_dca_predictions
        print "    mode: %s" % mapping_mode
        print "    function: %s" % cst_function
        print "    filter: %s" % filter_text
        check_file_existence(dca_prediction_filename)

        # load dca contacts from file
        dca = dcatools.parse_dca_data(dca_prediction_filename)

        # filter dca data
        if filter_text is not None:
            print "Filtering dca data:"
            dca_filter_chain = self.parse_dca_filter_string(filter_text)
            dcatools.filter_dca_data(dca_data=dca, dca_filter_chain=dca_filter_chain)

        # create constraints
        print "Creating constraints:"

        cst_info = dcatools.build_cst_info_from_dca_contacts(dca, sequence=self.config["sequence"], mapping_mode=mapping_mode, cst_function=cst_function, number_dca_predictions=number_dca_predictions)

        # write to file
        with open(output_filename, "w") as out:
            for c in cst_info:
                out.write("%s %d %s %d %s\n" % (c[0], c[1], c[2], c[3], " ".join(map(str, c[4]))))

    def edit_constraints(self, constraints, output_filename=None, cst_function="FADE -100 26 20 -2 2"):
        """
        Edit an existing .cst file, replacing the rosetta function.
        :param constraints: constraints selection
        :param output_filename: fixed output filename or None to automatically create
        :param cst_function: rosetta function and parameters as text string
        """
        cst_name, input_filename = self.parse_cst_name_and_filename(constraints)
        output_filename = self._create_constraints_output_filename(input_filename, output_filename, cst_function)
        print "Constraints editing:"
        print "    input_filename:  %s" % input_filename
        print "    output_filename: %s" % output_filename
        print "    function:       %s" % cst_function
        check_file_existence(input_filename)

        if input_filename == output_filename:
            raise SimulationException("Input and output filename cannot be the same")

        pattern = re.compile(r"^(\S+\s+\d+\s+\S+\s+\d+)\s+(\S+)\s+.+$")
        with open(output_filename, "w") as output_fd:
            for line in utils.read_file_line_by_line(input_filename):
                m = pattern.match(line)
                output_fd.write("%s %s\n" % (m.group(1), cst_function))

    def print_status(self):
        """
        Print summary of all predictions and their current status.

        Output format contains the columns P, M, A, and E
        P: preparation step
        M: motif generation:
        A: assembly:
        E: evaluation

        If a step is completed, "X" is shown, "-" otherwise.
        For motif generation a "*" may be shown to indicate that models from a different
        set of constraints are used.
        """
        self.check_config()
        cst_names = ["none"]
        cst_names += [splitext(basename(cst_file))[0] for cst_file in glob.glob("constraints/*.cst")]
        cst_names += [basename(cst_dir) for cst_dir in glob.glob("predictions/*")]

        print "Status:"
        print "  cst                                 P M A E"
        print "  -------------------------------------------"

        for cst_name in sorted(set(cst_names), key=natural_sort_key):
            done_preparation = "-"
            done_motifs = "-"
            done_assembly = "-"
            done_evaluation = "-"

            dir_prediction = "predictions/%s" % cst_name
            dir_assembly = "predictions/%s/assembly" % cst_name
            dir_motifs = "predictions/%s/motifs" % cst_name
            dir_output = "predictions/%s/output" % cst_name

            if os.path.isfile(dir_assembly + "/assembly.cst"):
                done_preparation = "X"

            if os.path.isfile(dir_prediction + "/MOTIF_OVERRIDE"):
                done_motifs = "*"
            elif os.path.isfile(dir_motifs + "/motif1.out"):
                done_motifs = "X"

            if glob.glob(dir_assembly + "/*.out"):
                done_assembly = "X"

            if os.path.isfile(dir_output + "/cluster_1.pdb"):
                done_evaluation = "X"

            # noinspection PyShadowingNames
            def print_status_line(cst_name, status_list):
                line = "  %-035s" % cst_name
                line += " " + " ".join(status_list)
                print line

            print_status_line(cst_name, [done_preparation, done_motifs, done_assembly, done_evaluation])
