'''
Created on Aug 13, 2014

@author: sebastian
'''

import time
import shutil
import pprint
import glob
import struct
import re
import os
import subprocess
import sys
import errno
import string
import pickle
from . import dcatools
from . import pdbtools
from utils import readFileLineByLine
from os.path import splitext, basename, abspath


def checkFileExistence(path, alternativeErrorText=None):
    if not os.path.isfile(path):
        raise SimulationException(alternativeErrorText if alternativeErrorText is not None else "Cannot find file: %s" % path)


def checkDirExistence(path, alternativeErrorText=None):
    if not os.path.isdir(path):
        raise SimulationException(alternativeErrorText if alternativeErrorText is not None else "Cannot find directory: %s" % path)


def deleteGlob(pattern, print_notice=True):
    for f in glob.glob(pattern):
        if print_notice:
            print "deleting %s..." % (f)
        try:
            if os.path.isdir(f):
                shutil.rmtree(f, ignore_errors=True)
            else:
                os.remove(f)
        except:
            pass


def makeDirectory(directory):
    try:
        os.mkdir(directory)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(directory):
            pass
        else:
            raise


def mergeSilentFiles(target, sources):
    n = 0
    pattern_header = re.compile("^(?:SEQUENCE:|SCORE:\s+score).*")
    pattern_normal = re.compile("^(.*)S_\d+$")
    # if target already exists, see how many structures we have already
    if os.path.isfile(target):
        with open(target, "r") as t:
            for line in t:
                pass
            if 'line' in locals():
                m = re.match(".*S_(\d+)$", line)
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
    deleteGlob(sources, print_notice=True)

    # return total number of structures
    return n


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]


def parseCstFile(constraints_file):
    cst_info = []
    for line in readFileLineByLine(constraints_file):
        if len(line) > 6 and line[0] != '[':
            cols = string.split(line)
            atom_name1 = cols[0]
            res1 = int(cols[1])
            atom_name2 = cols[2]
            res2 = int(cols[3])
            cst_info.append( [ atom_name1, res1, atom_name2, res2, cols[4:]] )
    return cst_info

# TODO: Correcting atom names IS MOST LIKELY WRONG! It is commented out for now and can probably be removed completely.
def fixAtomNamesInCst(cst_info, sequence):
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
    pass


class Command(object):
    def __init__(self, command, add_suffix=None, dry_run=False, print_commands=True, stdin=None, quiet=False):
        self.command = command
        self.add_suffix = add_suffix
        self.dry_run = dry_run
        self.print_commands = print_commands
        self.stdin = stdin
        self.quiet = quiet


class RNAPrediction(object):
    '''
    Base class used for prediction simulation
    '''

    CONFIG_FILE = ".config"

    def loadConfig(self):
        if os.path.exists(RNAPrediction.CONFIG_FILE):
            c = open(RNAPrediction.CONFIG_FILE)
            self.config = pickle.load(c)
            c.close()
        else:
            self.config = {}

    def saveConfig(self):
        c = open(RNAPrediction.CONFIG_FILE, "w")
        pickle.dump(self.config, c)
        c.close()

    def printConfig(self):
        print "Simulation configuration:"
        print "    path: %s" % (abspath(os.getcwd()))
        if self.config:
            for key, value in sorted(self.config.iteritems()):
                if key == "motif_res_maps" or key == "motif_stem_sets" or key == "motifs" or key == "stems" or key == "evaluate":
                    print "    %s: ..." % (key)
                else:
                    print "    %s: %s" %(key, value)
        else:
            print "    No configuration found."

    # TODO: support different data types?
    def modifyConfig(self, key, value):
        if self.config.has_key(key):
            self.config[key] = None if value == "-" else value
        else:
            raise SimulationException("No such config entry: %s" % (key))


    def checkConfig(self):
        if not self.config:
            raise SimulationException("No config file found. Please run --prepare first!")

    def __init__(self, sysconfig, path):
        '''
        Create or load a prediction simulation
        '''
        try:
            os.chdir(path)
        except:
            raise SimulationException("Invalid basepath: %s" % (path))
        self.sysconfig = sysconfig
        self.loadConfig()

    def executeCommand(self, command, add_suffix=None, dry_run=False, print_commands=True, stdin=None, quiet=False):
        self.executeCommands([Command(command, add_suffix, dry_run, print_commands, stdin, quiet)])

    def executeCommands(self, commands, threads=1):
        ':type commands: list'

        processes = []
        stdout = None
        try:
            while True:
                while commands and len(processes) < threads:
                    c = commands.pop(0)
                    if c.add_suffix == "rosetta":
                        c.command[0] = self.sysconfig.rosetta_exe_path + c.command[0] + self.sysconfig.rosetta_exe_suffix
                    elif c.add_suffix == "gromacs":
                        c.command[0] = self.sysconfig.gromacs_exe_path + c.command[0] + self.sysconfig.gromacs_exe_suffix
                    if c.print_commands:
                        print " ".join(c.command)
                    if self.sysconfig.subprocess_buffsize is not None:
                        c.command = ["stdbuf", "-o", self.sysconfig.subprocess_buffsize] + c.command
                    if not c.dry_run:
                        stdout = None
                        if c.quiet and stdout is None:
                            stdout = open(os.devnull, "w")
                        stderr = stdout
                        try:
                            p = subprocess.Popen(c.command, stdin=(subprocess.PIPE if c.stdin != None else None), stdout=stdout, stderr=stderr)
                            # store original command in process
                            p.command = c.command
                            processes.append(p)
                            if c.stdin != None:
                                p.stdin.write(c.stdin)
                        except OSError, e:
                            print "lol1"
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
                try:
                    p.kill()
                except:
                    pass
            raise
        finally:
            if stdout is not None:
                stdout.close()

    def _make_tag_with_dashes(self, int_vector ):
        tag = []

        start_res = int_vector[0]
        for i in range( 1, len(int_vector)+1 ):
            if i==len( int_vector)  or  int_vector[i] != int_vector[i-1]+1:

                stop_res = int_vector[i-1]
                if stop_res > start_res:
                    tag += ["%d-%d" % (start_res, stop_res )]
                else:
                    tag += ["%d" % (stop_res)]

                if ( i < len( int_vector) ): start_res = int_vector[i]

        return tag

    def prepare(self, fasta_file="sequence.fasta", params_file="secstruct.txt", native_pdb_file=None, data_file=None, torsions_file=None, name=None):
        if name is None:
            name = os.path.basename(os.path.abspath(os.getcwd()))
        for f in [fasta_file, params_file, native_pdb_file, data_file, torsions_file]:
            if f is not None:
                checkFileExistence(f)
        makeDirectory("constraints")
        makeDirectory("preparation")
        makeDirectory("dca")
        makeDirectory("predictions")
        self.config["fasta_file"] = fasta_file
        self.config["params_file"] = params_file
        self.config["native_pdb_file"] = native_pdb_file
        self.config["data_file"] = data_file
        self.config["torsions_file"] = torsions_file
        self.config["name"] = name

        print
        print "Preparation:"

        # Read in files
        for line in readFileLineByLine(fasta_file):
            if line[0] == ">":
                continue
            self.config["sequence"] = line.lower()
            break
        print self.config["sequence"]
        numres = len( self.config["sequence"] )

        # Read in data information
        self.config["data_info"] = []
        if data_file != None:
            self.config["backbone_burial_info"] = []
            for line in readFileLineByLine(data_file):
                if len( line ) > 6 and line[:6]=='EXPOSE':
                    cols = string.split( line )
                    for i in range( len(cols)/3 ):
                        self.config["data_info"].append( [int( cols[ 3*i+1] ), cols[3*i+2],cols[3*i+3]] )
                if len( line ) > 15 and line[:15]=='BACKBONE_BURIAL':
                    cols = string.split( line )
                    for i in range( 1,len(cols) ):
                        self.config["backbone_burial_info"].append( int( cols[i] ) )

        pair_map = {}
        all_pairs = []

        complement = {'a':['u'], 'u':['a','g'], 'c':['g'], 'g':['c','u']};

        # Parse out stems
        basepair_mismatch = []
        cutpoints_original = []
        for line in readFileLineByLine(params_file):
            if line[:4] == 'STEM':
                cols = string.split( line )
                for i in range( len( cols )):
                    if cols[i] == 'PAIR':
                        #Offset to get to python numbering (starts with zero)
                        res1 = int(cols[i+1])-1
                        res2 = int(cols[i+2])-1
                        pair_map[ res1 ] = res2
                        pair_map[ res2 ] = res1
                        all_pairs.append( [res1,res2] )
                        assert ( self.config["sequence"][res1] in complement[ self.config["sequence"][res2] ] )
            elif line.count( '(' ) > 0:         #maybe dot/bracket notation (((...)))
                print line
                self.config["secstruc"] = line
                left_brackets = []
                for i in range( len(line) ):
                    if line[i] == '(':  left_brackets.append( i )
                    if line[i] == ')':
                        res1 = left_brackets[-1]
                        res2 = i
                        del( left_brackets[-1] )
                        pair_map[ res1 ] = res2
                        pair_map[ res2 ] = res1
                        all_pairs.append( [res1,res2] )
                        # check if the secondary structure watson-crick pair is allowed
                        if not self.config["sequence"][res1] in complement[ self.config["sequence"][res2] ]:
                            basepair_mismatch += [[res1, res2]]
                assert( len (left_brackets) == 0 )
            else:
                try:
                    cols = string.split( line )
                    res1 = int( cols[ 0 ] ) - 1
                    res2 = int( cols[ 1 ] ) - 1
                    pair_map[ res1 ] = res2
                    pair_map[ res2 ] = res1
                    all_pairs.append( [res1,res2] )
                    assert ( self.config["sequence"][res1] in complement[ self.config["sequence"][res2] ] )
                except:
                    continue

            if line[:13] == 'CUTPOINT_OPEN':
                cutpoints_original = map( lambda x:int(x), string.split( line[14:] ) )

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

        #print pair_map

        # Parse out stems
        already_in_stem = {}
        for i in range( numres ): already_in_stem[ i ] = 0

        self.config["stems"] = []
        for i in range( numres ):
            if pair_map.has_key( i ) and not already_in_stem[ i ]:  # In a base pair
                k = i
                stem_res = []

                stem_res.append( [k, pair_map[k]] )
                already_in_stem[ k ] = 1
                already_in_stem[ pair_map[k] ] = 1

                # Can we extend in one direction?
                while( pair_map.has_key( k + 1 ) and  pair_map[ k+1 ] == pair_map[ k ] - 1  and  not already_in_stem[k+1] ):
                    k += 1
                    stem_res.append( [k, pair_map[k]] )
                    already_in_stem[ k ] = 1
                    already_in_stem[ pair_map[k] ] = 1

                # Do not allow single WC base pairs.
                if ( len( stem_res ) <2 ):
                    print 'All stems must have length > 1 bp '
                    print stem_res
                    exit()
                self.config["stems"].append( stem_res )

        # Parse out motifs
        already_in_motif = {}
        for i in range( numres ): already_in_motif[ i ] = 0

        motif_cutpoints = []
        self.config["motif_res_maps"] = []
        self.config["motif_stem_sets"] = []
        self.config["motifs"] = []
        for i in range( numres ):

            if not already_in_motif[i] and ( not already_in_stem[ i ] or
                                             ( i > 0 and already_in_stem[i-1] and \
                                                   pair_map[i]+1 != pair_map[i-1] ) ):

                motif_res = []
                motif_stem_set = []
                cutpoints = []

                if ( i > 1 ):

                    # back up to beginning of stem.
                    k = i - 1
                    motif_stem = []

                    # first base pair.
                    motif_stem.append( [ k, pair_map[k] ] )
                    motif_res.append( k )
                    motif_res.append( pair_map[k] )

                    k-= 1
                    while k >= 0 and already_in_stem[k] and \
                            (pair_map[k]-1 == pair_map[k+1]):
                        motif_stem.append( [ k, pair_map[k] ] )
                        motif_res.append( k )
                        motif_res.append( pair_map[k] )
                        k -= 1
                    motif_stem_set.append( motif_stem )
                    k += 1
                    cutpoints.append( pair_map[k] )

                #print 'AFTER FIRST HELIX: ', motif_res

                k = i
                while ( k not in motif_res and k < numres):
                    # Move forward to next stem:
                    while ( k < numres and not already_in_stem[ k ] ):
                        if (already_in_motif[ k ] ):
                            print 'Hey cant deal with pseudoknots!'
                            exit()
                        motif_res.append( k )
                        k += 1

                    stem_start = k

                    if k >= numres:
                        cutpoints.append( k-1 )
                        break

                    if k in motif_res : break

                    motif_stem = []
                    motif_stem.append( [ k, pair_map[k] ] )
                    motif_res.append( k )
                    motif_res.append( pair_map[ k ] )
                    k+=1

                    while ( k < numres and already_in_stem[ k ] and \
                                pair_map[k-1] == pair_map[k]+1 and not k in motif_res):
                        motif_stem.append( [ k, pair_map[k] ] )
                        motif_res.append( k )
                        motif_res.append( pair_map[ k ] )
                        k += 1
                    motif_stem_set.append( motif_stem )
                    cutpoints.append( k-1 )

                    # Next non-helical part..
                    k = pair_map[ stem_start ] + 1

                    #print 'AFTER NEXT HELIX: ', motif_res

                motif_res.sort()

                motif_res_map = {}
                for k in range( len( motif_res ) ):
                    motif_res_map[ motif_res[k] ] = k
                    already_in_motif[ motif_res[k] ] = 1

                self.config["motifs"].append( motif_res )
                self.config["motif_stem_sets"].append( motif_stem_set )
                motif_cutpoints.append( cutpoints )
                self.config["motif_res_maps"].append( motif_res_map )
                #print 'CUTPOINTS ', cutpoints

        for i in range( len(self.config["stems"]) ):

            # Fasta
            tag = 'preparation/stem%d.fasta' % (i+1)
            fid = open( tag , 'w' )
            fid.write( '>'+tag+'\n')

            stem_res = self.config["stems"][i]
            stem_length = len( stem_res )

            for k in range( stem_length ):
                fid.write( self.config["sequence"][stem_res[k][0]] )
                #print stem_res[k][0]+1,
            for k in range( stem_length ):
                fid.write( self.config["sequence"][stem_res[stem_length-k-1][1]] )
                #print stem_res[stem_length-k-1][1]+1,

            #print
            fid.write('\n')
            fid.close()
            print 'Created: ', tag

            # pdb_file
            if native_pdb_file != None:
                command = ["pdbslice.py",
                          native_pdb_file,
                          "-segment",
                          "%d" % (stem_res[0][0]+1),
                          "%d" % (stem_res[-1][0]+1),
                          "%d" % (stem_res[-1][-1]+1),
                          "%d" % (stem_res[0][-1]+1),
                          "preparation/stem%d_" % (i+1)]
                self.executeCommand(command)
                native_pdb_file_subset =  'preparation/stem%d_%s' % (i+1, native_pdb_file )
                print 'Created: ', native_pdb_file_subset


        # Output motif jobs
        for i in range( len(self.config["motifs"]) ):

            # Fasta
            motif_fasta_file = 'preparation/motif%d.fasta' % (i+1)
            fid = open( motif_fasta_file , 'w' )
            fid.write( '>'+motif_fasta_file+'\n')

            motif_res = self.config["motifs"][i]
            motif_length = len( motif_res )

            for k in range( motif_length ):
                fid.write( self.config["sequence"][motif_res[k]] )
            fid.write('\n')
            fid.close()
            print 'Created: ', motif_fasta_file

            # params file
            motif_stem_set = self.config["motif_stem_sets"][ i ]
            motif_res_map = self.config["motif_res_maps"][ i ]
            motif_cutpoint = motif_cutpoints[ i ]

            motif_params_file = 'preparation/motif%d.params' % (i+1)
            fid = open( motif_params_file , 'w' )

            for k in range( len(motif_stem_set) ):
                motif_stem = motif_stem_set[ k ]
                fid.write( 'STEM   ' )
                for n in range( len( motif_stem ) ):
                    fid.write( '  PAIR %d %d W W A'  % \
                               ( motif_res_map[ motif_stem[n][0] ]+1,
                                 motif_res_map[ motif_stem[n][1] ]+1 ) )
                fid.write( '\n' )

                fid.write( 'OBLIGATE PAIR %d %d W W A \n\n'  % \
                               ( motif_res_map[ motif_stem[-1][0] ]+1,
                                 motif_res_map[ motif_stem[-1][1] ]+1 ) )

            motif_cutpoint.sort()

            #print motif_res
            #print motif_cutpoint

            if ( len( motif_cutpoint ) > 1 ):
                fid.write( 'CUTPOINT_OPEN ' )
                for k in range( len( motif_cutpoint ) ):
                    if motif_res_map[ motif_cutpoint[k] ] < (len( motif_res )-1):
                        fid.write( ' %d' % (motif_res_map[ motif_cutpoint[k] ]+1) )
            fid.write('\n')
            fid.close()
            print 'Created: ', motif_params_file

            # pdb_file
            if native_pdb_file != None:
                command = ["pdbslice.py",
                           native_pdb_file,
                           "-subset"]
                for k in range( motif_length ):
                    command += ["%d" % (motif_res[k]+1)]
                command += ["preparation/motif%d_" % (i+1)]
                self.executeCommand(command)
                native_pdb_file_subset =  'preparation/motif%d_%s' % (i+1, native_pdb_file )
                print 'Created: ', native_pdb_file_subset

            if data_file != None:
                motif_data_file = 'preparation/motif%d.data' % ( i+1 )
                fid_data = open( motif_data_file, 'w' )
                fid_data.write( 'EXPOSE' )
                for data in self.config["data_info"]:
                    if data[0]-1 in motif_res_map.keys():
                        fid_data.write( '   %d %s %s ' % (motif_res_map[data[0]-1]+1,data[1],data[2]) )
                fid_data.write('\n')

                if len( self.config["backbone_burial_info"] ) > 0:
                    fid_data.write( 'BACKBONE_BURIAL ' )
                    for k in self.config["backbone_burial_info"]:
                        if k-1 in motif_res_map.keys():
                            fid_data.write( ' %d' % (motif_res_map[ k-1 ] + 1) )
                    fid_data.write( '\n' )
                fid_data.close()
                print 'Created: ', motif_data_file


        if len( cutpoints_original ) > 0:
            fid.write( 'CUTPOINT_OPEN ' )
            for cutpoint in cutpoints_original:
                fid.write( ' %d'  % (cutpoint+1) )
            fid.write( '\n' )


        cutpoints  = []
        for i in range( len(self.config["motifs"]) ):

            motif_stem_set = self.config["motif_stem_sets"][ i ]

            motif_stem = motif_stem_set[ 0 ]

            possible_cutpoints =  [ motif_stem[ 0 ][ 0 ], motif_stem[ 1 ][ 1 ] ]
            possible_cutpoints.sort()
            #print possible_cutpoints
            if ( possible_cutpoints[0] not in cutpoints):
                cutpoints.append( possible_cutpoints[ 0 ] )


        params_file = "preparation/sequence.params"
        fid = open( params_file, 'w')

        if len( cutpoints ) > 0:
            fid.write( 'CUTPOINT_CLOSED ' )
            #cutpoints.sort()
            for cutpoint in cutpoints:
                fid.write( ' %d'  % (cutpoint+1) )
            fid.write( '\n' )

        #for cutpoint in cutpoints:
        #    fid.write( 'OBLIGATE   PAIR %d %d W W A\n' % (cutpoint+1, pair_map[cutpoint]+1) )

        for i in range( len(self.config["stems"]) ):
            stem_res = self.config["stems"][i]
            fid.write( 'STEM ')
            for k in range( len( stem_res )):
                fid.write( ' PAIR %d %d W W A ' % \
                               ( stem_res[k][ 0 ]+1, stem_res[k][ 1 ]+1 ) )
            fid.write('\n')
            fid.write( 'OBLIGATE PAIR %d %d W W A \n\n'  % \
                           ( stem_res[-1][0] + 1,\
                             stem_res[-1][1] + 1 ) )

        fid.close()
        print 'Created: ', params_file

        # default cutpoint constraints
        assemble_cst_file = "preparation/cutpoints.cst"
        fid = open( assemble_cst_file,'w')
        fid.write('[ atompairs ]\n')
        for cutpoint in cutpoints:
            fid.write( "O3'  %d  P     %d  HARMONIC  1.619  2.0\n" % \
                           ( cutpoint+1, cutpoint+2 ) )
        fid.close()
        print 'Created: ', assemble_cst_file


    def create_helices(self, dry_run=False, threads=1):
        self.checkConfig()
        commands = list()
        deleteGlob("preparation/stem*.out")
        for i in range(len(self.config["stems"])):
            command = ["rna_helix",
                       "-fasta", "preparation/stem%d.fasta" % (i+1),
                       "-out:file:silent","preparation/stem%d.out" % (i+1)]
            commands.append(Command(command, add_suffix="rosetta", dry_run=dry_run))
        self.executeCommands(commands, threads=threads)

        # rna_helix dumps <sequence>.pdb files in the working directory.
        # These can be useful for us, so move them to the preparation directory with a proper stemX.pdb filename.
        # TODO: What happens if there are multiple helices with the same sequence? In that case we need to move them
        #       out of the way before running the next command.
        #       Just disable multithreading here (not really needed with stem generation anyways)?
        for i in range(len(self.config["stems"])):
            sequence = [self.config["sequence"][j[0]] for j in self.config["stems"][i]]
            sequence += [self.config["sequence"][j[1]] for j in reversed(self.config["stems"][i])]
            shutil.move("%s.pdb" % ("".join(sequence)), "preparation/stem%d.pdb" % (i + 1))


    def _parseCstNameAndFilename(self, constraints):
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
            cst_file = "constraints/%s.cst" % (cst_name)
            checkFileExistence(cst_file)
        return (cst_name, cst_file)

    # prepare constraints files for motif generation and assembly
    # this is done separately to prevent a race condition when starting multiple parallel assembly jobs
    def prepareCst(self, constraints=None):
        self.checkConfig()
        cst_name, cst_file = self._parseCstNameAndFilename(constraints)
        print "Constraints preparation:"
        print "    constraints: %s" % (cst_name)
        if cst_file is None:
            print "Nothing to do."
            return

        dir_assembly = "predictions/%s/assembly" % (cst_name)
        dir_motifs = "predictions/%s/motifs" % (cst_name)
        makeDirectory("predictions/%s" % (cst_name))
        makeDirectory(dir_assembly)
        makeDirectory(dir_motifs)

        # load constraints information
        cst_info = parseCstFile(cst_file)
        # probably wrong, commented out for now!
        #cst_info = fixAtomNamesInCst(cst_info, self.config["sequence"])

        # extract relevant constraints for motif generation from global cst file
        for i in range(len(self.config["motifs"])):
            motif_res_map = self.config["motif_res_maps"][i]
            motif_cst_file = '%s/motif%d.cst' % (dir_motifs, i+1)
            fid_cst = open(motif_cst_file, 'w')
            fid_cst.write( '[ atompairs ]\n' )
            for cst in cst_info:
                if cst[1]-1 in motif_res_map.keys() and cst[3]-1 in motif_res_map.keys():
                    fid_cst.write( '%s %d %s %d %s\n' % (cst[0], motif_res_map[cst[1]-1]+1,cst[2],motif_res_map[cst[3]-1]+1," ".join(map(str, cst[4]))) )
            fid_cst.close()
            print 'Created: ', motif_cst_file

        # assembly constraints: default harmonic cutpoints + tertiary constraints
        assembly_cst = "%s/assembly.cst" % (dir_assembly)
        shutil.copy("preparation/cutpoints.cst", assembly_cst)
        with open(assembly_cst, "a") as ft:
            for cst in cst_info:
                ft.write('%s %d %s %d %s \n' % (cst[0], cst[1], cst[2], cst[3], " ".join(map(str, cst[4]))))
        print 'Created: ', assembly_cst


    def create_motifs(self, nstruct=50000, cycles=20000, dry_run=False, seed=None, use_native_information=False, threads=1, constraints=None):
        self.checkConfig()
        cst_name, cst_file = self._parseCstNameAndFilename(constraints)
        print "Motif creation configuration:"
        print "    constraints: %s" % (cst_name)
        print "    cycles: %s" % (cycles)
        print "    nstruct: %s" % (nstruct)
        print "    dry_run: %s" % (dry_run)
        print "    random_seed: %s" % (seed)
        print "    threads: %s" % (threads)

        dir_motifs = "predictions/%s/motifs" % (cst_name)
        makeDirectory("predictions/%s" % (cst_name))
        makeDirectory(dir_motifs)

        n_motifs = len(self.config["motifs"])

        # check if motif constraints were created correctly
        if cst_file is not None:
            for i in range(n_motifs):
                checkFileExistence("%s/motif%d.cst" % (dir_motifs, i+1), "Motif cst files not found. Please run --prepare-cst step!")

        # merge all motifs abd check what we have so far
        completed = {}
        for i in range(n_motifs):
            completed[i] = mergeSilentFiles("%s/motif%d.out" % (dir_motifs, i + 1), "%s/motif%d_*.out" % (dir_motifs, i + 1))

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
                       "-fasta", "preparation/motif%d.fasta" % (i+1),
                       "-params_file", "preparation/motif%d.params" % (i+1),
                       "-cycles", "%d" % (cycles),
                       "-mute", "all",
                       "-close_loops",
                       "-close_loops_after_each_move",
                       "-minimize_rna"]

            if cst_file is not None:
                command += ["-cst_file", '%s/motif%d.cst' % (dir_motifs, i+1)]
            if self.config["native_pdb_file"] != None and use_native_information:
                command += ["-native", "preparation/motif%d_%s" % (i+1, self.config["native_pdb_file"])]
            if self.config["data_file"] != None:
                command += ["-data_file", "preparation/motif%d.data" % (i+1)]
            if self.config["torsions_file"] != None:
                command += ["-vall_torsions", "preparation/motif%d.torsions" % (i+1)]

            which_stems = []
            stem_chunk_res = []
            motif_stem_set = self.config["motif_stem_sets"][i]
            motif_res_map = self.config["motif_res_maps"][i]
            for k in range( len(motif_stem_set) ):
                motif_stem = motif_stem_set[k]
                #need to find in stems
                for n in range( len( self.config["stems"] )  ):
                    stem = self.config["stems"][n]
                    found_match = 0
                    for q in range( len( stem ) ) :
                        if motif_stem[0][0] in stem[q]:
                            found_match = 1
                            break
                    if found_match: break
                which_stems.append( n )

                for q in range( len( stem ) ):
                    stem_chunk_res.append( motif_res_map[ stem[q][0] ]+1 )
                for q in range( len( stem ) ):
                    stem_chunk_res.append( motif_res_map[ stem[ len(stem) - q - 1][1] ]+1 )

            command += ["-in:file:silent_struct_type", "rna",
                        "-in:file:silent"]
            for n in which_stems:
                command += ["preparation/stem%d.out" % (n+1)]

            command += ["-in:file:input_res"]
            command += self._make_tag_with_dashes(stem_chunk_res)

            for j in range(threads):
                if structs_threads[j] == 0:
                    continue
                command_full = command + ["-out:file:silent", "%s/motif%d_%d.out" % (dir_motifs, i + 1, j + 1),
                                          "-nstruct", "%d" % (structs_threads[j])]
                if seed != None:
                    command_full += ["-constant_seed", "-jran", "%d" % (seed)]
                    seed += 1

                commands.append(Command(command_full, add_suffix="rosetta", dry_run=dry_run))
        self.executeCommands(commands, threads=threads)

        # merge motifs
        for i in range(n_motifs):
            mergeSilentFiles("%s/motif%d.out" % (dir_motifs, i + 1), "%s/motif%d_*.out" % (dir_motifs, i + 1))


    # TODO: When documenting later, explain that with assemble, nstruct is used for each single thread, while with createMotifs, it is distributed.
    def assemble(self, nstruct=50000, cycles=20000, constraints=None, dry_run=False, seed=None, use_native_information=False, threads=1, motifsOverride=None):
        self.checkConfig()
        cst_name, cst_file = self._parseCstNameAndFilename(constraints)
        print "Assembly configuration:"
        print "    constraints: %s" % (cst_name)
        print "    cycles: %s" % (cycles)
        print "    nstruct: %s" % (nstruct)
        print "    dry_run: %s" % (dry_run)
        print "    random_seed: %s" % (seed)

        # Check if a different set of motifs should be used.
        dir_assembly = "predictions/%s/assembly" % (cst_name)
        if motifsOverride is not None:
            override_name, override_file = self._parseCstNameAndFilename(motifsOverride)
            dir_motifs = "predictions/%s/motifs" % (override_name)
        else:
            dir_motifs = "predictions/%s/motifs" % (cst_name)
        print "    using motifs: %s" % (dir_motifs)

        checkFileExistence("%s/motif1.out" % (dir_motifs))
        checkFileExistence("preparation/cutpoints.cst")
        checkDirExistence("predictions/%s" % (cst_name))

        makeDirectory(dir_assembly)

        # If an assembly was run before, perform sanity check to see if the motifs used last time match the ones that are going to be used now.
        file_assembly_motif_override = "%s/MOTIF_OVERRIDE" % (dir_assembly)
        if len(glob.glob("%s/*.out" % (dir_assembly))):
            dir_motifs_last = "predictions/%s/motifs" % (cst_name)
            if os.path.isfile(file_assembly_motif_override):
                with open(file_assembly_motif_override, "r") as f:
                    dir_motifs_last = f.read().strip()
            if dir_motifs_last != dir_motifs:
                raise SimulationException("Error: Motifs used last time don't match the current ones! %s <-> %s." % (dir_motifs_last, dir_motifs))

        # Check if assembly constraints are available.
        # Do NOT always try to create them here because that would cause a race condition when multiple assembly jobs are started in parallel!
        assembly_cst = "%s/assembly.cst" % (dir_assembly)
        if cst_file is not None:
            checkFileExistence(assembly_cst, "Assembly cst file '%s' not found. Please run --prepare-cst step!" % (assembly_cst))

        # In case of motifs override, put a hint .txt file in the assembly directory to indicate that non-matching motifs were used.
        if motifsOverride is not None:
            with open(file_assembly_motif_override, "w") as f:
                f.write(dir_motifs)

        commands = list()

        command = ["rna_denovo",
                   "-minimize_rna",
                   "-fasta", self.config["fasta_file"],
                   "-in:file:silent_struct_type", "binary_rna",
                   "-cycles", "%d" % (cycles),
                   "-nstruct", "%d" % (nstruct),
                   "-params_file", "preparation/sequence.params",
                   "-close_loops",
                   "-in:file:silent"]

        for i in range(len(self.config["stems"])):
            command += ["preparation/stem%d.out" % (i+1)]

        for i in range(len(self.config["motifs"])):
            command += ["%s/motif%d.out" % (dir_motifs, i+1)]

        if cst_file is not None:
            command += ["-cst_file", assembly_cst]

        chunk_res = []
        for n in range( len(self.config["stems"])  ):
            stem = self.config["stems"][n]
            for q in range( len( stem ) ):        chunk_res.append(stem[q][0] + 1)
            for q in range( len( stem ) ):        chunk_res.append(stem[ len(stem) - q - 1][1] + 1)

        for n in range( len(self.config["motifs"]) ):
            motif_res = self.config["motifs"][n]
            for m in motif_res: chunk_res.append(m+1)

        command += ["-in:file:input_res"]
        command += self._make_tag_with_dashes(chunk_res)

        if self.config["native_pdb_file"] != None and use_native_information:
            command += ["-native", self.config["native_pdb_file"]]


        if self.config["torsions_file"] != None:
            command += ["-vall_torsions", self.config["torsions_file"]]
        if self.config["data_file"] != None:
            command += ["-data_file", self.config["data_file"]]

        for j in range(threads):
            # In case no seed is specified, we do what rosetta does to get a random number, and seed it with that.
            # This way we know the number beforehand and can choose an appropriate output filename.
            # In case a seed was specified, we increment it here with the thread number.
            if seed is None:
                seed_used = struct.unpack("=i", os.urandom(4))[0]
            else:
                seed_used = seed + j

            command_full = command + ["-out:file:silent", "%s/assembly_%d.out" % (dir_assembly, seed_used),
                                      "-constant_seed", "-jran", "%d" % (seed_used)]

            commands.append(Command(command_full, add_suffix="rosetta", dry_run=dry_run))

        self.executeCommands(commands, threads=threads)


    def extract(self, constraints=None):
        self.checkConfig()
        cst_name, cst_file = self._parseCstNameAndFilename(constraints)
        print "Extraction:"
        print "    constraints: %s" % (cst_name)

        dir_assembly = "predictions/%s/assembly" % (cst_name)
        dir_output = "predictions/%s/output" % (cst_name)
        dir_tmp = "predictions/%s/temp" % (cst_name)
        file_evaldata = "%s/evaldata.dat" % (dir_output)

        if not os.path.isdir("predictions/%s" % (cst_name)):
            raise SimulationException("No prediction directory for constraint '%s' found! Maybe you should assemble first?" % (cst_name))

        files_assembly = glob.glob("%s/assembly_*.out" % (dir_assembly))
        if len(files_assembly) == 0:
            raise SimulationException("No assembly files for constraint '%s' found! Maybe you should assemble first?" % (cst_name))

        # cleanup
        deleteGlob("%s/assembly_p.pdb" % (dir_assembly))
        deleteGlob(dir_output)
        deleteGlob(dir_tmp)
        makeDirectory(dir_output)
        makeDirectory(dir_tmp)

        # create dict to store evaluation data
        evalData = {"models": {}, "clusters": {}}

        # model counter
        model = 0

        # loop over all out files matching the constraint
        for f in sorted(files_assembly, key=natural_sort_key):
            print "  processing rosetta silent file: %s..." % (f)

            description = splitext(basename(f))[0]

            # read out file and store a dict of the scores and the full score line
            scores = dict()
            score_regex = re.compile("^SCORE:\s+([-0-9.]+)\s.*(S_\d+)$")
            for line in readFileLineByLine(f):
                m = score_regex.match(line)
                if m:
                    scores[m.group(2)]={"line": line.rstrip(), "score": float(m.group(1))}


            # delete any pdb files if existing
            deleteGlob("S_*.pdb")

            # extract pdb files
            command = ["rna_extract", "-in:file:silent", abspath(f), "-in:file:silent_struct_type", "rna"]
            owd = os.getcwd()
            try:
                os.chdir(dir_tmp)
                self.executeCommand(command, add_suffix="rosetta")
            finally:
                os.chdir(owd)

            # loop over all extracted pdb files
            for fe in sorted(glob.glob("%s/S_*.pdb" % (dir_tmp))):
                model += 1
                name = basename(fe)[:-4]
                remark = "%s %s" % (description, scores[name]["line"])

                # extract all P atoms from the pdb file
                p_only = pdbtools.extractPOnly(fe)

                # append p only model to trajectory
                pdbtools.writePdb("%s/assembly_p.pdb" % (dir_assembly), data=p_only, model=model, remark=remark, append=True)

                # write p only model to temp file
                pdbtools.writePdb("%s/%09d_p.pdb" % (dir_tmp, model), data=p_only, model=model, remark=remark)

                # move original pdb to temp directory
                shutil.move(fe, "%s/%09d.pdb" % (dir_tmp, model))

                # create evaluation dict for model
                evalData["models"][model] = {"source_file": f, "source_name": name, "score": scores[name]["score"]}

        # save evaluation data to file
        with open(file_evaldata, "w") as f:
            pickle.dump(evalData, f)


    # TODO: set the cutoff value back to 0.40? It was set to 0.41 because the old bash script used integer comparison and even 0.409 was treated as 0.40
    def evaluate(self, constraints=None, cluster_limit=10, cluster_cutoff=0.41):
        self.checkConfig()
        cst_name, cst_file = self._parseCstNameAndFilename(constraints)
        print "Evaluation configuration:"
        print "    cluster_limit: %s" % (cluster_limit)
        print "    cluster_cutoff: %s" % (cluster_cutoff)
        print "    constraints: %s" % (cst_name)

        dir_assembly = "predictions/%s/assembly" % (cst_name)
        dir_output = "predictions/%s/output" % (cst_name)
        dir_tmp = "predictions/%s/temp" % (cst_name)
        file_evaldata = "%s/evaldata.dat" % (dir_output)

        # load evaluation data and check if extraction step was run before
        try:
            with open(file_evaldata, "r") as f:
                evalData = pickle.load(f)
            evalData["models"][1]["score"]
        except KeyError:
            raise SimulationException("Broken data file %s. Try running --extract again" % (file_evaldata))
        except:
            raise SimulationException("Data file %s for constraint not found. Did you forget to run --extract?" % (file_evaldata))

        if not os.path.isdir(dir_tmp) or not os.path.isfile("%s/%09d_p.pdb" % (dir_tmp, len(evalData["models"]))):
            raise SimulationException("No extracted pdb for constraint '%s' found. Did you delete temp/ files?" % (cst_name))

        # clear old evaluation clusters
        deleteGlob(dir_output)
        makeDirectory(dir_output)

        evalData["clusters"] = {}
        for m in evalData["models"].iteritems():
            m[1].pop("native_rmsd", None)
            m[1].pop("cluster", None)
            m[1].pop("rmsd_to_primary", None)


        # calculate native rmsd values for all models if native pdb available
        filename_rmsd = "%s/rmsd.xvg" % (dir_tmp)
        if self.config["native_pdb_file"] != None:
            # create native_p.pdb if needed
            native_p_only = "%s_p.pdb" % (self.config["native_pdb_file"][:-4])
            if not os.path.exists(native_p_only):
                print "  creating p only native pdb file: %s" % (native_p_only)
                pdbtools.writePdb(native_p_only, data=pdbtools.extractPOnly(self.config["native_pdb_file"]), model=1, remark="native p only")
            print "  caluculating rmsd values to native structure for all models..."
            sys.stdout.write("    ")
            self.executeCommand(["g_rms", "-quiet", "-s", "native_p.pdb", "-f", "%s/assembly_p.pdb" % (dir_assembly), "-o", filename_rmsd], add_suffix="gromacs", stdin="1\n1\n", quiet=True)
            for line in readFileLineByLine(filename_rmsd):
                if not re.match(r"^[\s\d-]", line):
                    continue
                model, native_rmsd = line.split()
                evalData["models"][int(float(model))]["native_rmsd"] = float(native_rmsd)
            deleteGlob(filename_rmsd, print_notice=False)

        # cluster counter
        cluster = 0

        log = open("%s/clusters.log" % (dir_output), "w")
        print "  sorting models into clusters..."

        # do we have native rmsd information?
        use_native = self.config["native_pdb_file"] != None
        # loop over all models sorted by their score
        for model, data in sorted(evalData["models"].items(), key=lambda x: x[1]["score"]):
            filename_pdb = "%s/%09d.pdb" % (dir_tmp, model)
            filename_pdb_p = "%s/%09d_p.pdb" % (dir_tmp, model)

            # check if the current structure matches a cluster
            matches_cluster = 0
            rmsd_to_cluster_primary = cluster_cutoff
            for c in range(cluster):
                # calculate rmsd between cluster and pdb
                self.executeCommand(["g_rms", "-quiet", "-s", "%s/cluster_%d_p.pdb" % (dir_output, c + 1), "-f", filename_pdb_p, "-o", filename_rmsd], add_suffix="gromacs", stdin="1\n1\n", quiet=True, print_commands=False)
                for line in readFileLineByLine(filename_rmsd):
                    pass
                new_rmsd = float(line.split()[1])
                if new_rmsd < rmsd_to_cluster_primary:
                    rmsd_to_cluster_primary = new_rmsd
                    matches_cluster = c + 1
                deleteGlob(filename_rmsd, print_notice=False)

            if matches_cluster == 0:
                cluster = cluster + 1
                shutil.copyfile(filename_pdb, "%s/cluster_%d.pdb" % (dir_output, cluster))
                shutil.copyfile(filename_pdb_p, "%s/cluster_%d_p.pdb" % (dir_output, cluster))
                evalData["clusters"][cluster] = {"primary_model": model}
                evalData["models"][model]["cluster"] = cluster
                evalData["models"][model]["rmsd_to_primary"] = 0
                output = "new cluster: %02s, model: %05s, native_rmsd: %.7f, score: %.3f" % (cluster,
                                                                                           model,
                                                                                           evalData["models"][model]["native_rmsd"] if use_native else -1,
                                                                                           evalData["models"][model]["score"])
                print "    %s" % (output)
                log.write(output + "\n")
            else:
                evalData["models"][model]["cluster"] = matches_cluster
                evalData["models"][model]["rmsd_to_primary"] = rmsd_to_cluster_primary
                output = "    cluster: %02s, model: %05s, native_rmsd: %.7f, score: %.3f, rmsd_to_cluster_primary: %.7f" % (matches_cluster,
                                                                                                                        model,
                                                                                                                        evalData["models"][model]["native_rmsd"] if use_native else -1,
                                                                                                                        evalData["models"][model]["score"],
                                                                                                                        evalData["models"][model]["rmsd_to_primary"])
                print "    %s" % (output)
                log.write(output + "\n")

            log.flush()

            if cluster == cluster_limit:
                print "    maximum number of clusters reached. stopping..."
                break

        # save evaluation data
        with open(file_evaldata, "w") as f:
            pickle.dump(evalData, f)


    def compare(self):
        self.checkConfig()
        if self.config["native_pdb_file"] is None:
            raise SimulationException("Cannot compare without native information.")
        print self.config["name"]

        def printComparisonLine(cst_name, comparisons):
            print "  %-035s %05s %05s %05s" % (cst_name, comparisons[0], comparisons[1], comparisons[2])

        # loop over all different constraint sets
        # that is either files in the "constraints" directory, or directories under "predictions", and always "none"
        cst_names = ["none"]
        cst_names += [splitext(basename(cst_file))[0] for cst_file in glob.glob("constraints/*.cst")]
        cst_names += [basename(cst_dir) for cst_dir in glob.glob("predictions/*")]
        for cst_name in sorted(set(cst_names), key=natural_sort_key):
            try:
                with open("predictions/%s/output/evaldata.dat" % (cst_name), "r") as f:
                    evalData = pickle.load(f)
                    evalData["models"][1]["native_rmsd"]
            except:
                printComparisonLine(cst_name, ["-", "-", "-"])
                continue
            comparisons = []
            for c in (1, 5, 10):
                if c > len(evalData["clusters"]):
                    comparisons.append("-")
                    continue
                min_rmsd = 999
                for c2 in range(1, c + 1):
                    model = evalData["clusters"][c2]["primary_model"]
                    rmsd = evalData["models"][model]["native_rmsd"]
                    if rmsd < min_rmsd:
                        min_rmsd = rmsd
                comparisons.append("%.2f" % (min_rmsd * 10))
            printComparisonLine(cst_name, comparisons)


    # retrieve a list for models by kind:
    # "model": internal numbering / order in which they were extracted from rosetta data
    # "top": number by score
    # "cluster": model representing the nth cluster
    def getModels(self, constraints, numbers, kind="model"):
        cst_name, cst_file = self._parseCstNameAndFilename(constraints)

        dir_output = "predictions/%s/output" % (cst_name)
        dir_tmp = "predictions/%s/temp" % (cst_name)

        with open("predictions/%s/output/evaldata.dat" % (cst_name), "r") as f:
            evalData = pickle.load(f)

        # initialize list
        models_sorted = None

        results = []

        for number in numbers:
            # decide which model we want
            if kind == "cluster":
                model = evalData["clusters"][number]["primary_model"]
            elif kind == "model":
                model = number
            elif kind == "top":
                # do we need to create a sorted list?
                if models_sorted is None:
                    models_sorted = sorted(evalData["models"].items(), key=lambda x: x[1]["score"])
                model = models_sorted[number - 1][0]

            m = evalData["models"][model]

            # add the model number and pdb paths to the model dict
            m["model"] = model
            m["pdbfile"] = "%s/%09d.pdb" % (dir_tmp, model)
            m["pdbfile_p"] = "%s/%09d_p.pdb" % (dir_tmp, model)
            results.append(m)

        return results


    # create a reasonable output filename from
    # output format uses placeholders for input name, number of predictions, and function
    # TODO: include mappingMode?
    def _createConstraintsOutputFilename(self, inputFileName, outputFileName, cstFunction, numberDcaPredictions=None, outputFormat="%n_%f"):
        sourceBasename = splitext(basename(inputFileName))[0]
        cstFunctionUnderscore = cstFunction.replace(" ", "_")

        # chose a default filename, is none given
        if outputFileName is None:
            outputFileName = "constraints/%s.cst" % (outputFormat)
        # check if a path was given, if not put it in the constraints dir
        elif basename(outputFileName) == outputFileName:
            outputFileName = "constraints/%s" % (outputFileName)

        # make sure the filename always ends in .cst
        if splitext(outputFileName)[1] != ".cst":
            outputFileName = outputFileName + ".cst"

        # replace placeholders with actual values
        outputFileName = outputFileName.replace("%f", cstFunctionUnderscore).replace("%n", sourceBasename).replace("%d", str(numberDcaPredictions))
        return outputFileName


    def makeConstraints(self, dcaPredictionFileName="dca/dca.txt", outputFileName=None, numberDcaPredictions=100, cstFunction="FADE -100 26 20 -2 2", filterText=None, mappingMode="allAtomWesthof"):
        # TODO: make this dependable on the prepare step? Or separate the whole constraints creation into an independent application?
        self.checkConfig()
        outputFileName = self._createConstraintsOutputFilename(dcaPredictionFileName, outputFileName, cstFunction, numberDcaPredictions, "%d%n_%f")
        print "Constraints creation:"
        print "    dcaPredictionFileName: %s" % (dcaPredictionFileName)
        print "    outputFileName: %s" % (outputFileName)
        print "    numberDcaPredictions: %d" % (numberDcaPredictions)
        print "    mode: %s" % (mappingMode)
        print "    function: %s" % (cstFunction)
        print "    filter: %s" % (filterText) # TODO: make filters a class and add a __str__ representation
        checkFileExistence(dcaPredictionFileName)

        # load dca contacts from file
        dca = dcatools.parseDcaData(dcaPredictionFileName)

        # filter dca data
        if filterText is not None:
            print "Filtering dca data:"
            dcaFilterChain = dcatools.parseFilterLine(filterText, self)
            dcatools.filterDcaData(dcaData=dca, dcaFilterChain=dcaFilterChain)

        # create constraints
        print "Creating constraints:"

        cst_info = dcatools.buildCstInfoFromDcaContacts(dca, sequence=self.config["sequence"], mappingMode=mappingMode, cstFunction=cstFunction, numberDcaPredictions=numberDcaPredictions)

        # write to file
        with open(outputFileName, "w") as out:
            for c in cst_info:
                out.write("%s %d %s %d %s\n" % (c[0], c[1], c[2], c[3], " ".join(map(str, c[4]))))


    def editConstraints(self, constraints, outputFileName=None, cstFunction="FADE -100 26 20 -2 2"):
        ':type cstFunction:string'
        cst_name, inputFileName = self._parseCstNameAndFilename(constraints)
        outputFileName = self._createConstraintsOutputFilename(inputFileName, outputFileName, cstFunction)
        print "Constraints editing:"
        print "    inputFileName:  %s" % (inputFileName)
        print "    outputFileName: %s" % (outputFileName)
        print "    function:       %s" % (cstFunction)
        checkFileExistence(inputFileName)

        if inputFileName == outputFileName:
            raise SimulationException("Input and output filename cannot be the same")

        pattern = re.compile(r"^(\S+\s+\d+\s+\S+\s+\d+)\s+(\S+)\s+.+$")
        with open(outputFileName, "w") as outputFd:
            for line in readFileLineByLine(inputFileName):
                m = pattern.match(line)
                outputFd.write("%s %s\n" % (m.group(1), cstFunction))
