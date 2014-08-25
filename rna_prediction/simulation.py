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
import ConfigParser
import pickle
from os.path import expanduser, splitext, basename, abspath


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]


class SimulationException(Exception):
    pass


class SysConfig(object):
    SYSCONFIG_FILE = expanduser("~/.rna_predict")

    def __init__(self):
        '''
        Load system configuration
        '''
        self.loadSysConfig()

    def loadSysConfig(self):
        #defaults
        self.sysconfig = {"rosetta_exe_path": "",
                          "rosetta_exe_suffix": ".linuxgccrelease",
                          "gromacs_exe_path": "",
                          "gromacs_exe_suffix": ""}

        config = ConfigParser.RawConfigParser()
        config.read(SysConfig.SYSCONFIG_FILE)
        if config.has_section("rosetta"):
            if config.has_option("rosetta", "exe_path"):
                self.sysconfig["rosetta_exe_path"] = re.sub("/+$", "", config.get("rosetta", "exe_path")) + "/"
            if config.has_option("rosetta", "exe_suffix"):
                self.sysconfig["rosetta_exe_suffix"] = config.get("rosetta", "exe_suffix")
        if config.has_section("gromacs"):
            if config.has_option("gromacs", "exe_path"):
                self.sysconfig["gromacs_exe_path"] = re.sub("/+$", "", config.get("gromacs", "exe_path")) + "/"
            if config.has_option("gromacs", "exe_suffix"):
                self.sysconfig["gromacs_exe_suffix"] = config.get("gromacs", "exe_suffix")

    def checkSysConfig(self):
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        progs = [self.sysconfig["rosetta_exe_path"] + "rna_helix" + self.sysconfig["rosetta_exe_suffix"],
                 self.sysconfig["gromacs_exe_path"] + "g_rms" + self.sysconfig["gromacs_exe_suffix"]]


        def isOk(prog):
            fpath, fname = os.path.split(prog)
            if fpath:
                if is_exe(prog):
                    return True
            else:
                for path in os.environ["PATH"].split(os.pathsep):
                    path = path.strip('"')
                    exe_file = os.path.join(path, prog)
                    if is_exe(exe_file):
                        return True

        fail = []
        success = []
        for prog in progs:
            if isOk(prog):
                success += [prog]
            else:
                fail += [prog]
        return {    "fail": fail, "success": success}

    def printSysConfig(self):
        print "System configuration:"
        for key, value in sorted(self.sysconfig.iteritems()):
            print "    %s: %s" %(key, value)


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
        self.config[key] = value


    def checkConfig(self):
        if not self.config:
            raise SimulationException("No config file found. Please run --prepare first!")

    def checkFileExistence(self, path):
        if not os.path.isfile(path):
            raise SimulationException("Cannot find file: %s" % path)

    def __init__(self, sysconfig, path):
        '''
        Create or load a prediction simulation
        '''
        try:
            os.chdir(path)
        except:
            raise SimulationException("Invalid basepath: %s" % (path))
        self.sysconfig = sysconfig.sysconfig
        self.loadConfig()

    def makeDirectory(self, directory):
        try:
            os.mkdir(directory)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(directory):
                pass
            else:
                raise

    def executeCommand(self, command, add_suffix=None, dry_run=False, print_commands=True, stdin=None, quiet=False):
        if add_suffix == "rosetta":
            command[0] = self.sysconfig["rosetta_exe_path"] + command[0] + self.sysconfig["rosetta_exe_suffix"]
        elif add_suffix == "gromacs":
            command[0] = self.sysconfig["gromacs_exe_path"] + command[0] + self.sysconfig["gromacs_exe_suffix"]
        if print_commands:
            print " ".join(command)
        if not dry_run:
            stdout = None
            if quiet:
                stdout = open(os.devnull, "w")
            stderr = stdout
            try:
                p = subprocess.Popen(command, stdin=(subprocess.PIPE if stdin != None else None), stdout=stdout, stderr=stderr)
                if stdin != None:
                    p.communicate(input=stdin)
                p.wait()
                if p.returncode > 0:
                    raise SimulationException("Non-zero return code from executed command: %s" % " ".join(command))
            except OSError, e:
                raise SimulationException("Failed to execute command: %s, Reason: %s" % (" ".join(command), e))
            finally:
                if quiet:
                    stdout.close()

    def make_tag_with_dashes(self, int_vector ):
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

    def prepare(self, fasta_file="sequence.fasta", params_file="secstruct.txt", native_pdb_file=None, data_file=None, cst_file=None, torsions_file=None, name=None):
        self.makeDirectory("stems_and_motifs")
        self.makeDirectory("assembly")
        self.makeDirectory("constraints")
        self.makeDirectory("output")
        self.makeDirectory("temp")
        self.config["fasta_file"] = fasta_file
        self.config["params_file"] = params_file
        self.config["native_pdb_file"] = native_pdb_file
        self.config["data_file"] = data_file
        self.config["cst_file"] = cst_file
        self.config["torsions_file"] = torsions_file
        self.config["name"] = name

        print
        print "Preparation:"

        # Read in files
        lines = open( fasta_file ).readlines()
        self.config["sequence"] = lines[1][:-1]
        print self.config["sequence"]
        numres = len( self.config["sequence"] )

        # Read in data information
        self.config["data_info"] = []
        if data_file != None:
            self.config["backbone_burial_info"] = []
            lines = open( data_file ).readlines()
            for line in lines:
                if len( line ) > 6 and line[:6]=='EXPOSE':
                    cols = string.split( line )
                    for i in range( len(cols)/3 ):
                        self.config["data_info"].append( [int( cols[ 3*i+1] ), cols[3*i+2],cols[3*i+3]] )
                if len( line ) > 15 and line[:15]=='BACKBONE_BURIAL':
                    cols = string.split( line )
                    for i in range( 1,len(cols) ):
                        self.config["backbone_burial_info"].append( int( cols[i] ) )

        pyrimidines = ['c','u']
        purines = ['a','g']
        cst_info = []
        if cst_file != None:
            lines = open( cst_file ).readlines()
            for line in lines:
                if len( line ) > 6 and line[0]!='[':
                    cols = string.split( line )
                    atom_name1 = cols[0]
                    res1 = int( cols[1] )
                    atom_name2 = cols[2]
                    res2 = int( cols[3] )
                    if ( self.config["sequence"][ res1 - 1 ] in pyrimidines and atom_name1=='N1'):
                        atom_name1 = 'N3'
                        print 'correcting atom name for ', res1
                    if ( self.config["sequence"][ res2 - 1 ] in pyrimidines and atom_name2=='N1'):
                        atom_name2 = 'N3'
                        print 'correcting atom name for ', res2
                    if ( self.config["sequence"][ res1 - 1 ] in purines and atom_name1=='N3'):
                        atom_name1 = 'N1'
                        print 'correcting atom name for ', res1
                    if ( self.config["sequence"][ res2 - 1 ] in purines and atom_name2=='N3'):
                        atom_name2 = 'N1'
                        print 'correcting atom name for ', res2

                    cst_info.append( [ atom_name1, res1, atom_name2, res2, string.join( cols[4:] )] )

        pair_map = {}
        all_pairs = []

        complement = {'a':['u'], 'u':['a','g'], 'c':['g'], 'g':['c','u']};

        # Parse out stems
        lines = open( params_file ).readlines()
        cutpoints_original = []
        for line in lines:
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
                self.config["secstruc"] = line[:-1]
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
                        assert ( self.config["sequence"][res1] in complement[ self.config["sequence"][res2] ] )
                assert( len (left_brackets) == 0 )
            else:
                try:
                    cols = string.split( line[:-1] )
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
            tag = 'stems_and_motifs/stem%d.fasta' % (i+1)
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
                          "stems_and_motifs/stem%d_" % (i+1)]
                self.executeCommand(command)
                native_pdb_file_subset =  'stems_and_motifs/stem%d_%s' % (i+1, native_pdb_file )
                print 'Created: ', native_pdb_file_subset


        # Output motif jobs
        for i in range( len(self.config["motifs"]) ):

            # Fasta
            motif_fasta_file = 'stems_and_motifs/motif%d.fasta' % (i+1)
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

            motif_params_file = 'stems_and_motifs/motif%d.params' % (i+1)
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
                command += ["stems_and_motifs/motif%d_" % (i+1)]
                self.executeCommand(command)
                native_pdb_file_subset =  'stems_and_motifs/motif%d_%s' % (i+1, native_pdb_file )
                print 'Created: ', native_pdb_file_subset

            if data_file != None:
                motif_data_file = 'stems_and_motifs/motif%d.data' % ( i+1 )
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

            if cst_file != None:
                motif_cst_file = 'stems_and_motifs/motif%d.cst' % ( i+1 )
                fid_cst = open( motif_cst_file, 'w' )
                fid_cst.write( '[ atompairs ]\n' )
                for cst in cst_info:
                    if cst[1]-1 in motif_res_map.keys() and cst[3]-1 in motif_res_map.keys():
                        fid_cst.write( '%s %d %s %d %s\n' % (cst[0], motif_res_map[cst[1]-1]+1,cst[2],motif_res_map[cst[3]-1]+1,cst[4]) )
                fid_cst.close()
                print 'Created: ', motif_cst_file


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


        params_file = "stems_and_motifs/sequence.params"
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

        ########
        assemble_cst_file = "constraints/default.cst"
        fid = open( assemble_cst_file,'w')
        fid.write('[ atompairs ]\n')
        for cutpoint in cutpoints:
            fid.write( "O3'  %d  P     %d  HARMONIC  1.619  2.0\n" % \
                           ( cutpoint+1, cutpoint+2 ) )
        if cst_file != None:
            for cst in cst_info:
                    fid.write( '%s %d %s %d %s \n' % (cst[0], cst[1], cst[2], cst[3], cst[4]) )

        fid.close()
        print 'Created: ', assemble_cst_file


        #########

    def create_helices(self, dry_run=False):
        self.checkConfig()
        for i in range(len(self.config["stems"])):
            command = ["rna_helix",
                       "-fasta", "stems_and_motifs/stem%d.fasta" % (i+1),
                       "-out:file:silent","stems_and_motifs/stem%d.out" % (i+1)]
            self.executeCommand(command, add_suffix="rosetta", dry_run=dry_run)

    def create_motifs(self, nstruct=20000, cycles=50000, dry_run=False, seed=-1, use_native_information=False):
        self.checkConfig()
        print "Assembly configuration:"
        print "    cycles: %s" % (cycles)
        print "    nstruct: %s" % (nstruct)
        print "    dry_run: %s" % (dry_run)
        print "    random_seed: %s" % (seed)

        for i in range(len(self.config["motifs"])):
            command = ["rna_denovo",
                       "-fasta", "stems_and_motifs/motif%d.fasta" % (i+1),
                       "-params_file", "stems_and_motifs/motif%d.params" % (i+1),
                       "-nstruct", "%d" % (nstruct),
                       "-out:file:silent", "stems_and_motifs/motif%d.out" % (i+1),
                       "-cycles", "%d" % (cycles),
                       "-mute", "all",
                       "-close_loops",
                       "-close_loops_after_each_move",
                       "-minimize_rna"]

            if self.config["native_pdb_file"] != None and use_native_information:
                command += ["-native", "stems_and_motifs/motif%d_%s" % (i+1, self.config["native_pdb_file"])]
            if self.config["data_file"] != None:
                command += ["-data_file", "stems_and_motifs/motif%d.data" % (i+1)]
            if self.config["cst_file"] != None:
                command += ["-cst_file", "stems_and_motifs/motif%d.cst" % (i+1)]
            if self.config["torsions_file"] != None:
                command += ["-vall_torsions", "stems_and_motifs/motif%d.torsions" % (i+1)]

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
                command += ["stems_and_motifs/stem%d.out" % (n+1)]

            command += ["-chunk_res"]
            command += self.make_tag_with_dashes(stem_chunk_res)

            if seed != -1:
                command += ["-constant_seed", "-jran", "%d" % (seed)]

            self.executeCommand(command, add_suffix="rosetta", dry_run=dry_run)

    def assemble(self, nstruct=20000, cycles=50000, constraints_file="constraints/default.cst", dry_run=False, seed=-1, use_native_information=False):
        self.checkConfig()
        # In case the seed isn't specify, we do what rosetta does to get a random number, and seed it with that
        # this way we know the number beforehand and can choose an appropriate output filename.
        if seed == -1:
            seed = struct.unpack("=i", os.urandom(4))[0]

        print "Assembly configuration:"
        print "    cycles: %s" % (cycles)
        print "    nstruct: %s" % (nstruct)
        print "    constraints: %s" % (constraints_file)
        print "    dry_run: %s" % (dry_run)
        print "    random_seed: %s" % (seed)
        self.checkFileExistence(constraints_file)

        command = ["rna_denovo",
                   "-minimize_rna",
                   "-fasta", self.config["fasta_file"],
                   "-in:file:silent_struct_type", "binary_rna",
                   "-cycles", "%d" % (cycles),
                   "-nstruct", "%d" % (nstruct),
                   "-out:file:silent", "assembly/%s_%d.out" % (splitext(basename(constraints_file))[0], seed),
                   "-params_file", "stems_and_motifs/sequence.params",
                   "-cst_file", constraints_file,
                   "-close_loops",
                   "-in:file:silent"]

        for i in range(len(self.config["stems"])):
            command += ["stems_and_motifs/stem%d.out" % (i+1)]

        for i in range(len(self.config["motifs"])):
            command += ["stems_and_motifs/motif%d.out" % (i+1)]

        chunk_res = []
        for n in range( len(self.config["stems"])  ):
            stem = self.config["stems"][n]
            for q in range( len( stem ) ):        chunk_res.append(stem[q][0] + 1)
            for q in range( len( stem ) ):        chunk_res.append(stem[ len(stem) - q - 1][1] + 1)

        for n in range( len(self.config["motifs"]) ):
            motif_res = self.config["motifs"][n]
            for m in motif_res: chunk_res.append(m+1)

        command += ["-chunk_res"]
        command += self.make_tag_with_dashes(chunk_res)

        if self.config["native_pdb_file"] != None and use_native_information:
            command += ["-native", self.config["native_pdb_file"]]


        if self.config["torsions_file"] != None:
            command += ["-vall_torsions", self.config["torsions_file"]]
        if self.config["data_file"] != None:
            command += ["-data_file", self.config["data_file"]]

        command += ["-constant_seed", "-jran", "%d" % (seed)]

        self.executeCommand(command, add_suffix="rosetta", dry_run=dry_run)


    def deleteGlob(self, pattern, print_notice=True):
        for f in glob.glob(pattern):
            if print_notice:
                print "deleting %s..." % (f)
            try:
                os.remove(f)
            except:
                pass

    def extractPOnly(self, filename):
        p_only = ""
        with open(filename, "r") as fd:
            for line in fd:
                fields = line.split()
                if fields[0] == "ATOM" and "P" in fields[2]:
                    p_only += line
            p_only += "TER\n"
            return p_only

    def writePdb(self, filename, data, model=1, remark=None, append=False):
        with open(filename, "a" if append else "w") as output_pdb_fd:
            output_pdb_fd.write("MODEL %d\n" % (model))
            output_pdb_fd.write("REMARK %s\n" % (remark))
            output_pdb_fd.write(data)
            output_pdb_fd.write("ENDMDL\n")

    def extract(self, constraints_file="constraints/default.cst"):
        self.checkConfig()
        cst_name = splitext(basename(constraints_file))[0]
        print "Extraction:"
        print "    constraints: %s" % (cst_name)
        self.checkFileExistence(constraints_file)

        # create dict if needed
        if not self.config.has_key("evaluate"):
            self.config["evaluate"] = {}

        # delete old files
        self.deleteGlob("assembly/%s.pdb" % (cst_name))
        self.deleteGlob("output/%s*.pdb" % (cst_name))
        self.deleteGlob("output/%s.log" % (cst_name))

        # clear old evaluation data
        self.config["evaluate"][cst_name] = {"models": {}, "clusters": {}}

        # model counter
        model = 0

        # cleanup
        self.deleteGlob("temp/%s*.pdb" % (cst_name))

        # loop over all out files matching the constraint
        for f in sorted(glob.glob("assembly/%s_*.out" % (cst_name)), key=natural_sort_key):
            print "  processing rosetta silent file: %s..." % (f)

            description = splitext(basename(f))[0]

            # read out file and store a dict of the scores and the full score line
            scores = dict()
            score_regex = re.compile("^SCORE:\s+([-0-9.]+)\s.*(S_\d+)$")
            with open(f, "r") as fd:
                for line in fd:
                    m = score_regex.match(line)
                    if m:
                        scores[m.group(2)]={"line": line.rstrip(), "score": float(m.group(1))}


            # delete any pdb files if existing
            self.deleteGlob("S_*.pdb")

            # extract pdb files
            command = ["rna_extract", "-in:file:silent", f, "-in:file:silent_struct_type", "rna"]
            self.executeCommand(command, add_suffix="rosetta")

            # loop over all extracted pdb files
            for fe in sorted(glob.glob("S_*.pdb")):
                model += 1
                name = fe[:-4]
                remark = "%s %s" % (description, scores[name]["line"])

                # extract all P atoms from the pdb file
                p_only = self.extractPOnly(fe)

                # append p only model to trajectory
                self.writePdb("assembly/%s.pdb" % (cst_name), data=p_only, model=model, remark=remark, append=True)

                # write p only model to temp file
                self.writePdb("temp/%s_%09d_p.pdb" % (cst_name, model), data=p_only, model=model, remark=remark)

                # move original pdb to temp directory
                shutil.move(fe, "temp/%s_%09d.pdb" % (cst_name, model))

                # create evaluation dict for model
                self.config["evaluate"][cst_name]["models"][model] = {"source_file": f, "source_name": name, "score": scores[name]["score"]}


    # TODO: set the cutoff value back to 0.40? It was set to 0.41 because the old bash script used integer comparison and even 0.409 was treated as 0.40
    def evaluate(self, constraints_file="constraints/default.cst", cluster_limit=10, cluster_cutoff=0.41):
        self.checkConfig()
        cst_name = splitext(basename(constraints_file))[0]
        print "Evaluation configuration:"
        print "    cluster_limit: %s" % (cluster_limit)
        print "    cluster_cutoff: %s" % (cluster_cutoff)
        print "    constraints: %s" % (cst_name)
        self.checkFileExistence(constraints_file)

        # check if extraction step was run before
        try:
            self.config["evaluate"][cst_name]["models"][1]["score"]
        except:
            raise SimulationException("No extracted data for constraint '%s' found. Did you run with --extract?\n" % (cst_name))


        # clear old evaluation clusters
        self.config["evaluate"][cst_name]["clusters"] = {}
        for m in self.config["evaluate"][cst_name]["models"].iteritems():
            m[1].pop("native_rmsd", None)
            m[1].pop("cluster", None)
            m[1].pop("rmsd_to_primary", None)


        # calculate native rmsd values for all models if native pdb available
        if self.config["native_pdb_file"] != None:
            # create native_p.pdb if needed
            native_p_only = "%s_p.pdb" % (self.config["native_pdb_file"][:-4])
            if not os.path.exists(native_p_only):
                print "  creating p only native pdb file: %s" % (native_p_only)
                self.writePdb(native_p_only, data=self.extractPOnly(self.config["native_pdb_file"]), model=1, remark="native p only")
            print "  caluculating rmsd values to native structure for all models..."
            sys.stdout.write("    ")
            self.executeCommand(["g_rms", "-quiet", "-s", "native_p.pdb", "-f", "assembly/%s.pdb" % (cst_name)], add_suffix="gromacs", stdin="1\n1\n", quiet=True)
            with open("rmsd.xvg", "r") as r:
                for line in r:
                    if not re.match(r"^[\s\d-]", line):
                        continue
                    model, native_rmsd = line.split()
                    self.config["evaluate"][cst_name]["models"][int(float(model))]["native_rmsd"] = float(native_rmsd)
            self.deleteGlob("rmsd.xvg", print_notice=False)

        # cluster counter
        cluster = 0

        log = open("output/%s.log" % (cst_name), "w")
        print "  sorting models into clusters..."
        # loop over all models sorted by their score
        for model, data in sorted(self.config["evaluate"][cst_name]["models"].items(), key=lambda x: x[1]["score"]):
            filename_pdb = "temp/%s_%09d.pdb" % (cst_name, model)
            filename_pdb_p = "temp/%s_%09d_p.pdb" % (cst_name, model)

            # check if the current structure matches a cluster
            matches_cluster = 0
            rmsd_to_cluster_primary = cluster_cutoff
            for c in range(cluster):
                # calculate rmsd between cluster and pdb
                self.executeCommand(["g_rms", "-quiet", "-s", "output/%s_%d_p.pdb" % (cst_name, c + 1), "-f", filename_pdb_p], add_suffix="gromacs", stdin="1\n1\n", quiet=True, print_commands=False)
                with open("rmsd.xvg", "r") as r:
                    for line in r:
                        pass
                    new_rmsd = float(line.split()[1])
                    if new_rmsd < rmsd_to_cluster_primary:
                        rmsd_to_cluster_primary = new_rmsd
                        matches_cluster = c + 1
                self.deleteGlob("rmsd.xvg", print_notice=False)

            if matches_cluster == 0:
                cluster = cluster + 1
                shutil.copyfile(filename_pdb, "output/%s_%d.pdb" % (cst_name, cluster))
                shutil.copyfile(filename_pdb_p, "output/%s_%d_p.pdb" % (cst_name, cluster))
                self.config["evaluate"][cst_name]["clusters"][cluster] = {"primary_model": model}
                self.config["evaluate"][cst_name]["models"][model]["cluster"] = cluster
                self.config["evaluate"][cst_name]["models"][model]["rmsd_to_primary"] = 0
                output = "new cluster: %02s, model: %05s, native_rmsd: %.7f, score: %.3f" % (cluster,
                                                                                           model,
                                                                                           self.config["evaluate"][cst_name]["models"][model]["native_rmsd"],
                                                                                           self.config["evaluate"][cst_name]["models"][model]["score"])
                print "    %s" % (output)
                log.write(output + "\n")
            else:
                self.config["evaluate"][cst_name]["models"][model]["cluster"] = matches_cluster
                self.config["evaluate"][cst_name]["models"][model]["rmsd_to_primary"] = rmsd_to_cluster_primary
                output = "    cluster: %02s, model: %05s, native_rmsd: %.7f, score: %.3f, rmsd_to_cluster_primary: %.7f" % (matches_cluster,
                                                                                                                        model,
                                                                                                                        self.config["evaluate"][cst_name]["models"][model]["native_rmsd"],
                                                                                                                        self.config["evaluate"][cst_name]["models"][model]["score"],
                                                                                                                        self.config["evaluate"][cst_name]["models"][model]["rmsd_to_primary"])
                print "    %s" % (output)
                log.write(output + "\n")

            if cluster == cluster_limit:
                print "    maximum number of clusters reached. stopping..."
                break

    # TODO: change output formatting to something else?
    def compare(self):
        self.checkConfig()
        sys.stdout.write(self.config["name"])
        # loop over all different constraint sets
        for cst_file in sorted(glob.glob("constraints/*.cst")):
            cst_name = splitext(basename(cst_file))[0]
            sys.stdout.write("   " + cst_name)

            for c in (1, 5, 10):
                min_rmsd = 999
                for c2 in range(1, c + 1):
                    model = self.config["evaluate"][cst_name]["clusters"][c2]["primary_model"]
                    rmsd = self.config["evaluate"][cst_name]["models"][model]["native_rmsd"]
                    if rmsd < min_rmsd:
                        min_rmsd = rmsd
                sys.stdout.write(" %.2f" % (min_rmsd * 10))
        sys.stdout.write("\n")

