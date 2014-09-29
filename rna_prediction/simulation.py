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
import dcatools
from os.path import splitext, basename, abspath


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]


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
        self.sysconfig = sysconfig
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
        if name is None:
            name = os.path.basename(os.path.abspath(os.getcwd()))
        for f in [fasta_file, params_file, native_pdb_file, data_file, cst_file, torsions_file]:
            if f is not None:
                self.checkFileExistence(f)
        self.makeDirectory("stems_and_motifs")
        self.makeDirectory("assembly")
        self.makeDirectory("constraints")
        self.makeDirectory("output")
        self.makeDirectory("temp")
        self.makeDirectory("dca")
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

    def create_helices(self, dry_run=False, threads=1):
        self.checkConfig()
        commands = list()
        for i in range(len(self.config["stems"])):
            command = ["rna_helix",
                       "-fasta", "stems_and_motifs/stem%d.fasta" % (i+1),
                       "-out:file:silent","stems_and_motifs/stem%d.out" % (i+1)]
            commands.append(Command(command, add_suffix="rosetta", dry_run=dry_run))
        self.executeCommands(commands, threads=threads)

        # rna_helix dumps <sequence>.pdb files in the working directory.
        # These can be useful for us, so move them to the stems_and_motifs directory with a proper stemX.pdb filename.
        # TODO: What happens if there are multiple helices with the same sequence? In that case we need to move them
        #       out of the way before running the next command.
        #       Just disable multithreading here (not really needed with stem generation anyways)?
        for i in range(len(self.config["stems"])):
            with open("stems_and_motifs/stem%d.fasta" % (i + 1), "r") as f:
                l = f.readlines()
                sequence = l[1].strip()
                shutil.move("%s.pdb" % (sequence), "stems_and_motifs/stem%d.pdb" % (i + 1))


    def mergeMotifs(self, target, sources):
        n = 0
        pattern_header = re.compile("^(?:SEQUENCE:|SCORE:\s+score).*")
        pattern_normal = re.compile("^(.*)S_\d+$")
        with open(target, "a+") as t:
            for line in t:
                pass
            if 'line' in locals():
                m = re.match(".*S_(\d+)$", line)
                if m:
                    n = int(m.group(1))
            for source in glob.glob(sources):
                print "merging %s into %s" % (source, target)
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
        self.deleteGlob(sources, print_notice=True)
        return n


    def create_motifs(self, nstruct=50000, cycles=20000, dry_run=False, seed=None, use_native_information=False, threads=1):
        self.checkConfig()
        print "Assembly configuration:"
        print "    cycles: %s" % (cycles)
        print "    nstruct: %s" % (nstruct)
        print "    dry_run: %s" % (dry_run)
        print "    random_seed: %s" % (seed)
        print "    threads: %s" % (threads)


        # merge all motifs abd check what we have so far
        n_motifs = len(self.config["motifs"])
        completed = {}
        for i in range(n_motifs):
            completed[i] = self.mergeMotifs("stems_and_motifs/motif%d.out" % (i + 1), "stems_and_motifs/motif%d_*.out" % (i + 1))

        commands = list()
        for i in range(len(self.config["motifs"])):
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
                       "-fasta", "stems_and_motifs/motif%d.fasta" % (i+1),
                       "-params_file", "stems_and_motifs/motif%d.params" % (i+1),
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

            command += ["-in:file:input_res"]
            command += self.make_tag_with_dashes(stem_chunk_res)

            for j in range(threads):
                if structs_threads[j] == 0:
                    continue
                command_full = command + ["-out:file:silent", "stems_and_motifs/motif%d_%d.out" % (i + 1, j + 1),
                                          "-nstruct", "%d" % (structs_threads[j])]
                if seed != None:
                    command_full += ["-constant_seed", "-jran", "%d" % (seed)]
                    seed += 1

                commands.append(Command(command_full, add_suffix="rosetta", dry_run=dry_run))
        self.executeCommands(commands, threads=threads)

        # merge motifs
        for i in range(n_motifs):
            self.mergeMotifs("stems_and_motifs/motif%d.out" % (i + 1), "stems_and_motifs/motif%d_*.out" % (i + 1))


    # TODO: When documenting later, explain that with assemble, nstruct is used for each single thread, while with createMotifs, it is distributed.
    def assemble(self, nstruct=50000, cycles=20000, constraints_file="constraints/default.cst", dry_run=False, seed=None, use_native_information=False, threads=1):
        self.checkConfig()
        print "Assembly configuration:"
        print "    cycles: %s" % (cycles)
        print "    nstruct: %s" % (nstruct)
        print "    constraints: %s" % (constraints_file)
        print "    dry_run: %s" % (dry_run)
        print "    random_seed: %s" % (seed)
        self.checkFileExistence(constraints_file)

        commands = list()

        command = ["rna_denovo",
                   "-minimize_rna",
                   "-fasta", self.config["fasta_file"],
                   "-in:file:silent_struct_type", "binary_rna",
                   "-cycles", "%d" % (cycles),
                   "-nstruct", "%d" % (nstruct),
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

        command += ["-in:file:input_res"]
        command += self.make_tag_with_dashes(chunk_res)

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

            command_full = command + ["-out:file:silent", "assembly/%s_%d.out" % (splitext(basename(constraints_file))[0], seed_used),
                                      "-constant_seed", "-jran", "%d" % (seed_used)]

            commands.append(Command(command_full, add_suffix="rosetta", dry_run=dry_run))

        self.executeCommands(commands, threads=threads)


    def deleteGlob(self, pattern, print_notice=True):
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

    def extractPOnly(self, filename):
        p_only = ""
        with open(filename, "r") as fd:
            # only extract first chain
            chain_id = None
            for line in fd:
                fields = line.split()
                if fields[0] == "ATOM" and "P" in fields[2]:
                    if chain_id is None:
                        chain_id = fields[4]
                    elif chain_id != fields[4]:
                        break
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

        # delete old files
        self.deleteGlob("assembly/%s.pdb" % (cst_name))
        self.deleteGlob("output/%s*.pdb" % (cst_name))
        self.deleteGlob("output/%s.log" % (cst_name))
        self.deleteGlob("output/%s.dat" % (cst_name))

        # create dict to store evaluation data
        evalData = {"models": {}, "clusters": {}}

        # model counter
        model = 0

        # cleanup
        tmp_dir = "temp/%s" % (cst_name)
        self.deleteGlob(tmp_dir)
        self.makeDirectory(tmp_dir)

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
            command = ["rna_extract", "-in:file:silent", abspath(f), "-in:file:silent_struct_type", "rna"]
            owd = os.getcwd()
            try:
                os.chdir(tmp_dir)
                self.executeCommand(command, add_suffix="rosetta")
            finally:
                os.chdir(owd)

            # loop over all extracted pdb files
            for fe in sorted(glob.glob("%s/S_*.pdb" % (tmp_dir))):
                model += 1
                name = basename(fe)[:-4]
                remark = "%s %s" % (description, scores[name]["line"])

                # extract all P atoms from the pdb file
                p_only = self.extractPOnly(fe)

                # append p only model to trajectory
                self.writePdb("assembly/%s.pdb" % (cst_name), data=p_only, model=model, remark=remark, append=True)

                # write p only model to temp file
                self.writePdb("%s/%09d_p.pdb" % (tmp_dir, model), data=p_only, model=model, remark=remark)

                # move original pdb to temp directory
                shutil.move(fe, "%s/%09d.pdb" % (tmp_dir, model))

                # create evaluation dict for model
                evalData["models"][model] = {"source_file": f, "source_name": name, "score": scores[name]["score"]}

        # save evaluation data to file
        with open("output/%s.dat" % (cst_name), "w") as f:
            pickle.dump(evalData, f)


    # TODO: set the cutoff value back to 0.40? It was set to 0.41 because the old bash script used integer comparison and even 0.409 was treated as 0.40
    def evaluate(self, constraints_file="constraints/default.cst", cluster_limit=10, cluster_cutoff=0.41):
        self.checkConfig()
        cst_name = splitext(basename(constraints_file))[0]
        print "Evaluation configuration:"
        print "    cluster_limit: %s" % (cluster_limit)
        print "    cluster_cutoff: %s" % (cluster_cutoff)
        print "    constraints: %s" % (cst_name)
        self.checkFileExistence(constraints_file)

        # load evaluation data and check if extraction step was run before
        evalDataFilename = "output/%s.dat" % (cst_name)
        try:
            with open(evalDataFilename, "r") as f:
                evalData = pickle.load(f)
            evalData["models"][1]["score"]
        except KeyError:
            raise SimulationException("Broken data file %s. Try running --extract again" % (evalDataFilename))
        except:
            raise SimulationException("Data file %s for constraint not found. Did you forget to run --extract?" % (evalDataFilename))

        tmp_dir = "temp/%s" % (cst_name)
        if not os.path.isdir(tmp_dir) or not os.path.isfile("%s/%09d_p.pdb" % (tmp_dir, len(evalData["models"]))):
            raise SimulationException("No extracted pdb for constraint '%s' found. Did you delete temp/ files?" % (cst_name))

        # clear old evaluation clusters
        self.deleteGlob("output/%s*.pdb" % (cst_name))
        self.deleteGlob("output/%s.log" % (cst_name))
        evalData["clusters"] = {}
        for m in evalData["models"].iteritems():
            m[1].pop("native_rmsd", None)
            m[1].pop("cluster", None)
            m[1].pop("rmsd_to_primary", None)


        # calculate native rmsd values for all models if native pdb available
        filename_rmsd = "%s/rmsd.xvg" % (tmp_dir)
        if self.config["native_pdb_file"] != None:
            # create native_p.pdb if needed
            native_p_only = "%s_p.pdb" % (self.config["native_pdb_file"][:-4])
            if not os.path.exists(native_p_only):
                print "  creating p only native pdb file: %s" % (native_p_only)
                self.writePdb(native_p_only, data=self.extractPOnly(self.config["native_pdb_file"]), model=1, remark="native p only")
            print "  caluculating rmsd values to native structure for all models..."
            sys.stdout.write("    ")
            self.executeCommand(["g_rms", "-quiet", "-s", "native_p.pdb", "-f", "assembly/%s.pdb" % (cst_name), "-o", filename_rmsd], add_suffix="gromacs", stdin="1\n1\n", quiet=True)
            with open(filename_rmsd, "r") as r:
                for line in r:
                    if not re.match(r"^[\s\d-]", line):
                        continue
                    model, native_rmsd = line.split()
                    evalData["models"][int(float(model))]["native_rmsd"] = float(native_rmsd)
            self.deleteGlob(filename_rmsd, print_notice=False)

        # cluster counter
        cluster = 0

        log = open("output/%s.log" % (cst_name), "w")
        print "  sorting models into clusters..."

        # do we have native rmsd information?
        use_native = self.config["native_pdb_file"] != None
        # loop over all models sorted by their score
        for model, data in sorted(evalData["models"].items(), key=lambda x: x[1]["score"]):
            filename_pdb = "%s/%09d.pdb" % (tmp_dir, model)
            filename_pdb_p = "%s/%09d_p.pdb" % (tmp_dir, model)

            # check if the current structure matches a cluster
            matches_cluster = 0
            rmsd_to_cluster_primary = cluster_cutoff
            for c in range(cluster):
                # calculate rmsd between cluster and pdb
                self.executeCommand(["g_rms", "-quiet", "-s", "output/%s_%d_p.pdb" % (cst_name, c + 1), "-f", filename_pdb_p, "-o", filename_rmsd], add_suffix="gromacs", stdin="1\n1\n", quiet=True, print_commands=False)
                with open(filename_rmsd, "r") as r:
                    for line in r:
                        pass
                    new_rmsd = float(line.split()[1])
                    if new_rmsd < rmsd_to_cluster_primary:
                        rmsd_to_cluster_primary = new_rmsd
                        matches_cluster = c + 1
                self.deleteGlob(filename_rmsd, print_notice=False)

            if matches_cluster == 0:
                cluster = cluster + 1
                shutil.copyfile(filename_pdb, "output/%s_%d.pdb" % (cst_name, cluster))
                shutil.copyfile(filename_pdb_p, "output/%s_%d_p.pdb" % (cst_name, cluster))
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
        with open(evalDataFilename, "w") as f:
            pickle.dump(evalData, f)


    def compare(self):
        self.checkConfig()
        if self.config["native_pdb_file"] is None:
            raise SimulationException("Cannot compare without native information.")
        print self.config["name"]

        def printComparisonLine(cst_name, comparisons):
            print "  %-035s %05s %05s %05s" % (cst_name, comparisons[0], comparisons[1], comparisons[2])

        # loop over all different constraint sets
        for cst_file in sorted(glob.glob("constraints/*.cst"), key=natural_sort_key):
            cst_name = splitext(basename(cst_file))[0]
            try:
                with open("output/%s.dat" % (cst_name), "r") as f:
                    evalData = pickle.load(f)
                    evalData["models"][1]["native_rmsd"]
            except:
                printComparisonLine(cst_name, ["-", "-", "-"])
                continue
            comparisons = []
            for c in (1, 5, 10):
                min_rmsd = 999
                for c2 in range(1, c + 1):
                    model = evalData["clusters"][c2]["primary_model"]
                    rmsd = evalData["models"][model]["native_rmsd"]
                    if rmsd < min_rmsd:
                        min_rmsd = rmsd
                comparisons.append("%.2f" % (min_rmsd * 10))
            printComparisonLine(cst_name, comparisons)


    def makeConstraints(self, pdbMapping=None, dcaPredictionFileName="dca/dca.txt", outputFileName=None, numberDcaPredictions=100, cstFunction="FADE -100 26 20 -2 2"):
        # TODO: make this dependable on the prepare step? Or separate the whole constraints creation into an independent application?
        self.checkConfig()
        if outputFileName is None:
            outputFileName = "constraints/%s.cst" % (splitext(basename(dcaPredictionFileName))[0])
        print "Constraints creation:"
        print "    dcaPredictionFileName: %s" % (dcaPredictionFileName)
        print "    outputFileName: %s" % (outputFileName)
        print "    numberDcaPredictions: %d" % (numberDcaPredictions)
        print "    function: %s" % (cstFunction)
        print "  %s pdbMapping (user): %s" % (" " if pdbMapping is None else "*", pdbMapping)
        self.checkFileExistence(dcaPredictionFileName)

        if pdbMapping is not None:
            pdbMapping = dcatools.createPdbMapping(self.config["sequence"], pdbMapping)
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
                            print "  %s pdbMapping (file): %s" % (" " if pdbMapping is not None else "*", m.group(2))
                            if pdbMapping is not None:
                                continue
                            pdbMapping = dcatools.createPdbMapping(self.config["sequence"], m.group(2))
                    continue
                # data line
                if len(dca) >= numberDcaPredictions:
                    break
                parts = line.split(" ")
                if pdbMapping is not None:
                    try:
                        dca.append([pdbMapping[int(parts[0])], pdbMapping[int(parts[1])]])
                    except:
                        raise SimulationException("Invalid PDB mapping. Could not access residue: %s" % (parts[:2]))
                else:
                    dca.append([int(parts[0]), int(parts[1])])

        atoms = []
        first = True
        for res in self.config["sequence"]:
            res = res.upper()
            atoms.append((res, dcatools.getAtomsForRes(res, termPhosphate=(not first))))
            first = False

        distanceMap = dcatools.getContactDistanceMap()
        distanceMapMean = dcatools.getMeanDistanceMapMean(distanceMap, meanCutoff=6.0, stdCutoff=3.0)

        print "Creating constraints:"
        shutil.copy("constraints/default.cst", outputFileName)
        with open(outputFileName, "a") as out:
            for residueContact in dca:
                res1 = atoms[residueContact[0] - 1]
                res2 = atoms[residueContact[1] - 1]
                contactKey = res1[0] + res2[0]

                for atom1 in res1[1]:
                    for atom2 in res2[1]:
                        atomContactKey = atom1 + '-' + atom2
                        if atomContactKey in distanceMapMean[contactKey]:
                            distance = distanceMapMean[contactKey][atomContactKey][0] / 10.0
                            print "%s %s %s %s" % (residueContact, contactKey, atomContactKey, distance)
                            out.write("%s %s %s %s %s\n" % (atom1, residueContact[0], atom2, residueContact[1], cstFunction))

    def editConstraints(self, inputFileName, outputFileName=None, cstFunction="FADE -100 26 20 -2 2"):
        ':type cstFunction:string'
        cstName = splitext(basename(inputFileName))[0]
        cstFunctionUnderscore = cstFunction.replace(" ", "_")
        if outputFileName is None:
            outputFileName = "constraints/%s_%s.cst" % (cstName, cstFunctionUnderscore)
        else:
            outputFileName = outputFileName.replace("%f", cstFunctionUnderscore)
        print "Constraints editing:"
        print "    inputFileName:  %s" % (inputFileName)
        print "    outputFileName: %s" % (outputFileName)
        print "    function:       %s" % (cstFunction)
        self.checkFileExistence(inputFileName)

        if inputFileName == outputFileName:
            raise SimulationException("Input and output filename cannot be the same")

        pattern = re.compile(r"^(\S+\s+\d+\s+\S+\s+\d+)\s+(\S+)\s+.+$")
        with open(inputFileName, "r") as inputFd:
            with open(outputFileName, "w") as outputFd:
                for line in inputFd:
                    m = pattern.match(line)
                    # TODO: is this a good idea? what if we actually want to use HARMONIC for our custom constraints at some point?
                    if m and m.group(2) != "HARMONIC":
                        outputFd.write("%s %s\n" % (m.group(1), cstFunction))
                    else:
                        outputFd.write(line)
