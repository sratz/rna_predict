#!/usr/local/bin/python2.7
# encoding: utf-8
'''
rna_prediction.main -- predict tertiary rna structure

rna_prediction.main is a tool to predict tertiary rna structure based on secondary structure information and a set of constraints

It defines classes_and_methods

@author:     sra

@copyright:  2014 SCC/MBS @ KIT. All rights reserved.

@license:    license

@contact:    sebastian.ratz@student.kit.edu
'''

import sys
import os

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

from rna_prediction.simulation import RNAPrediction
from rna_prediction.simulation import SysConfig
from rna_prediction.simulation import SimulationException

__all__ = []
__version__ = 0.1
__date__ = '2014-08-13'
__updated__ = '2014-08-13'

DEBUG = 1
TESTRUN = 0
PROFILE = 0


def main(argv=None): # IGNORE:C0111

    def printToStderr(message, prefix="error"):
        sys.stderr.write("%s: %s: %s\n" % (program_name, prefix, message))

    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __doc__.split("\n")[1]
    program_license = '''%s

  Created by sra on %s.
  Copyright 2014 organization_name. All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    # Setup argument parser
    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("--prepare", dest="prepare", action="store_true", help="prepare stems and motifs [default: %(default)s]")
    parser.add_argument("--create-helices", dest="create_helices", action="store_true", help="create ideal a-helices [default: %(default)s]")
    parser.add_argument("--create-motifs", dest="create_motifs", action="store_true", help="create motifs [default: %(default)s]")
    parser.add_argument("--assemble", dest="assemble", action="store_true", help="assemble [default: %(default)s]")
    parser.add_argument("--extract", dest="extract", action="store_true", help="extract pdb data and scrore [default: %(default)s]")
    parser.add_argument("--evaluate", dest="evaluate", action="store_true", help="evaluate data (clusters) [default: %(default)s]")
    parser.add_argument("--compare", dest="compare", action="store_true", help="print comparison of prediction to native structure [default: %(default)s]")
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument('-q', "--quiet", dest="quiet", action="store_true", help="don't print config on start [default: %(default)s]")
    parser.add_argument("--config", dest="config", help="modify config variable", metavar=("KEY", "VALUE"), nargs=2)
    parser.add_argument(dest="basepaths", help="paths to base directories [default: %(default)s]", metavar="basepaths", nargs="*", default=".")
    parser.add_argument("--name", dest="name", help="simulation name (for example the analyzed pdb). basename of directory will be used if not given [default: %(default)s]")
    parser.add_argument("--native", dest="native", help="native pdb [default: %(default)s]")
    parser.add_argument("--cycles", dest="cycles", help="number of cycles [default: %(default)s]", default=20000, type=int)
    parser.add_argument("--nstruct", dest="nstruct", help="number of structures to create [default: %(default)s]", default=50000, type=int)
    parser.add_argument("--seed", dest="seed", help="random seed [default: %(default)s]", default=-1, type=int)
    parser.add_argument("--use-native", dest="use_native", action="store_true", help="use native information for motif generation and assembly [default: %(default)s]")
    parser.add_argument("--cluster-cutoff", dest="cluster_cutoff", help="cluster cutoff in nm [default: %(default)s]", default=0.41, type=float)
    parser.add_argument("--cluster-limit", dest="cluster_limit", help="maximum number of clusters to create [default: %(default)s]", default=10, type=int)
    parser.add_argument("-n", "--dry-run", dest="dry_run", action="store_true", help="don't execute anything, just print commands [default: %(default)s]")

    # system configuration
    sysconfig = SysConfig()
    if sysconfig.isSysConfigOk():
        sysconfig.printSysConfig()
    else:
        printToStderr("Cannot access rosetta binaries. Check sconfig file!")
        return 1

    # Process arguments
    args = parser.parse_args()

    # TODO: Is it really useful to allow more than one path argument? This makes error / exit code handling complicated.
    # TODO: For now, at least try to make this work reliably if only one path is given.
    exit_code = 0
    for path in args.basepaths:
        try:
            os.chdir(path)
        except:
            printToStderr("Invalid basepath: " + path)
            continue

        try:
            p = RNAPrediction(sysconfig)

            if args.config:
                p.modifyConfig(args.config[0], args.config[1])
                p.saveConfig()

            name = args.name
            if name is None:
                name = os.path.basename(os.path.abspath(path))

            if not args.quiet:
                p.printConfig()

            if args.prepare:
                p.prepare(native_pdb_file=args.native, name=name)
                p.saveConfig()
            if args.create_helices:
                p.create_helices(dry_run=args.dry_run)
                p.saveConfig()
            if args.create_motifs:
                p.create_motifs(dry_run=args.dry_run, nstruct=args.nstruct, cycles=args.cycles, seed=args.seed, use_native_information=args.use_native)
                p.saveConfig()
            if args.assemble:
                p.assemble(dry_run=args.dry_run, nstruct=args.nstruct, cycles=args.cycles, seed=args.seed, use_native_information=args.use_native)
                p.saveConfig()
            if args.extract:
                p.extract()
                p.saveConfig()
            if args.evaluate:
                p.evaluate(cluster_limit=args.cluster_limit, cluster_cutoff=args.cluster_cutoff)
                p.saveConfig()
            if args.compare:
                p.compare()
        except SimulationException, e:
            printToStderr(e)
            exit_code = 1
            continue
        except KeyboardInterrupt:
            exit_code = 1
            continue
        except Exception, e:
            if DEBUG or TESTRUN:
                raise
            else:
                printToStderr(e)
            # TODO: Save the config here, so we don't lose anything even in case of errors?
            exit_code = 1
            continue


    return exit_code

if __name__ == "__main__":
    sys.exit(main())
