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

import argparse
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

from rna_prediction.simulation import RNAPrediction
from rna_prediction.simulation import SimulationException
from rna_prediction.sysconfig import SysConfig
from .dcatools import DcaException
from . import dcatools
from . import tools

__all__ = []
__version__ = 0.1
__date__ = '2014-08-13'
__updated__ = '2014-08-13'

DEBUG = 1
TESTRUN = 0
PROFILE = 0


def main(argv=None): # IGNORE:C0111

    # use line-buffered stdout and stderr
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)
    sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 1)

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

    class SubcommandHelpFormatter(RawDescriptionHelpFormatter):
        def _format_action(self, action):
            parts = super(RawDescriptionHelpFormatter, self)._format_action(action)
            if action.nargs == argparse.PARSER:
                parts = "\n".join(parts.split("\n")[1:])
            return parts

    # Setup argument parser
    parser = ArgumentParser(description=program_license, formatter_class=SubcommandHelpFormatter)
    subparsers = parser.add_subparsers(dest="subcommand", title="subcommands", metavar="subcommand")

    parser.add_argument('-j', '--threads', dest="threads", help="maximum number of parallel subprocesses [default: %(default)s]", default=1, type=int)
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument('-q', "--quiet", dest="quiet", action="store_true", help="don't print config on start [default: %(default)s]")
    parser.add_argument("-n", "--dry-run", dest="dry_run", action="store_true", help="don't execute and only print external commands [default: %(default)s]")

    parser_prepare = subparsers.add_parser("prepare", help="prepare stems and motifs")
    parser_preparecst = subparsers.add_parser("prepare-cst", help="prepare constraints file for motif generation and assembly")
    parser_helices= subparsers.add_parser("create-helices", help="create ideal a-helices")
    parser_motifs = subparsers.add_parser("create-motifs", help="create motifs")
    parser_assemble = subparsers.add_parser("assemble", help="assemble models")
    parser_evaluate = subparsers.add_parser("evaluate", help="evaluate data (clustering and rmsd calculation)")
    parser_compare = subparsers.add_parser("compare", help="print comparison of prediction to native structure")
    parser_makeconstraints = subparsers.add_parser("make-constraints", help="create a constraints file from a dca prediction")
    parser_editconstraints = subparsers.add_parser("edit-constraints", help="replace rosetta function in a constraints file")
    parser_status = subparsers.add_parser("status", help="print status summary")
    parser_config = subparsers.add_parser("config", help="modify config variable")
    parser_tools = subparsers.add_parser("tools", help="various plot generation tools")

    parser_prepare.add_argument("--name", dest="name", help="simulation name (for example the analyzed pdb) [default: basename of directory]")
    parser_prepare.add_argument("--native", dest="native", help="native pdb [default: %(default)s]")
    parser_prepare.add_argument("--sequence", dest="sequence", help="sequence fasta file [default: %(default)s]", default="sequence.fasta")
    parser_prepare.add_argument("--secstruct", dest="secstruct", help="secondary structure file [default: %(default)s]", default="secstruct.txt")
    parser_preparecst.add_argument("--override-motifs-cst", dest="motifs_cst", help="use motifs from a different constraints set [default: %(default)s]", const="none", nargs="?")
    parser_motifs.add_argument("--cycles-motifs", dest="cycles_motifs", help="number of cycles for motif generation [default: %(default)s]", default=5000, type=int)
    parser_motifs.add_argument("--nstruct-motifs", dest="nstruct_motifs", help="number of motif structures to create [default: %(default)s]", default=4000, type=int)
    parser_assemble.add_argument("--cycles-assembly", dest="cycles_assembly", help="number of cycles for assembly [default: %(default)s]", default=20000, type=int)
    parser_assemble.add_argument("--nstruct-assembly", dest="nstruct_assembly", help="number of assembly structures to create [default: %(default)s]", default=50000, type=int)
    parser_evaluate.add_argument("--cluster-cutoff", dest="cluster_cutoff", help="cluster cutoff in angström [default: %(default)s]", default=4.1, type=float)
    parser_evaluate.add_argument("--cluster-limit", dest="cluster_limit", help="maximum number of clusters to create [default: %(default)s]", default=10, type=int)
    parser_makeconstraints.add_argument("--dca-file", dest="dca_file", help="dca file to use as input [default: %(default)s]", default="dca/dca.txt")
    parser_makeconstraints.add_argument("--dca-count", dest="dca_count", help="maximum number o dca predictions to use [default: %(default)s]", default=100, type=int)
    parser_makeconstraints.add_argument("--mapping-mode", dest="mapping_mode", help="mapping mode to use for constraints creation [default: %(default)s", default="allAtomWesthof")
    parser_makeconstraints.add_argument("--filter", dest="filter", help="run dca contacts though (a) filter(s). For syntax information refer to the documentation. [default: %(default)s", default=None)
    parser_config.add_argument("key", help="config key")
    parser_config.add_argument("value", help="config value")
    parser_tools.add_argument("positional", nargs="*")

    for p in (parser_motifs, parser_assemble):
        p.add_argument("--seed", dest="seed", help="force random seed (when multithreading, this will be incremented for each process) [default: %(default)s]", default=None, type=int)
        p.add_argument("--use-native", dest="use_native", action="store_true", help="use native information for motif generation and assembly [default: %(default)s]")

    for p in (parser_preparecst, parser_motifs, parser_assemble, parser_editconstraints, parser_evaluate):
        p.add_argument("--cst", dest="cst", help="constraint selection [default: %(default)s]", default=None)

    for p in (parser_makeconstraints, parser_editconstraints):
        p.add_argument("--cst-function", dest="cst_function", help="rosetta function to use for the constraints [default: '%(default)s']", default="FADE -100 26 20 -2 2")
        p.add_argument("--cst-out-file", dest="cst_out_file", help="output cst file [default: inferred from input file]", default=None)


    # Process arguments
    args = parser.parse_args()

    # system configuration
    sysconfig = SysConfig()
    check = sysconfig.checkSysConfig()
    failed = check["fail"]
    if failed:
        for f in failed:
            printToStderr("Sysconfig error: Cannot access %s" % (f))
        return 1
    else:
        if not args.quiet:
            sysconfig.printSysConfig()

    try:
        if args.subcommand == "tools":
            # TODO: hack! introduce proper argument parsing for tools
            sys.argv = sys.argv[1:]
            tools.tools()
            return 0

        # for all the other subcommannds we need a simulation
        p = RNAPrediction(sysconfig, os.getcwd())

        # treat "config" special, because we want to print the configuartion after we change something
        if args.subcommand == "config":
            p.modifyConfig(args.key, args.value)
            p.saveConfig()

        if not args.quiet:
            p.printConfig()

        if args.subcommand == "status":
            p.printStatus()
        elif args.subcommand == "prepare":
            p.prepare(fasta_file=args.sequence, params_file=args.secstruct, native_pdb_file=args.native, name=args.name)
            p.saveConfig()
        elif args.subcommand == "prepare-cst":
            p.prepareCst(constraints=args.cst, motifsOverride=args.motifs_cst)
        elif args.subcommand == "make-constraints":
            p.makeConstraints(dcaPredictionFileName=args.dca_file, outputFileName=args.cst_out_file, numberDcaPredictions=args.dca_count, cstFunction=args.cst_function, filterText=args.filter, mappingMode=args.mapping_mode)
        elif args.subcommand == "edit-constraints":
            p.editConstraints(constraints=args.cst, outputFileName=args.cst_out_file,cstFunction=args.cst_function)
        elif args.subcommand == "create-helices":
            p.create_helices(dry_run=args.dry_run, threads=args.threads)
        elif args.subcommand == "create-motifs":
            p.create_motifs(dry_run=args.dry_run, nstruct=args.nstruct_motifs, cycles=args.cycles_motifs, seed=args.seed, use_native_information=args.use_native, threads=args.threads, constraints=args.cst)
        elif args.subcommand == "assemble":
            p.assemble(dry_run=args.dry_run, constraints=args.cst, nstruct=args.nstruct_assembly, cycles=args.cycles_assembly, seed=args.seed, use_native_information=args.use_native, threads=args.threads)
        elif args.subcommand == "evaluate":
            p.evaluate(constraints=args.cst, cluster_limit=args.cluster_limit, cluster_cutoff=args.cluster_cutoff)
        elif args.subcommand == "compare":
            p.compare()

        return 0

    except (SimulationException, DcaException) as e:
        printToStderr(e)
        return 1
    except KeyboardInterrupt:
        return 1
    except Exception, e:
        if DEBUG or TESTRUN:
            raise
        else:
            printToStderr(e)
        return 1

if __name__ == "__main__":
    sys.exit(main())
