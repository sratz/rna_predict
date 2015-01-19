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
from rna_prediction.simulation import SimulationException
from rna_prediction.sysconfig import SysConfig
from .dcatools import DcaException
from . import dcatools

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

    # Setup argument parser
    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
    group_steps = parser.add_argument_group(title="main arguments (processing steps)")
    group_steps.add_argument("--prepare", dest="prepare", action="store_true", help="prepare stems and motifs [default: %(default)s]")
    group_steps.add_argument("--prepare-cst", dest="prepare_cst", action="store_true", help="prepare constraints file for motif generation and assembly [default: %(default)s]")
    group_steps.add_argument("--create-helices", dest="create_helices", action="store_true", help="create ideal a-helices [default: %(default)s]")
    group_steps.add_argument("--create-motifs", dest="create_motifs", action="store_true", help="create motifs [default: %(default)s]")
    group_steps.add_argument("--make-constraints", dest="make_constraints", action="store_true", help="create a constraints file from a dca prediction [default: %(default)s]")
    group_steps.add_argument("--edit-constraints", dest="edit_constraints", action="store_true", help="replace rosetta function in a constraints file [default: %(default)s]")
    group_steps.add_argument("--assemble", dest="assemble", action="store_true", help="assemble [default: %(default)s]")
    group_steps.add_argument("--extract", dest="extract", action="store_true", help="extract pdb data and scrore [default: %(default)s]")
    group_steps.add_argument("--evaluate", dest="evaluate", action="store_true", help="evaluate data (clusters) [default: %(default)s]")
    group_steps.add_argument("--compare", dest="compare", action="store_true", help="print comparison of prediction to native structure [default: %(default)s]")
    group_prepare = parser.add_argument_group(title="options for the --prepare step")
    group_prepare.add_argument("--name", dest="name", help="simulation name (for example the analyzed pdb) [default: basename of directory]")
    group_prepare.add_argument("--native", dest="native", help="native pdb [default: %(default)s]")
    group_prepare.add_argument("--sequence", dest="sequence", help="sequence fasta file [default: %(default)s]", default="sequence.fasta")
    group_prepare.add_argument("--secstruct", dest="secstruct", help="secondary structure file [default: %(default)s]", default="secstruct.txt")
    group_preparecst = parser.add_argument_group(title="options for the --prepare-cst step")
    group_preparecst.add_argument("--override-motifs-cst", dest="motifs_cst", help="use motifs from a different constraints set [default: %(default)s]", const="none", nargs="?")
    group_motifs = parser.add_argument_group(title="options for the --create-motifs step")
    group_motifs.add_argument("--cycles-motifs", dest="cycles_motifs", help="number of cycles for motif generation [default: %(default)s]", default=5000, type=int)
    group_motifs.add_argument("--nstruct-motifs", dest="nstruct_motifs", help="number of motif structures to create [default: %(default)s]", default=4000, type=int)
    group_assembly = parser.add_argument_group(title="options for the --assemble step")
    group_assembly.add_argument("--cycles-assembly", dest="cycles_assembly", help="number of cycles for assembly [default: %(default)s]", default=20000, type=int)
    group_assembly.add_argument("--nstruct-assembly", dest="nstruct_assembly", help="number of assembly structures to create [default: %(default)s]", default=50000, type=int)
    group_simulaion = parser.add_argument_group(title="options for the --create-motifs and --assemble steps")
    group_simulaion.add_argument("--seed", dest="seed", help="force random seed (when multithreading, this will be incremented for each process) [default: %(default)s]", default=None, type=int)
    group_simulaion.add_argument("--use-native", dest="use_native", action="store_true", help="use native information for motif generation and assembly [default: %(default)s]")
    group_evaluate = parser.add_argument_group(title="options for the --evaluate step")
    group_evaluate.add_argument("--cluster-cutoff", dest="cluster_cutoff", help="cluster cutoff in nm [default: %(default)s]", default=0.41, type=float)
    group_evaluate.add_argument("--cluster-limit", dest="cluster_limit", help="maximum number of clusters to create [default: %(default)s]", default=10, type=int)
    group_cst = parser.add_argument_group(title="constraint selection (for --assemble, --extract, --evaluate, --edit-constraint)")
    group_cst.add_argument("--cst", dest="cst", help="constraint file or name to use in assembly, extraction or evaluation steps [default: %(default)s]", default=None)
    group_makecst = parser.add_argument_group(title="options for --make-constraints")
    group_makecst.add_argument("--dca-file", dest="dca_file", help="dca file to use as input [default: %(default)s]", default="dca/dca.txt")
    group_makecst.add_argument("--dca-count", dest="dca_count", help="maximum number o dca predictions to use [default: %(default)s]", default=100, type=int)
    group_makecst.add_argument("--mapping-mode", dest="mapping_mode", help="mapping mode to use for constraints creation [default: %(default)s", default="allAtomWesthof")
    group_makecst.add_argument("--filter", dest="filter", help="run dca contacts though (a) filter(s). For syntax information refer to the documentation. [default: %(default)s", default=None)
    group_editmakecst = parser.add_argument_group(title="options for --make-constraints, --edit-constraints")
    group_editmakecst.add_argument("--cst-function", dest="cst_function", help="rosetta function to use for the constraints [default: '%(default)s']", default="FADE -100 26 20 -2 2")
    group_editmakecst.add_argument("--cst-out-file", dest="cst_out_file", help="output cst file [default: inferred from input file]", default=None)
    group_other = parser.add_argument_group(title="other arguments")
    group_other.add_argument('-j', '--threads', dest="threads", help="maximum number of parallel subprocesses [default: %(default)s]", default=1, type=int)
    group_other.add_argument('-V', '--version', action='version', version=program_version_message)
    group_other.add_argument('-q', "--quiet", dest="quiet", action="store_true", help="don't print config on start [default: %(default)s]")
    group_other.add_argument("-n", "--dry-run", dest="dry_run", action="store_true", help="don't execute and only print external commands [default: %(default)s]")
    group_other.add_argument("--config", dest="config", help="modify config variable", metavar=("KEY", "VALUE"), nargs=2)

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
        p = RNAPrediction(sysconfig, os.getcwd())

        if args.config:
            p.modifyConfig(args.config[0], args.config[1])
            p.saveConfig()

        if not args.quiet:
            p.printConfig()

        if not args.config and not args.prepare and not args.prepare_cst and not args.make_constraints and not args.edit_constraints and not args.create_helices and not args.create_motifs and not args.assemble and not args.extract and not args.evaluate and not args.compare:
            p.printStatus()

        if args.prepare:
            p.prepare(fasta_file=args.sequence, params_file=args.secstruct, native_pdb_file=args.native, name=args.name)
            p.saveConfig()
        if args.prepare_cst:
            p.prepareCst(constraints=args.cst, motifsOverride=args.motifs_cst)
        if args.make_constraints:
            p.makeConstraints(dcaPredictionFileName=args.dca_file, outputFileName=args.cst_out_file, numberDcaPredictions=args.dca_count, cstFunction=args.cst_function, filterText=args.filter, mappingMode=args.mapping_mode)
        if args.edit_constraints:
            p.editConstraints(constraints=args.cst, outputFileName=args.cst_out_file,cstFunction=args.cst_function)
        if args.create_helices:
            p.create_helices(dry_run=args.dry_run, threads=args.threads)
        if args.create_motifs:
            p.create_motifs(dry_run=args.dry_run, nstruct=args.nstruct_motifs, cycles=args.cycles_motifs, seed=args.seed, use_native_information=args.use_native, threads=args.threads, constraints=args.cst)
        if args.assemble:
            p.assemble(dry_run=args.dry_run, constraints=args.cst, nstruct=args.nstruct_assembly, cycles=args.cycles_assembly, seed=args.seed, use_native_information=args.use_native, threads=args.threads)
        if args.extract:
            p.extract(constraints=args.cst)
        if args.evaluate:
            p.evaluate(constraints=args.cst, cluster_limit=args.cluster_limit, cluster_cutoff=args.cluster_cutoff)
        if args.compare:
            p.compare()
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

    return 0

if __name__ == "__main__":
    sys.exit(main())
