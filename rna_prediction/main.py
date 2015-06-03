#!/usr/bin/env python
# coding: utf-8

import argparse
from argparse import ArgumentParser
from argparse import HelpFormatter
import functools
import os
import sys

from .dcatools import DcaException
from .simulation import RNAPrediction
from .simulation import SimulationException
from .sysconfig import SysConfig
from . import tools
from . import utils
from . import __version__


def main():
    """rna_prediction main commandline parser."""

    # use line-buffered stdout and stderr
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)
    sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 1)

    def print_to_stderr(message, prefix="error"):
        sys.stderr.write("%s: %s\n" % (prefix, message))

    class SubcommandHelpFormatter(HelpFormatter):
        def _format_action(self, action):
            parts = super(SubcommandHelpFormatter, self)._format_action(action)
            if action.nargs == argparse.PARSER:
                parts = "\n".join(parts.split("\n")[1:])
            return parts

    def comma_separated_ranges_to_list_parser(s):
        try:
            return utils.comma_separated_ranges_to_list(s)
        except ValueError as ve:
            raise argparse.ArgumentTypeError(ve.message)

    def comma_separated_entries_to_dict_parser(s, type_key, type_value):
        try:
            return utils.comma_separated_entries_to_dict(s, type_key, type_value)
        except ValueError as ve:
            raise argparse.ArgumentTypeError(ve.message)

    # Setup argument parser
    parser = ArgumentParser(formatter_class=SubcommandHelpFormatter)
    subparsers = parser.add_subparsers(dest="subcommand", title="subcommands", metavar="subcommand")

    parser.add_argument('-j', '--threads', dest="threads", help="maximum number of parallel subprocesses [default: %(default)s]", default=1, type=int)
    parser.add_argument('-V', '--version', action='version', version="%%(prog)s version %s" % __version__)
    parser.add_argument('-q', "--quiet", dest="quiet", action="store_true", help="don't print config on start [default: %(default)s]")
    parser.add_argument("-n", "--dry-run", dest="dry_run", action="store_true", help="don't execute and only print external commands [default: %(default)s]")  # TODO: actually make dry-run work, not just for external subcommands

    parser_prepare = subparsers.add_parser("prepare", help="prepare stems and motifs")
    parser_preparecst = subparsers.add_parser("prepare-cst", help="prepare constraints file for motif generation and assembly")
    subparsers.add_parser("create-helices", help="create ideal a-helices")
    parser_motifs = subparsers.add_parser("create-motifs", help="create motifs")
    parser_assemble = subparsers.add_parser("assemble", help="assemble models")
    parser_evaluate = subparsers.add_parser("evaluate", help="evaluate data (clustering and rmsd calculation)")
    parser_evaluatecustom = subparsers.add_parser("evaluate-custom", help="evaluate using custom scoring")
    parser_makeconstraints = subparsers.add_parser("make-constraints", help="create a constraints file from a dca prediction")
    parser_editconstraints = subparsers.add_parser("edit-constraints", help="replace rosetta function in a constraints file")
    parser_status = subparsers.add_parser("status", help="print status summary")
    parser_config = subparsers.add_parser("config", help="modify config variable")
    parser_printmodels = subparsers.add_parser("print-models", help="print models")
    parser_extractmodels = subparsers.add_parser("extract-models", help="extract models as pdb files")

    parser_prepare.add_argument("--name", dest="name", help="simulation name (for example the analyzed pdb) [default: basename of directory]")
    parser_prepare.add_argument("--native", dest="native", help="native pdb [default: %(default)s]")
    parser_prepare.add_argument("--sequence", dest="sequence", help="sequence fasta file [default: %(default)s]", default="sequence.fasta")
    parser_prepare.add_argument("--secstruct", dest="secstruct", help="secondary structure file [default: %(default)s]", default="secstruct.txt")
    parser_preparecst.add_argument("--override-motifs-cst", dest="motifs_cst", help="use motifs from a different constraints set [default: %(default)s]", const="none", nargs="?")
    parser_motifs.add_argument("--cycles", dest="cycles", help="number of cycles for motif generation [default: %(default)s]", default=5000, type=int)
    parser_motifs.add_argument("--nstruct", dest="nstruct", help="number of motif structures to create [default: %(default)s]", default=4000, type=int)
    parser_motifs.add_argument("--motif-subset", dest="motif_subset", help="list of motifs to create models for [default: all motifs]", default=None, type=comma_separated_ranges_to_list_parser)
    parser_assemble.add_argument("--cycles", dest="cycles", help="number of cycles for assembly [default: %(default)s]", default=20000, type=int)
    parser_assemble.add_argument("--nstruct", dest="nstruct", help="number of assembly structures to create [default: %(default)s]", default=50000, type=int)
    parser_evaluate.add_argument("--cluster-cutoff", dest="cluster_cutoff", help="cluster cutoff in angström [default: %(default)s]", default=4.1, type=float)
    parser_evaluate.add_argument("--cluster-limit", dest="cluster_limit", help="maximum number of clusters to create [default: %(default)s]", default=10, type=int)
    parser_evaluate.add_argument("--full-eval", dest="full_eval", action="store_true", help="force full evaulation (scores and rmsd) instead of clustering only [default: %(default)s]")
    parser_status.add_argument("--compare", dest="compare", action="store_true", help="print rmsd comparison to native structure [default: %(default)s]")
    parser_status.add_argument("csts", metavar="csts", help="predictions to display [default: all]", nargs="*")

    for p in (parser_makeconstraints, parser_evaluatecustom):
        p.add_argument("--dca-file", dest="dca_file", help="dca file to use as input [default: %(default)s]", default="dca/dca.txt")
        p.add_argument("--dca-count", dest="dca_count", help="maximum number of dca predictions to use [default: %(default)s]", default=100, type=int)

    parser_makeconstraints.add_argument("--mapping-mode", dest="mapping_mode", help="mapping mode to use for constraints creation [default: %(default)s", default="allAtomWesthof")
    parser_makeconstraints.add_argument("--filter", dest="filter", help="run dca contacts though (a) filter(s). For syntax information refer to the documentation. [default: %(default)s", default=None)
    parser_config.add_argument("key", help="config key")
    parser_config.add_argument("value", help="config value")

    parser_evaluatecustom.add_argument("--threshold", help="threshold to be treated as contact [default: %(default)s]", default=7.5, type=float)
    parser_evaluatecustom.add_argument("--radius", help="number of neighboring residues in each direction to take into account [default: %(default)s]", default=2, type=int)
    parser_evaluatecustom.add_argument("--full-eval", dest="full_eval", action="store_true", help="force full evaulation (scores and rmsd) instead of clustering only [default: %(default)s]")

    for p in (parser_motifs, parser_assemble):
        p.add_argument("--seed", dest="seed", help="force random seed (when multithreading, this will be incremented for each process) [default: %(default)s]", default=None, type=int)
        p.add_argument("--use-native", dest="use_native", action="store_true", help="use native information for motif generation and assembly [default: %(default)s]")

    for p in (parser_preparecst, parser_motifs, parser_assemble, parser_editconstraints, parser_evaluate, parser_printmodels, parser_extractmodels, parser_evaluatecustom):
        p.add_argument("--cst", dest="cst", help="constraint selection [default: %(default)s]", default=None)

    for p in (parser_makeconstraints, parser_editconstraints):
        p.add_argument("--cst-function", dest="cst_function", help="rosetta function to use for the constraints [default: '%(default)s']", default="FADE -100 26 20 -2 2")
        p.add_argument("--cst-out-file", dest="cst_out_file", help="output cst file [default: inferred from input file]", default=None)

    for p in (parser_printmodels, parser_extractmodels):
        p.add_argument("--mode", dest="model_mode", help="model selection mode (tag, top, ntop, cluster) [default: %(default)s]", default="tag")
        p.add_argument("model_list", metavar="model", help="model number or tag, depending on mode argument", nargs="+")

    # tools sub-parser and sub-sub-parser
    parser_tools = subparsers.add_parser("tools", help="various plot generation tools", formatter_class=SubcommandHelpFormatter)
    subparser_tools = parser_tools.add_subparsers(dest="subcommand_tool", title="tools", metavar="tool")
    parser_tools_plot_contact_distances = subparser_tools.add_parser("plot-contact-distances", help="plot histogram for each nucleotide pair contact containing the distances of the atoms involved")
    parser_tools_plot_contact_atoms = subparser_tools.add_parser("plot-contact-atoms", help="plot atoms involved in forming nucleotide contacts")
    parser_tools_plot_constraint_quality = subparser_tools.add_parser("plot-constraint-quality", help="plot constraint quality")
    parser_tools_plot_dca_contacts_in_pdb = subparser_tools.add_parser("plot-dca-contacts-in-pdb", help="visualize how well DCA contacts are fullfiled in PDB files")
    parser_tools_plot_clusters = subparser_tools.add_parser("plot-clusters", help="plot score over native rmsd")
    parser_tools_plot_pdb_comparison = subparser_tools.add_parser("plot-pdb-comparison", help="compare PDB files by plotting the distance of the residues")
    parser_tools_plot_gdt = subparser_tools.add_parser("plot-gdt", help="create a global distance test plot")

    parser_tools_plot_contact_atoms.add_argument("--mean-cutoff", dest="mean_cutoff", help="limit for average distance [default: '%(default)s']", type=float, default=6.0)
    parser_tools_plot_contact_atoms.add_argument("--std-cutoff", dest="std_cutoff", help="limit for the standard deviation [default: '%(default)s']", type=float, default=3.0)

    parser_tools_plot_constraint_quality.add_argument("--dca-mode", dest="dca_mode", action="store_true", help="visualize residue-residue DCA instead of atom-atom constraints [default: '%(default)s']")
    parser_tools_plot_constraint_quality.add_argument("comparison_pdb", metavar="comparison-pdb", help="pdb filename to compare to")
    parser_tools_plot_constraint_quality.add_argument("sources", help="list of dca files, constraints, or 'filter:...'", nargs="+")

    parser_tools_plot_dca_contacts_in_pdb.add_argument("dca_file", metavar="dca-file",  help="dca file to use as input [default: %(default)s]", default="dca/dca.txt")
    parser_tools_plot_dca_contacts_in_pdb.add_argument("pdb_files", metavar="pdb-files", help="list of pdb files", nargs="+")

    parser_tools_plot_clusters.add_argument("cst", help="constraint selection")
    parser_tools_plot_clusters.add_argument("--max-models", dest="max_models", help="limit to number of models if > 1, or relative percentage if <= 1 [default: %(default)s]", type=float, default=0.99)
    parser_tools_plot_clusters.add_argument("--score-weights", dest="score_weights", help="alternative rosetta score weighting [default: %(default)s]", type=functools.partial(comma_separated_entries_to_dict_parser, type_key=str, type_value=float), default=None)

    parser_tools_plot_pdb_comparison.add_argument("ref_pdb", metavar="ref-pdb", help="reference PDB filename")
    parser_tools_plot_pdb_comparison.add_argument("sample_pdbs", metavar="sample-pdbs", help="list of sample PDB filenames", nargs="+")

    parser_tools_plot_gdt.add_argument("ref_pdb", metavar="ref-pdb", help="reference PDB filename")
    parser_tools_plot_gdt.add_argument("sample_pdbs", metavar="sample-pdbs", help="list of sample PDB filenames", nargs="+")

    # Process arguments
    args = parser.parse_args()

    # system configuration
    sysconfig = SysConfig()
    check = sysconfig.check_sysconfig()
    failed = check["fail"]
    if failed:
        for f in failed:
            print_to_stderr("Sysconfig error: Cannot access %s" % f)
        return 1
    else:
        if not args.quiet:
            sysconfig.print_sysconfig()

    try:
        if args.subcommand == "tools":
            if args.subcommand_tool == "plot-contact-distances":
                tools.plot_contact_distances()
            elif args.subcommand_tool == "plot-contact-atoms":
                tools.plot_contact_atoms(mean_cutoff=args.mean_cutoff, std_cutoff=args.std_cutoff)
            elif args.subcommand_tool == "plot-constraint-quality":
                tools.plot_constraint_quality(comparison_pdb=args.comparison_pdb, sources=args.sources, dca_mode=args.dca_mode)
            elif args.subcommand_tool == "plot-dca-contacts-in-pdb":
                tools.plot_dca_contacts_in_pdb(dca_prediction_filename=args.dca_file, pdb_files=args.pdb_files)
            elif args.subcommand_tool == "plot-clusters":
                tools.plot_clusters(cst=args.cst, max_models=args.max_models, score_weights=args.score_weights)
            elif args.subcommand_tool == "plot-pdb-comparison":
                tools.plot_pdb_comparison(pdb_ref_filename=args.ref_pdb, pdbs_sample_filenames=args.sample_pdbs)
            elif args.subcommand_tool == "plot-gdt":
                tools.plot_gdt(pdb_ref_filename=args.ref_pdb, pdbs_sample_filenames=args.sample_pdbs)
            return 0

        # for all the other subcommannds we need a simulation
        p = RNAPrediction(sysconfig)

        # treat "config" special, because we want to print the configuartion after we change something
        if args.subcommand == "config":
            p.modify_config(args.key, args.value)
            p.save_config()

        if not args.quiet:
            p.print_config()

        if args.subcommand == "status":
            p.print_status(native_compare=args.compare, csts=args.csts)
        elif args.subcommand == "prepare":
            p.prepare(fasta_file=args.sequence, params_file=args.secstruct, native_pdb_file=args.native, name=args.name)
            p.save_config()
        elif args.subcommand == "prepare-cst":
            p.prepare_cst(constraints=args.cst, motifs_override=args.motifs_cst)
        elif args.subcommand == "make-constraints":
            p.make_constraints(dca_prediction_filename=args.dca_file, output_filename=args.cst_out_file, number_dca_predictions=args.dca_count, cst_function=args.cst_function, filter_text=args.filter, mapping_mode=args.mapping_mode)
        elif args.subcommand == "edit-constraints":
            p.edit_constraints(constraints=args.cst, output_filename=args.cst_out_file, cst_function=args.cst_function)
        elif args.subcommand == "create-helices":
            p.create_helices(dry_run=args.dry_run, threads=args.threads)
        elif args.subcommand == "create-motifs":
            p.create_motifs(dry_run=args.dry_run, nstruct=args.nstruct, cycles=args.cycles, seed=args.seed, use_native_information=args.use_native, threads=args.threads, constraints=args.cst, motif_subset=args.motif_subset)
        elif args.subcommand == "assemble":
            p.assemble(dry_run=args.dry_run, constraints=args.cst, nstruct=args.nstruct, cycles=args.cycles, seed=args.seed, use_native_information=args.use_native, threads=args.threads)
        elif args.subcommand == "evaluate":
            p.evaluate(constraints=args.cst, cluster_limit=args.cluster_limit, cluster_cutoff=args.cluster_cutoff, full_evaluation=args.full_eval)
        elif args.subcommand == "print-models":
            p.print_models(constraints=args.cst, model_list=args.model_list, kind=args.model_mode)
        elif args.subcommand == "extract-models":
            p.extract_models(constraints=args.cst, model_list=args.model_list, kind=args.model_mode)
        elif args.subcommand == "evaluate-custom":
            p.evaluate_custom(constraints=args.cst, dca_prediction_filename=args.dca_file, full_evaluation=args.full_eval, threshold=args.threshold, radius=args.radius, number_dca_predictions=args.dca_count, threads=args.threads)

        return 0

    except (SimulationException, DcaException) as e:
        print_to_stderr(e)
        return 1
    except KeyboardInterrupt:
        return 1


if __name__ == "__main__":
    sys.exit(main())
