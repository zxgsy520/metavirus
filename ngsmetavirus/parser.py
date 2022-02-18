#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ngsmetavirus.config import *

__all__ = ["add_mngs_qc_args", "add_rm_host_args", "add_mngs_asm_args", 
    "add_annotation_args", "add_mngs_all_args"]

def add_workflow_args(parser):
    """
    add workflow arguments to parser
    :param parser: argparse object
    :return: parser
    """

    workflow_group = parser.add_argument_group(title="Workflow arguments", )
    workflow_group.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    workflow_group.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    workflow_group.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    workflow_group.add_argument("--work_dir", metavar="DIR", default=".",
        help="Work directory (default: current directory)")
    workflow_group.add_argument("--out_dir", metavar="DIR", default=".",
        help="Output directory (default: current directory)")

    return parser


def add_mngs_qc_args(parser):

    parser.add_argument("-p", "--prefix", metavar='STR', type=str, default='MNGS',
        help="Input the name of the sample, default=MNGS.")
    parser.add_argument("-r1", "--read1", metavar="FILE", type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-r2", "--read2", metavar="FILE", type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-ref", "--reference", metavar="FILE", type=str, default="",
        help="Input the host's reference database.")
    parser.add_argument('--nohost', action='store_true',
        help='Input the reference database is not the host.')
    parser.add_argument("-dt", "--dtype", metavar='STR', type=str,
        choices=["mgi", "illumina", "other"], default="illumina",
        help="Set up the sequencing platform of the data, default=illumina.")
    parser.add_argument("-at", "--atype", metavar='STR', type=str,
        choices=["metagenome", "metaviral", "rnaviral"], default="metagenome",
        help="""Set the type of analysis(metagenome, metaviral, rnaviral),\
              default=metagenome.""")
    parser.add_argument("--trim", metavar="INT", type=int, default=5,
        help="Set trim length, default=5")
    parser.add_argument("--qvalue", metavar="INT", type=int, default=20,
        help="The quality value that a base is qualified, default=20")
    parser.add_argument("--thread", metavar="INT", type=int, default=1,
        help="Analysis of the number of threads used, default=1")
    parser = add_workflow_args(parser)

    return parser


def add_rm_host_args(parser):

    parser.add_argument("-r1", "--read1", metavar="FILE", nargs='+', type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-r2", "--read2", metavar="FILE", nargs='+', type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-ref", "--reference", metavar="FILE", type=str, required=True,
        help="Input the host's reference database.")
    parser.add_argument('--stat', metavar='FILE', type=str, required=True,
        help='Input the statistical results of the original data.')
    parser.add_argument("--thread", type=int, default=1,
        help="threads used in genome mapping")
    parser = add_workflow_args(parser)

    return parser


def add_mngs_tax_args(parser):

    parser.add_argument("-p", "--prefix", metavar='STR', type=str, default='MNGS',
        help="Input the name of the sample, default=MNGS.")
    parser.add_argument("-r1", "--read1", metavar="FILE", nargs='+', type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-r2", "--read2", metavar="FILE", nargs='+', type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("--thread", type=int, default=1,
        help="threads used in genome mapping")
    parser = add_workflow_args(parser)

    return parser


def add_mngs_asm_args(parser):

    parser.add_argument("-p", "--prefix", metavar='STR', type=str, default='MNGS',
        help="Input the name of the sample, default=MNGS.")
    parser.add_argument("-r1", "--read1", metavar="FILE", type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-r2", "--read2", metavar="FILE", type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-at", "--atype", metavar='STR', type=str,
        choices=["metagenome", "metaviral", "rnaviral"], default="metagenome",
        help="""Set the type of analysis(metagenome, metaviral, rnaviral),\
              default=metagenome.""")
    parser.add_argument("-m", "--memory", metavar="INT", type=int, default=10,
        help="The memory that the program runs on ,default=10(G)")
    parser.add_argument("-t", "--thread", metavar="INT", type=int, default=1,
        help="threads used in genome assembly, default=1")
    parser = add_workflow_args(parser)

    return parser


def add_annotation_args(parser):

    parser.add_argument("-g", "--genomes", metavar="FILE", nargs='+', type=str, required=True,
        help="Input the assembled genome file.")
    parser.add_argument("--group", metavar="FILE", type=str, required=True,
        help="Input the grouping file of each sample.")
    parser.add_argument("-d", "--depths", metavar="FILE", nargs='+', type=str, required=True,
        help="Input genome depth statistics file.")
    parser.add_argument("-m", "--model", metavar="FILE", type=str,
        default="/export/personal/software/software/MetaGeneMark/v.3.38//MetaGeneMark_v1.mod",
        help="Model used for gene prediction.")
    parser.add_argument("--thread", metavar="INT", type=int, default=6,
        help="Set the number of threads used by the process, default=6")

    parser = add_workflow_args(parser)

    return parser


def add_mngs_all_args(parser):

    parser.add_argument("-m", "--memory", metavar="INT", type=int, default=10,
        help="The memory that the program runs on ,default=10(G)")
    parser = add_mngs_qc_args(parser)

    return parser
