#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import glob
import logging
import argparse

from ngsmetavirus.config import *
from ngsmetavirus.common import check_path, mkdir, read_tsv
from dagflow import DAG, Task, do_dag

LOG = logging.getLogger(__name__)
__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def create_merge_data_task(prefix, read1, read2, work_dir, job_type="local"):

    if len(read1)==1:
        comm1 = """\
ln -s {read1} {prefix}.raw.R1.fq.gz
ln -s {read2} {prefix}.raw.R2.fq.gz
""".format(read1=read1[0], read2=read2[0], prefix=prefix)
    else:
        comm1 = """\
cat {read1} >{prefix}.raw.R1.fq.gz
cat {read2} >{prefix}.raw.R2.fq.gz
""".format(read1=" ".join(read1), read2=" ".join(read2), prefix=prefix)

    task = Task(
        id="merge_data_%s" % prefix,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
{comm1}
""".format(comm1=comm1)
    )
    read1 = os.path.join(work_dir, "%s.raw.R1.fq.gz" % prefix)
    read2 = os.path.join(work_dir, "%s.raw.R2.fq.gz" % prefix)

    return task, read1, read2


def create_mngs_task(prefix, read1, read2, reference, nohost, dtype, atype,
                     job_type, work_dir, out_dir, trim=0, project="", id=""):

    if nohost:
        nohost = "--nohost"
    else:
        nohost = ""
    if reference:
        rx = "--reference %s" % reference
    else:
        rx = ""

    task = Task(
        id="mngs_%s" % prefix,
        work_dir=work_dir,
        type="local",
        option="-pe smp 1",
        script="""
{root}/ngsmetavirus.py all \\
  --prefix {prefix} --read1 {read1} --read2 {read2} \\
  --dtype {dtype} --atype {atype} {rx} \\
  --project {project} --id {id} \\
  --trim {trim} --thread 6 --job_type {job_type} {nohost} \\
  --work_dir {work}  --out_dir {out}
""".format(root=ROOT,
            prefix=prefix,
            read1=read1,
            read2=read2,
            rx=rx,
            nohost=nohost,
            dtype=dtype,
            atype=atype,
            trim=trim,
            project=project,
            id=id,
            job_type=job_type,
            work=work_dir,
            out=out_dir
        )
    )

    return task


def run_mngs_multi(input, reference, nohost, dtype, atype, trim, job_type,
                   concurrent, refresh, work_dir, out_dir, project="", id=""):

    input = check_path(input)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    if not reference:
        pass
    else:
        try:
            reference = check_path(check_path)
        except:
            reference = check_path("%s.3.ht2" % reference)
            reference = reference.strip(".3.ht2")
        else:
            raise Exception("Reference genome %s does not exist" % reference)

    data = {}

    for line in read_tsv(input):
        if line[0] not in data:
            data[line[0]] = [[], []]
        data[line[0]][0].append(check_path(line[1]))
        data[line[0]][1].append(check_path(line[2]))

    dag = DAG("mngs_multi")

    for prefix in data:

        prefix_work = mkdir(os.path.join(work_dir, prefix))
        data_task, read1, read2 = create_merge_data_task(
            prefix=prefix,
            read1=data[prefix][0],
            read2=data[prefix][1],
            work_dir=prefix_work
        )
        dag.add_task(data_task)
        task = create_mngs_task(
            prefix=prefix,
            read1=read1,
            read2=read2,
            reference=reference,
            nohost=nohost,
            dtype=dtype,
            atype=atype,
            trim=trim,
            job_type=job_type,
            project=project,
            id=id,
            work_dir=prefix_work,
            out_dir=mkdir(os.path.join(out_dir, prefix))
        )
        dag.add_task(task)
        task.set_upstream(data_task)

    do_dag(dag, concurrent, refresh)

    return 0


def mngs_multi(args):

    run_mngs_multi(
        input=args.input,
        reference=args.reference,
        nohost=args.nohost,
        dtype=args.dtype,
        atype=args.atype,
        trim=args.trim,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        concurrent=args.concurrent,
        refresh=args.refresh,
        job_type=args.job_type,
        project=args.project,
        id=args.id
    )


def add_mngs_multi_args(parser):

    parser.add_argument("input", metavar='FILE', type=str,
        help="Input the reads list.")
    parser.add_argument("-ref", "--reference", metavar="FILE", type=str, default="",
        help="Input the host's reference database.")
    parser.add_argument('--nohost', action='store_true',
        help='Input the reference database is not the host.')
    parser.add_argument("-dt", "--dtype", metavar='STR', type=str,
        choices=["mgi", "illumina", "other"], default="illumina",
        help="Set up the sequencing platform of the data, default=illumina.")
    parser.add_argument("-at", "--atype", metavar='STR', type=str,
        choices=["metagenome", "metaviral", "rnaviral"], default="metagenome",
        help="""Set the type of analysis(metagenome, metavirus, rnaviral),\
              default=metagenome.""")
    parser.add_argument("--trim", metavar="INT", type=int, default=5,
        help="Set trim length, default=5")
    parser.add_argument("--project", metavar="STR", type=str, required=True,
        help="Input project name.")
    parser.add_argument("--id", metavar="STR", type=str, required=True,
        help="Input project id.")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10).")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30).")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local).")
    parser.add_argument("--work_dir", metavar="DIR", type=str, default=".",
        help="Work directory (default: current directory).")
    parser.add_argument("--out_dir", metavar="DIR", type=str, default=".",
        help="Output directory (default: current directory).")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
attention:
    ngsmetavirus.py.py multi input.list
File format：
name    R1  R2

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_mngs_multi_args(parser)
    args = parser.parse_args()
    mngs_multi(args)


if __name__ == "__main__":
    main()
