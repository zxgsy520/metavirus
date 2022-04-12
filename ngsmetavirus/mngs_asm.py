#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import logging
import argparse

from collections import OrderedDict
from ngsmetavirus.config import *
from ngsmetavirus.parser import add_mngs_asm_args
from dagflow import DAG, Task, do_dag
from ngsmetavirus.common import check_path, mkdir, get_version


LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.0.0"


def create_megahit_task(prefix, read1, read2, thread, memory, job_type,
                        work_dir, out_dir):

    option = OrderedDict()
    option["megahit"] = {
        "version": get_version(SOFTWARE_VERSION["megahit"]),
        "option:": "--memory %s --num-cpu-threads %s" % (memory, thread)
    }
    option["spades"] = {
        "version": get_version(SOFTWARE_VERSION["spades"]),
        "option": "off"
    }

    task = Task(
        id="megahit",
        work_dir=work_dir,
        type=job_type,
        option="-l vf=%dG -pe smp %s %s" % (memory, thread, QUEUE),
        script="""
export PATH={megahit}:{bin}:$PATH
megahit -1 {read1} -2 {read2} \\
    --min-contig-len 500 --k-list 31,41,63,87,113,141 \\
    --num-cpu-threads {thread} --out-dir {prefix}_megahit
biotool sort_genome {prefix}_megahit/final.contigs.fa --minlen 50 > {prefix}.genome.fasta
biotool stats {prefix}.genome.fasta >{prefix}.stat_genome.tsv
#rm -rf {prefix}_megahit
cp {prefix}.genome.fasta {prefix}.stat_genome.tsv {out_dir}
""".format(megahit=MEGAHIT_BIN,
            bin=BIN,
            scripts=SCRIPTS,
            read1=read1,
            read2=read2,
            prefix=prefix,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, option, os.path.join(work_dir, "%s.genome.fasta" % prefix)


def create_spades_task(prefix, read1, read2, atype, thread, job_type,
                        work_dir, out_dir):

    option = OrderedDict()
    option["megahit"] = {
        "version": get_version(SOFTWARE_VERSION["megahit"]),
        "option:": "off"
    }
    option["spades"] = {
        "version": get_version(SOFTWARE_VERSION["spades"]),
        "option": "--%s --threads %s" % (atype, thread)
    }

    task = Task(
        id="spades",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={spades}:$PATH
spades.py --{atype} -1 {read1} -2 {read2} --threads {thread} -o {prefix}_spades
{bin}/biotool sort_genome {prefix}_spades/scaffolds.fasta --minlen 50 >{prefix}.genome.fasta
biotool stats {prefix}.genome.fasta >{prefix}.stat_genome.tsv
mv {prefix}_spades/contigs.fasta {prefix}.contigs.fasta
#mv {prefix}_spades/assembly_graph.fastg {prefix}.graph.fastg
#rm -rf {prefix}_spades
cp {prefix}.genome.fasta {prefix}.stat_genome.tsv {out_dir}
""".format(spades=SPADES_BIN,
            bin=BIN,
            scripts=SCRIPTS,
            read1=read1,
            read2=read2,
            atype=atype,
            prefix=prefix,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, option, os.path.join(work_dir, "%s.genome.fasta" % prefix)


def run_find_virus_task(genome, prefix, thread, job_type, work_dir, out_dir):

    task = Task(
        id="find_virus",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
python {root}/subpipe/find_virus.py {genome} --name {prefix} \\
  --thread {thread} --job_type {job_type} --work_dir {work_dir} --out_dir {out_dir}
""".format(root=ROOT,
            genome=genome,
            prefix=prefix,
            thread=thread,
            job_type=job_type,
            work_dir=work_dir,
            out_dir=out_dir
        )
    )

    return task


def run_mngs_asm(prefix, read1, read2, atype, thread, memory, job_type,
            concurrent, refresh, work_dir, out_dir):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    read1 = check_path(read1)
    read2 = check_path(read2)

    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }
    dag = DAG("mngs_asm")

    if atype in ["metaviral", "metaplasmid", "rnaviral", "plasmid"]:
        task, option, genome = create_spades_task(
            prefix=prefix,
            read1=read1,
            read2=read2,
            atype=atype,
            thread=thread,
            job_type=job_type,
            work_dir=work_dir,
            out_dir=out_dir,
        )
    else:
        task, option, genome = create_megahit_task(
           prefix=prefix,
           read1=read1,
           read2=read2,
           memory=memory,
           thread=thread,
           job_type=job_type,
           work_dir=work_dir,
           out_dir=out_dir,
        )
    options["software"].update(option)
    dag.add_task(task)
    
    if atype in ["metaviral", "rnaviral"]:
        find_virus_task = run_find_virus_task(
           genome=genome,
           prefix=prefix,
           thread=thread,
           job_type=job_type,
           work_dir=work_dir,
           out_dir=out_dir
        )
        dag.add_task(find_virus_task)
        find_virus_task.set_upstream(task)
    do_dag(dag, concurrent, refresh)

    return options, genome


def mngs_asm(args):

    options, genome = run_mngs_asm(
        prefix=args.prefix,
        read1=args.read1,
        read2=args.read2,
        atype=args.atype,
        thread=args.thread,
        memory=args.memory,
        job_type=args.args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
    )

    with open(os.path.join(args.out_dir, "mngs_asm.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    return 0


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_mngs_asm_args(parser)
    args = parser.parse_args()
    mngs_asm(args)


if __name__ == "__main__":
    main()
