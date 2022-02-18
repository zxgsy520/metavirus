#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import os.path
import shutil
import logging
import argparse

sys.path.append(os.path.join("/Work/user/zhangxg/pipeline/metavirus/v1.1.1/"))
from ngsmetavirus.config import *
from ngsmetavirus.common import check_path, mkdir, read_files, get_version
from dagflow import DAG, Task, ParallelTask, do_dag


LOG = logging.getLogger(__name__)
__version__ = "1.0.1"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []

SOFTWARE_VERSION = {
    "biotool":{
        "GETVER": "%s/biotool -h |grep 'version'" % BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "1.0.0",
    },
    "blastn":{
        "GETVER": "%s/blastn -version 2>&1| grep 'blastn'" % BLAST_BIN,
        "REGEXP": "\d+\.\d+\.\d+\+",
        "MINVER": "2.7.1+",
    },
}


def split_genome(genome, prefix, job_type, work_dir, size="5m", window="3kb"):

    dag = DAG("split_genome")
    task = Task(
        id="split_genome",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
python {script}/genome2short.py {genome} --window {window} --minlen 2000 --maxlen 2500 >short.fasta
{bin}/biotool seqsplit short.fasta --size {size} \\
  --prefix {prefix} --workdir {work_dir}
rm -rf short.fasta
""".format(script=SCRIPTS,
            bin=BIN,
            genome=genome,
            prefix=prefix,
            size=size,
            window=window,
            work_dir=work_dir,
        )
    )
    dag.add_task(task)
    do_dag(dag, 8, 10)

    genomes = read_files(work_dir, "%s.part*.fa" % prefix)
    tig2seq = os.path.join(work_dir, "contig2sequence.tsv")

    return genomes, tig2seq


def create_blstn_tasks(genomes, prefix, db, threads, job_type,
                       work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in genomes]
    id = "blstn"
    option = OrderedDict()
    option["blastn"] = {
        "version": get_version(SOFTWARE_VERSION["blastn"]),
        "option": ""
    }

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={blast}:$PATH
blastn -query {{genomes}} -db {db} \\
  -outfmt "6 std staxid sskingdom staxids" -max_target_seqs 1 \\
  -num_threads {threads} -out {{prefixs}}.blast.m6
""".format(blast=BLAST_BIN,
           db=db,
           threads=threads),
        prefixs=prefixs,
        genomes=genomes,
    )

    join_task = Task(
        id="merge_%s" % id,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.blast.m6 > {prefix}.m6
python {scripts}/obtain_taxonomy.py {prefix}.m6 -t {taxonomy} >{prefix}.species_annotation.txt
python {scripts}/contig_taxonomy.py {prefix}.species_annotation.txt >{prefix}.stat_contig_taxonomy.tsv
""".format(id=id,
           scripts=SCRIPTS,
           taxonomy=TAXONOMY,
           prefix=prefix,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task, option


def run_genome_blstn(prefix, genome, threads, job_type, work_dir, out_dir,
                  concurrent, refresh):

    genome = check_path(genome)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    work_dict = {
        "split": "00_data",
        "blstn": "blstn",
    }
    for k, v in work_dict.items():
        work_dict[k] = mkdir(os.path.join(work_dir, v))
    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    genomes, tig2seq = split_genome(genome=genome,
        prefix=prefix,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["split"]),
        size="5m",
        window="5kb"
    )
    dag = DAG("blastn")
    blastn_tasks, blastn_join, option = create_blstn_tasks(
        genomes=genomes,
        prefix=prefix,
        db=NT_DATABASE,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["blstn"],
        out_dir=out_dir
    )

    options["software"].update(option)
    dag.add_task(*blastn_tasks)
    dag.add_task(blastn_join)

    do_dag(dag, concurrent, refresh)

    return options


def genome_blstn(args):

    options = run_genome_blstn(
         genome=args.genome,
         prefix=args.prefix,
         threads=args.threads,
         job_type=args.job_type,
         work_dir=args.work_dir,
         out_dir=args.out_dir,
         concurrent=args.concurrent,
         refresh=args.refresh
    )
    with open(os.path.join(args.out_dir, "blastn.json"), "w") as fh:
         json.dump(options, fh, indent=2)

    return 0


def add_blastn_args(parser):

    parser.add_argument("genome", metavar="FILE", type=str,
        help="Genome in fasta file")
    parser.add_argument("-p", "--prefix", metavar="STR", type=str, default="G1",
        help="Name used in pipeline")
    parser.add_argument("-t", "--threads", metavar="INT", type=int, default=4,
        help="Threads used to run blastn (default: 4)")
    parser.add_argument("-c", "--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    parser.add_argument("-r", "--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    parser.add_argument("-j", "--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    parser.add_argument("-w", "--work_dir", metavar="DIR", type=str, default="work",
        help="Work directory (default: work)")
    parser.add_argument("-o", "--out_dir", metavar="DIR", type=str, default="out",
        help="Output directory (default: out)")


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
name:
    blastn: Genome and database comparison

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_blastn_args(parser)
    args = parser.parse_args()
    genome_blstn(args)


if __name__ == "__main__":
    main()
