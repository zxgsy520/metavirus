#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import os.path
import logging
import argparse


from ngsmetavirus import __author__, __email__, __version__
from ngsmetavirus.common import get_version, read_files, mkdir, check_path
from dagflow import Task, ParallelTask, DAG, do_dag
from ngsmetavirus.config import *
from seqkit.split import seq_split

LOG = logging.getLogger(__name__)
__all__ = []


def create_virsorter_tasks(genomes, name, thread=10, job_type="sge",
                          work_dir=".", out_dir="."):

    prefix = [os.path.basename(i) for i in genomes]
    option = OrderedDict()
    option["virsorter"] = {
        "version": get_version(SOFTWARE_VERSION["virsorter"]),
        "option": 'run --include-groups "dsDNAphage,lavidaviridae,NCLDV,RNA,ssDNA"'
    }
    id = "virsorter"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""\
export PATH={virsorter}:$PATH
virsorter run --seqfile {{genome}} \\
  --include-groups "dsDNAphage,lavidaviridae,NCLDV,RNA,ssDNA" \\
  --db-dir {virsorter}/db \\
  -j {thread} --min-score 0.6 --min-length 500 --working-dir ./
""".format(virsorter=VIRSORTER_BIN,
           thread=thread
           ),
        genome=genomes,
        prefix=prefix
    )

    join_task = Task(
        id="merge_virsorter",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 -V ",
        script="""\
cat {id}*/final-viral-score.tsv > {name}.viral-score.tsv
cat {id}*/final-viral-combined.fa > {name}.raw_viral.fasta
#cp {name}.raw_viral.fasta {name}.viral-score.tsv {out_dir}
""".format(id=id,
           script=SCRIPTS,
           name=name,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task, option, os.path.join(work_dir, "%s.viral-score.tsv" % name)


def create_deepvirfinder_tasks(genomes, name, thread=10, job_type="sge",
                               work_dir=".", out_dir="."):

    prefix = [os.path.basename(i) for i in genomes]
    option = OrderedDict()
    option["deepvirfinder"] = {
        "version": get_version(SOFTWARE_VERSION["deepvirfinder"]),
        "option": '--len 500'
    }
    id = "deepvirfinder"
    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""\
export MKL_THREADING_LAYER=GNU
export PATH={deepvirfinder}:$PATH
dvf.py -i {{genome}} --len 500 -c {thread} -o {{prefix}}
""".format(deepvirfinder=DEEPVIRFINDER_BIN,
           thread=thread
           ),
        genome=genomes,
        prefix=prefix
    )

    join_task = Task(
        id="merge_deepvirfinder",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 -V ",
        script="""\
cat {id}*/*/*_dvfpred.txt > {name}.dvfpred.txt
""".format(id=id,
           name=name,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task, option, os.path.join(work_dir, "%s.dvfpred.txt" % name)


def run_blast_task(genome, prefix, thread=10, job_type="sge",
                   work_dir=".", out_dir="."):

    task = Task(
        id="blast",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s " % thread,
        script="""
{root}/subpipe/genome_blastn.py {genome} --threads {thread} --job_type {job_type} \\
  --prefix {prefix} --work_dir ./ --out_dir ./
""".format(root=ROOT,
            thread=thread,
            genome=genome,
            prefix=prefix,
            job_type=job_type
        )
    )

    return task, os.path.join(work_dir, "01_blstn/%s.stat_contig_taxonomy.tsv" % prefix)


def create_checkv_task(genome, name, virsorter, deepvirfinder, tigtax,
                       thread=10, job_type="sge", work_dir=".", out_dir="."):

    option = OrderedDict()
    option["checkv"] = {
        "version": get_version(SOFTWARE_VERSION["checkv"]),
        "option": "default"
    }

    task = Task(
        id="checkv",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s " % thread,
        script="""
python {script}/get_merge_virus.py {genome} \\
  --virsorter {virsorter} --deepvirfinder {deepvirfinder} \\
  --blstn {tigtax} >{name}.raw_virus.fasta 2>{name}.virus_summary.tsv
export PATH={checkv}:$PATH
export CHECKVDB={checkv}/../../checkv-db-v1.0
rm -rf checkv_out
checkv end_to_end {name}.raw_virus.fasta checkv_out -t {thread}
""".format(checkv=CHECKV_BIN,
            script=SCRIPTS,
            genome=genome,
            virsorter=virsorter,
            deepvirfinder=deepvirfinder,
            tigtax=tigtax,
            name=name,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, option


def run_find_virus(genome, name, thread, job_type, concurrent, refresh,
                   work_dir, out_dir):

    work_dir = mkdir(os.path.abspath(work_dir))
    out_dir = mkdir(os.path.abspath(out_dir))
    genome = check_path(genome)
    work_dict = {
        "split": "00_data",
        "sorter": "01_virsorter",
        "finder": "02_deepvirfinder",
        "blast": "03_blast"
    }
    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    for k, v in work_dict.items():
        work_dict[k] = mkdir(os.path.join(work_dir, v))

    genomes = seq_split(
        filenames=[genome],
        mode="length",
        num=10000000,
        output_dir=os.path.join(work_dir, work_dict["split"]),
        minlen=1000
    )

    dag = DAG("find_virus")

    sorter_tasks, sorter_join, option, raw_viral = create_virsorter_tasks(
        genomes=genomes,
        name=name,
        thread=thread,
        job_type=job_type,
        work_dir=work_dict["sorter"],
        out_dir=out_dir
    )
    options["software"].update(option)
    dag.add_task(*sorter_tasks)
    dag.add_task(sorter_join)

    finder_tasks, finder_join, option, dvfpred =create_deepvirfinder_tasks(
        genomes=genomes,
        name=name,
        thread=thread,
        job_type=job_type,
        work_dir=work_dict["finder"],
        out_dir=out_dir
    )
    options["software"].update(option)
    dag.add_task(*finder_tasks)
    dag.add_task(finder_join)

    blast_task, tigtax = run_blast_task(
        genome=genome,
        prefix=name,
        thread=thread,
        job_type=job_type,
        work_dir=work_dict["blast"],
        out_dir=out_dir
    )
    dag.add_task(blast_task)

    checkv_task, option = create_checkv_task(
        genome=genome,
        name=name,
        virsorter=raw_viral,
        deepvirfinder=dvfpred,
        tigtax=tigtax,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    options["software"].update(option)
    dag.add_task(checkv_task)
    checkv_task.set_upstream(sorter_join)
    checkv_task.set_upstream(finder_join)
    checkv_task.set_upstream(blast_task)

    do_dag(dag, concurrent, refresh)

    return options


def find_virus(args):

    options = run_find_virus(
         genome=args.genome,
         name=args.name,
         thread=args.thread,
         job_type=args.job_type,
         work_dir=args.work_dir,
         out_dir=args.out_dir,
         concurrent=args.concurrent,
         refresh=args.refresh
    )
    with open(os.path.join(args.out_dir, "find_virus.json"), "w") as fh:
         json.dump(options, fh, indent=2)

    return 0


def add_find_virus_args(parser):

    parser.add_argument("genome", metavar="FILE",
        help="Genome in fasta file")
    parser.add_argument("-n", "--name", metavar="STR", required=True,
        help="Name used in pipeline")
    parser.add_argument("-t", "--thread", metavar="INT", type=int, default=4,
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
    find_virus: Predict virus sequence

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_find_virus_args(parser)
    args = parser.parse_args()
    find_virus(args)


if __name__ == "__main__":
    main()
