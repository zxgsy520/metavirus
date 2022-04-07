#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import os.path
import logging
import argparse

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../"))
from ngsmetavirus.common import get_version, read_files, mkdir, check_path
from dagflow import Task, ParallelTask, DAG, do_dag
from ngsmetavirus.config import *
from seqkit.split import seq_split

LOG = logging.getLogger(__name__)
__version__ = "1.2.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []

VIBRANT_BIN = "/Work/pipeline/software/meta/vibrant/v1.2.1/bin/"
VIBRANT_DB = "/Work/pipeline/software/meta/vibrant/v1.2.1/share/vibrant-1.2.1/db/"
VIRFINDER_R = "/Work/pipeline/software/meta/virfinder/v1.1/bin/"
VIRKRAKEN_BIN = "/Work/pipeline/software/Base/miniconda/v4.10.3/bin/"
KRAKEN2_DB = "/Work/user/liyilin/database/kraken"

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


def create_vibrant_tasks(genomes, prefix, thread=10, job_type="sge",
                               work_dir=".", out_dir="."):

    prefixs = [os.path.basename(i) for i in genomes]
    option = OrderedDict()
    option["vibrant"] = {
        "version": "v1.2.1",
        "option": "default"
    }
    id = "vibrant"
    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""\
export PATH={vibrant}:$PATH
VIBRANT_run.py -i {{genome}} -folder ./ \\
  -d {vibrant_db}/databases -m {vibrant_db}/files -t {thread}
""".format(vibrant=VIBRANT_BIN,
            vibrant_db=VIBRANT_DB,
            thread=thread
            ),
        genome=genomes,
    )

    join_task = Task(
        id="merge_vibrant",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""\
python {script}/read_vibrant_phage.py vibrant_*/VIBRANT_*/VIBRANT_phages_*/*.fna > {prefix}.vibrant.tsv
""".format(script=SCRIPTS,
            prefix=prefix
        )
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task, option, os.path.join(work_dir, "%s.vibrant.tsv" % prefix)



def create_virfinder_tasks(genomes, prefix, thread=10, job_type="sge",
                               work_dir=".", out_dir="."):

    prefixs = [os.path.basename(i) for i in genomes]
    option = OrderedDict()
    option["virfinder"] = {
        "version": "v1.1",
        "option": 'default'
    }
    id = "virfinder"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""\
export PATH={virfinder}:$PATH
if [ ! -e {{prefix}}.virfinder.tsv ]; then
  {virfinder}/Rscript {script}/virfinder.R {{genome}} |grep -v score >{{prefix}}.virfinder.tsv
fi
""".format(virfinder=VIRFINDER_R,
            script=SCRIPTS,
            ),
        genome=genomes,
        prefix=prefixs
    )

    join_task = Task(
        id="merge_virfinder",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 -V ",
        script="""\
cat {id}*/*.virfinder.tsv > {prefix}.virfinder.txt
""".format(id=id,
           prefix=prefix,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task, option, os.path.join(work_dir, "%s.virfinder.txt" % prefix)


def create_virkraken_task(genome, prefix, thread=10, job_type="sge",
                               work_dir=".", out_dir="."):
    option = OrderedDict()
    option["virkraken"] = {
        "version": "-",
        "option": "default"
    }
    if genome.endswith(".gz"):
        x = "--gzip-compressed"
    else:
        x = ""

    task = Task(
        id="virkraken",
        work_dir= work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""\
export PATH={kraken}:{virkraken}:$PATH
kraken2 -db {kraken_db} --threads {thread} --report {prefix}.kraken.report \\
  --output  {prefix}.kraken.out {x} {genome}
virkraken -f {prefix}.kraken.out -c {genome} -o {prefix}.virus
cut -d "," -f 1 {prefix}.virus.csv >{prefix}.virus.id
""".format(kraken=KRAKEN2_BIN,
            virkraken=VIRKRAKEN_BIN,
            kraken_db=KRAKEN2_DB,
            genome=genome,
            thread=thread,
            prefix=prefix,
            x=x
        )
    )

    return task, option, os.path.join(work_dir, "%s.virus.id" % prefix)


def create_checkv_task(genome, prefix,  virsorter, deepvirfinder, vibrant, virfinder, virkraken,
                       minlen="10kb", score=0.9, pvalue=0.01, thread=10,
                       job_type="sge", work_dir=".", out_dir="."):

    option = OrderedDict()
    option["checkv"] = {
        "version": get_version(SOFTWARE_VERSION["checkv"]),
        "option": "default"
    }
    if virkraken:
        virkraken = "--virkraken %s" % virkraken
    else:
        virkraken = ""
    task = Task(
        id="checkv",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s " % thread,
        script="""
python {script}/get_merge_virus.py {genome} \\
  --virsorter {virsorter} --deepvirfinder {deepvirfinder} \\
  --vibrant {vibrant} --virfinder {virfinder} \\
  {virkraken} --minlen {minlen} --score {score} \\
  --pvalue {pvalue} >{prefix}.virus.fasta 2>{prefix}.virus_summary.tsv
export PATH={checkv}:{bin}:$PATH
export CHECKVDB={checkv}/../../checkv-db-v1.0
rm -rf checkv_out
checkv end_to_end {prefix}.virus.fasta checkv_out -t {thread}
biotool stats {prefix}.virus.fasta >{prefix}.stat_virus.tsv
cp checkv_out/quality_summary.tsv {out_dir}/{prefix}.virus_quality_summary.tsv
cp {prefix}.virus.fasta {prefix}.virus_summary.tsv {prefix}.stat_virus.tsv {out_dir}
""".format(checkv=CHECKV_BIN,
            script=SCRIPTS,
            bin=BIN,
            genome=genome,
            virsorter=virsorter,
            deepvirfinder=deepvirfinder,
            vibrant=vibrant,
            virfinder=virfinder,
            virkraken=virkraken,
            minlen=minlen,
            score=score,
            pvalue=pvalue,
            prefix=prefix,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, option


def run_find_virus(genome, name, thread, job_type, concurrent, refresh,
                   work_dir, out_dir, minlen="10kb", score=0.9, pvalue=0.01):

    work_dir = mkdir(os.path.abspath(work_dir))
    out_dir = mkdir(os.path.abspath(out_dir))
    genome = check_path(genome)
    work_dict = {
        "split": "00_data",
        "sorter": "01_virsorter",
        "finder": "02_deepvirfinder",
        "vibrant": "03_vibrant",
        "virfinder": "04_virfinder",
        "virkraken": "05_virkraken"
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
    kraken_task, option, virkraken_tsv = create_virkraken_task(
        genome=genome,
        prefix=name,
        thread=thread,
        job_type=job_type,
        work_dir=work_dict["virkraken"],
        out_dir=out_dir
    )
    options["software"].update(option)
    dag.add_task(kraken_task)

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

    vibrant_tasks, vibrant_join, option, vibrant_tsv = create_vibrant_tasks(
        genomes=genomes,
        prefix=name,
        thread=thread,
        job_type=job_type,
        work_dir=work_dict["vibrant"],
        out_dir=out_dir
    )
    options["software"].update(option)
    dag.add_task(*vibrant_tasks)
    dag.add_task(vibrant_join)

    virfinder_tasks, virfinder_join, option, virfinder =create_virfinder_tasks(
        genomes=genomes,
        prefix=name,
        thread=thread,
        job_type=job_type,
        work_dir=work_dict["virfinder"],
        out_dir=out_dir
    )
    options["software"].update(option)
    dag.add_task(*virfinder_tasks)
    dag.add_task(virfinder_join)


    checkv_task, option = create_checkv_task(
        genome=genome,
        prefix=name,
        virsorter=raw_viral,
        deepvirfinder=dvfpred,
        vibrant=vibrant_tsv,
        virfinder=virfinder,
        virkraken=virkraken_tsv,
        minlen=minlen,
        score=score,
        pvalue=pvalue,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    options["software"].update(option)
    dag.add_task(checkv_task)
    checkv_task.set_upstream(kraken_task)
    checkv_task.set_upstream(sorter_join)
    checkv_task.set_upstream(finder_join)
    checkv_task.set_upstream(virfinder_join)
    checkv_task.set_upstream(vibrant_join)

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
         refresh=args.refresh,
         minlen=args.minlen,
        score=args.score,
        pvalue=args.pvalue,
    )
    with open(os.path.join(args.out_dir, "find_virus.json"), "w") as fh:
         json.dump(options, fh, indent=2)

    return 0


def add_find_virus_args(parser):

    parser.add_argument("genome", metavar="FILE",
        help="Genome in fasta file")
    parser.add_argument("-n", "--name", metavar="STR", required=True,
        help="Name used in pipeline")
    parser.add_argument("-ml", "--minlen", metavar="SRT", type=str, default="10kb",
        help="Set the minimum length of sequence filtering, default=10kb")
    parser.add_argument("--score", metavar="FLOAT", type=float, default=0.9,
        help="set score, default=0.9")
    parser.add_argument("--pvalue", metavar="FLOAT", type=float, default=0.01,
        help="Set pvalue, default=0.01")
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
