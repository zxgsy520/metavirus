#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import os.path
import shutil
import logging
import argparse

sys.path.append(os.path.join("/project/hgptemp/pipeline/v1.1.1/"))
from ngsmetavirus.config import *
from ngsmetavirus.common import check_path, mkdir, read_files, get_version
from dagflow import DAG, Task, ParallelTask, do_dag

LOG = logging.getLogger(__name__)
__version__ = "1.0.1"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []

TAXONONMY = "/project/hgptemp/database/taxonomy/202201/kraken.taxonomy.gz"
KRAKEN2_BIN = "/project/hgptemp/software/kraken/v2.1.2/bin/"
VIRKRAKEN_BIN = "/project/hgptemp/software/miniconda/v4.10.3/bin/"
SEQKIT_BIN = "/project/software/seqkit_dit/"
FLYE_BIN = "/project/hgptemp/software/flye/v2.9/bin/"

DATATYPE = OrderedDict([
    ("hifi", {
        "flye": "--meta --plasmids --iterations 2 --min-overlap 2500 --pacbio-hifi",
        "minimap2": "--secondary=no -ax asm5",
        "minimap2paf": "--secondary=no -x asm5"}
     ),
    ("clr", {
        "flye": "--meta --plasmids --iterations 2 --min-overlap 2500 --pacbio-raw",
        "minimap2": "--secondary=no -ax map-pb",
        "minimap2paf": "--secondary=no -x map-pb"}
     ),
    ("ont", {
        "flye": "--meta --plasmids --iterations 2 --min-overlap 2500 --nano-raw",
        "minimap2": "--secondary=no -ax map-ont",
        "minimap2paf": "--secondary=no -x map-ont"}
     ),
    ])

SOFTWARE_VERSION = {
    "kraken2": {
        "GETVER": "%s/kraken2 --version|grep 'version'" % KRAKEN2_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2.0.8",
    },
    "flye": {
        "GETVER": "%s/flye --version 2>&1" % FLYE_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "2.6",
    },
}

DB = "/project/hgptemp/database/Standard/20210517/standard"

def create_get_tgs_virus_task(reads, prefix, job_type, work_dir, out_dir, threads=10):

    if reads.endswith(".gz"):
        x = "--gzip-compressed"
    else:
        x = ""
    if reads.endswith("fq.gz") or reads.endswith("fq") or reads.endswith("fastq.gz") or reads.endswith("fastq"):
        f = "fastq"
    else:
        f = "fastq"

    task = Task(
        id="get_get_tgs_virus",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={kraken}:{virkraken}:$PATH
#kraken2 -db {db} --threads {threads} --report {prefix}.kraken.report \\
#  --output {prefix}.kraken.out {x} {reads}
#virkraken -f {prefix}.kraken.out -c {reads} -o {prefix}.virus
#cut -d "," -f 1 {prefix}.virus.csv >{prefix}.virus.id
python {scripts}/kraken2species.py {prefix}.kraken.out --taxonomy {taxonomy} >{prefix}.kraken.tax
grep Viruses {prefix}.kraken.tax |cut -f1 >{prefix}.virus.id
export PATH={seqkit}:$PATH
seqkit grep -f {prefix}.virus.id {reads} >{prefix}.virus_reads.{f}
""".format(kraken=KRAKEN2_BIN,
            virkraken=VIRKRAKEN_BIN,
            seqkit=SEQKIT_BIN,
            scripts=SCRIPTS, 
            db=DB,
            taxonomy=TAXONONMY,
            x=x,
            f=f,
            reads=reads,
            prefix=prefix,
            threads=threads,
        )
    )

    return task, os.path.join(work_dir, "%s.virus_reads.%s" % (prefix, f))


def create_flye_task(reads, prefix, data_type, size, threads, job_type, work_dir, out_dir):

    option = OrderedDict()
    option["flye"] = {
        "version": get_version(SOFTWARE_VERSION["flye"]),
        "option:": DATATYPE[data_type]["flye"]
    }
    memory = threads*4
    if memory >= 100:
        memory = 50

    task = Task(
        id="flye_%s" % prefix,
        work_dir=work_dir,
        type=job_type,
        option="-l vf=%dG -pe smp %s %s" % (memory, threads, QUEUE),
        script="""
export PATH={flye}:$PATH
flye {x} {reads} --threads {threads} --out-dir {prefix}_flye --genome-size {size}
cp {prefix}_flye/assembly_info.txt {prefix}.flye_info.txt
cp {prefix}_flye/assembly.fasta {prefix}.assembly.fasta
cp {prefix}_flye/assembly_graph.gfa {prefix}.flye.gfa
#rm -rf {prefix}_flye
""".format(flye=FLYE_BIN,
            prefix=prefix,
            reads=reads,
            threads=threads,
            size=size,
            x=DATATYPE[data_type]["flye"],
            out_dir=out_dir)
    )

    return task, option, os.path.join(work_dir, "%s.assembly.fasta" % prefix)


def run_asm_tgs_virus(reads, data_type, prefix, job_type, work_dir, out_dir,
    threads=10, size="10m", concurrent=10, refresh=30):

    reads = check_path(reads)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    work_dict = {
        "data": "00_data",
        "asm": "01_assemble",
    }
    for k, v in work_dict.items():
        work_dict[k] = mkdir(os.path.join(work_dir, v))
    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    dag = DAG("asm_virus")
    get_virus_task, reads = create_get_tgs_virus_task(
        reads=reads,
        prefix=prefix,
        job_type=job_type,
        work_dir=work_dict["data"],
        out_dir=out_dir,
        threads=10
    )
    dag.add_task(get_virus_task)

    flye_task, option, genome = create_flye_task(
        reads=reads,
        prefix=prefix,
        data_type=data_type,
        size=size,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["asm"],
        out_dir=out_dir
    )
    options["software"] = option
    dag.add_task(flye_task)
    flye_task.set_upstream(get_virus_task)

    do_dag(dag, concurrent, refresh)

    return options


def asm_tgs_virus(args):

    options = run_asm_tgs_virus(
         reads=args.reads,
         data_type=args.data_type,
         prefix=args.prefix,
         job_type=args.job_type,
         work_dir=args.work_dir,
         out_dir=args.out_dir,
         threads=args.threads,
         size=args.size,
         concurrent=args.concurrent,
         refresh=args.refresh
    )
    with open(os.path.join(args.out_dir, "asm_tgs_virus.json"), "w") as fh:
         json.dump(options, fh, indent=2)

    return 0


def add_help_args(parser):

    parser.add_argument("reads", metavar="FILE", type=str,
        help="Input third-generation sequencing data.")
    parser.add_argument("-dt", "--data_type", choices=["ont", "clr", "hifi"], default="ont",
        help="Input data type(ont, clr, hifi), default=ont.")
    parser.add_argument("-p", "--prefix", metavar="STR", type=str, default="A1",
        help="Input used in pipeline, default=A1")
    parser.add_argument("-s", "--size", metavar="STR", type=str, default="10m",
        help="Enter genome size, default=10m")
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
    asm_tgs_virus: Assembling Metaviral Genomes Using Three Generations of Data.

attention:
    asm_tgs_virus.py 

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_help_args(parser)
    args = parser.parse_args()
    asm_tgs_virus(args)


if __name__ == "__main__":
    main()
