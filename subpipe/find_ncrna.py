#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
import os.path
from collections import OrderedDict

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../"))
from ngsmetavirus.common import get_version, mkdir, cd
from ngsmetavirus.config import *
from seqkit.split import seq_split
from dagflow import Task, ParallelTask, DAG, do_dag


LOG = logging.getLogger(__name__)
__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


QUEUE = "-q all.q,s01"
TRNASCAN_BIN = "/Work/pipeline/software/Base/trnascan-se/v2.0.9/bin/"
BARRNAP_BIN = "/Work/pipeline/software/Base/barrnap/v0.9/bin/"
INFERNAL_BIN = "/Work/pipeline/software/Base/infernal/v1.1.4/bin/"
RFAM = "/Work/database/Rfam/"

SOFTWARE_VERSION = {
    "tRNAscan": {
        "GETVER": "%s/tRNAscan-SE 2>&1|grep -i '^tRNAscan-SE'" % TRNASCAN_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "2.0"
    },
    "barrnap": {
        "GETVER": "%s/barrnap -v 2>&1" % BARRNAP_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "0.9"
    },
    "infernal": {
        "GETVER": "%s/cmscan -h 2>&1| grep 'INFERNAL'" % INFERNAL_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "1.1.2"
    },
}

def create_trnascan_task(genomes, prefix, thread, kingdom="Bacteria",
                         job_type="sge", work_dir="", out_dir=""):

    if kingdom == "Bacteria":
        trna_mode = "-B -I"
    elif kingdom == "Archaea":
        trna_mode = "-A -I"
    elif kingdom == "Eukaryote":
        trna_mode = "-E -I"
    else:
        raise Exception()

    option = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }
    version = get_version(SOFTWARE_VERSION["tRNAscan"])
    option["software"]["tRNAscan"] = {
        "version": version,
        "option": "%s -m lsu,ssu,tsu" % trna_mode
    }
    prefixs  = [os.path.basename(i) for i in genomes]
    tasks = ParallelTask(
        id="tRNA",
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""\
export PATH={tRNAscan}:$PATH
# delete files created by previous run if exists, avoid tRNAscan overwriting check
rm -rf {{prefix}}.tRNA.*

time tRNAscan-SE --thread {thread} {trna_mode} \\
-o {{prefix}}.tRNA.tsv -f {{prefix}}.tRNA.structure.txt \\
{{genome}}
""".format(tRNAscan=TRNASCAN_BIN,
           thread=thread,
           trna_mode=trna_mode
           ),
        genome=genomes,
        prefix=prefixs
    )

    join_task = Task(
        id="merge_tRNA",
        work_dir=work_dir,
        type="local",
        option="-pe smp 1 %s" % QUEUE,
        script="""\
echo "#VERSION: tRNAscan-SE:{version}\n#ARGS:{trna_mode}\n" >{prefix}.tRNA.tsv
cat */*.tRNA.tsv >> {prefix}.tRNA.tsv
cat */*.tRNA.structure.txt > {prefix}.tRNA.structure.txt

cp {prefix}.tRNA.tsv {prefix}.tRNA.structure.txt {out_dir}
    """.format(prefix=prefix,
               version=version,
               trna_mode=trna_mode,
               out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task, option, os.path.join(work_dir, "%s.tRNA.tsv" % prefix)



def create_barrnap_task(genomes, prefix, thread, kingdom="Bacteria",
                         job_type="sge", work_dir="", out_dir=""):

    if kingdom == "Bacteria":
        rrna_mode ="--kingdom bac"
    elif kingdom == "Archaea":
        rrna_mode ="--kingdom arc"
    elif kingdom == "Eukaryote":
        rrna_mode ="--kingdom euk"
    else:
        raise Exception()

    option = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }
    version = get_version(SOFTWARE_VERSION["barrnap"])
    option["software"]["barrnap"] = {
        "version": version,
        "option": "%s --quiet" % rrna_mode
    }
    prefixs = [os.path.basename(i) for i in genomes]

    tasks = ParallelTask(
        id="rRNA",
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""\
export PATH={barrnap}:$PATH
time barrnap {rrna_mode} --threads {thread} --quiet {{genome}} >{{prefix}}.gff
""".format(barrnap=BARRNAP_BIN,
           rrna_mode=rrna_mode,
           thread=thread,
           ),
        genome=genomes,
        prefix=prefixs
    )

    join_task = Task(
        id="merge_rRNA",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 -V %s" % QUEUE,
        script="""\
echo "#VERSION: Barrnap:{version}\n#ARGS:{rrna_mode}\n" > {prefix}.rRNA.tsv
cat */*.gff >> {prefix}.rRNA.tsv

cp {prefix}.rRNA.tsv {out_dir}
""".format(prefix=prefix,
           version=version,
           rrna_mode=rrna_mode,
           out_dir=out_dir
           )
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task, option, os.path.join(work_dir, "%s.rRNA.tsv" % prefix)


def create_rfam_task(genomes, prefix, thread, job_type="sge",
                     work_dir=".", out_dir=""):

    prefixs = [os.path.basename(i) for i in genomes]
    option = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }
    version = get_version(SOFTWARE_VERSION["infernal"])
    option["software"]["infernal"] = {
        "version": version,
        "option": "--cut_ga --rfam --nohmmonly"
    }

    tasks = ParallelTask(
        id="rfam" ,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""\
export PATH={infernal}:$PATH
# for rfam 14.0
# from http://rfam.readthedocs.io/en/latest/genome-annotaion.html

time cmscan --cpu {thread} --cut_ga --rfam --nohmmonly --fmt 2 \\
--clanin {rfam}/Rfam.clanin \\
--tblout {{prefix}}.rfam.tsv \\
{rfam}/Rfam.cm {{genome}} > {{prefix}}.cmscan
""".format(rfam=RFAM,
           infernal=INFERNAL_BIN,
           thread=thread,
           ),
        genome=genomes,
        prefix=prefixs,
    )

    join_task = Task(
        id="merge_rfam",
        work_dir=work_dir,
        type="local",
        option="-pe smp 1 %s" % QUEUE,
        script="""\
echo "#VERSION: INFERNAL:{version}\n#ARGS:--cut_ga --rfam --nohmmonly \n" > {prefix}.rfam.tsv
cat */*.rfam.tsv >> {prefix}.rfam.tsv
cp {prefix}.rfam.tsv {out_dir}
""".format(prefix=prefix,
           version=version,
           out_dir=out_dir
           )
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task, option, os.path.join(work_dir, "%s.rfam.tsv" % prefix)


def create_rna_join(prefix, rrna, trna, rfam, work_dir, out_dir):

    task = Task(
        id="rna_join",
        work_dir=work_dir,
        type="local",
        option="-pe smp 1 %s" % QUEUE,
        script="""
cd {out_dir}
{script}/process_ncRNA.py \\
--rRNA {rrna} \\
--tRNA {trna} \\
--rfam {rfam} --rfamily {rfam_bd}/family.txt > {prefix}.ncRNA.gff3
#{script}/stat_gff.py {prefix}.ncRNA.gff3 -t tRNA rRNA ncRNA regulatory > {prefix}.ncRNA_stat.tsv
""".format(script=SCRIPTS,
            prefix=prefix,
            rrna=rrna,
            trna=trna,
            rfam=rfam,
            rfam_bd=RFAM,
            out_dir=out_dir)
    )

    return task


def run_ncrna_annotation(genome, prefix, kingdom, thread, job_type, work_dir, out_dir,
                         concurrent=10, refresh=30):

    work_dict = {
        'split': "00_data",
        'trna': "01_tRNAscan",
        'rrna': "02_barrnap",
        'rfam': "03_Rfam"
    }

    work_dir = mkdir(os.path.abspath(work_dir))
    out_dir = mkdir(os.path.abspath(out_dir))
    genome = os.path.abspath(genome)

    for k, v in work_dict.items():
        work_dict[k] = mkdir(os.path.join(work_dir, v))

    genomes = seq_split([genome], mode="length", num=10000000, output_dir=work_dict["split"])

    dag = DAG("ncRNA-Annotation")

    trna_tasks, trna_join, option, trna = create_trnascan_task(
        genomes=genomes,
        prefix=prefix,
        thread=thread,
        job_type=job_type,
        work_dir=work_dict["trna"],
        out_dir=out_dir
    )

    options = option
    dag.add_task(*trna_tasks)
    dag.add_task(trna_join)

    rrna_tasks, rrna_join, option, rrna = create_barrnap_task(
        genomes=genomes,
        prefix=prefix,
        thread=thread,
        job_type=job_type,
        work_dir=work_dict["rrna"],
        out_dir=out_dir
    )
    options["software"].update(option["software"])
    dag.add_task(*rrna_tasks)
    dag.add_task(rrna_join)


    rfam_tasks, rfam_join, option, rfam = create_rfam_task(
        genomes=genomes,
        prefix=prefix,
        thread=thread,
        job_type=job_type,
        work_dir=work_dict["rfam"],
        out_dir=out_dir
    )
    options["software"].update(option["software"])
    dag.add_task(*rfam_tasks)
    dag.add_task(rfam_join)


    rna_join = create_rna_join(
        prefix=prefix,
        rrna=rrna,
        trna=trna,
        rfam=rfam,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(rna_join)
    trna_join.set_downstream(rna_join)
    rrna_join.set_downstream(rna_join)
    rfam_join.set_downstream(rna_join)

    do_dag(dag, concurrent_tasks=concurrent, refresh_time=refresh)

    return options


def find_ncrna(args):

    options = run_ncrna_annotation(
        genome=args.genome,
        prefix=args.prefix,
        kingdom=args.kingdom,
        thread=args.thread,
        job_type=args.job_type,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        concurrent=args.concurrent,
        refresh=args.refresh)
    with open(os.path.join(args.out_dir, "find_ncrna.json"), "w") as fh:
         json.dump(options, fh, indent=2)

    return 0


def add_help_args(parser):

    parser.add_argument("genome", metavar="FILE",
        help="Genome in fasta file")
    parser.add_argument("-p", "--prefix", metavar="STR", required=True,
        help="Name used in pipeline")
    parser.add_argument("--kingdom", choices=["Archaea", "Bacteria", "Mitochondria", "Viruses", "Eukaryote"],
        help="Kingdom of the sample (default: Bacteria)", default="Bacteria"),
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
    description:

    version: %s
    contact:  %s <%s>\
        """ % (__version__, " ".join(__author__), __email__))

    parser = add_help_args(parser)
    args = parser.parse_args()
    find_ncrna(args)


if __name__ == "__main__":
    main()
