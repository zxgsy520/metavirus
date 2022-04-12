#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import logging
import argparse

from collections import OrderedDict
from ngsmetavirus.config import *
from ngsmetavirus.parser import add_mngs_ann_args
from dagflow import DAG, Task, do_dag
from ngsmetavirus.common import check_path, mkdir, get_version


LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.0.0"


def create_gmhmmp_task(genome, prefix, job_type, work_dir, out_dir):

    task = Task(
        id="gmhmmp",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
#{biotool}/biotool sort_genome {genome} --minlen 250 >{prefix}.genome.fasta
ln -sf {genome} {prefix}.genome.fasta
#cp {gmhmmp}/gm_key_64 ~/.gm_key
{gmhmmp}/gmhmmp -a -d -f G -m {gmhmmp}/MetaGeneMark_v1.mod -o {prefix}.mgm {prefix}.genome.fasta
python {script}/mgm2gff.py {prefix}.mgm --prefix {prefix} >{prefix}.genome.gff3
python {script}/plot_cds2pep.py {prefix}.gene.fasta --prefix {prefix}
{biotool}/biotool stats {prefix}.gene.fasta >{prefix}.stat_gene.tsv
cp {prefix}.genome.fasta {prefix}.protein_length.p* {prefix}.stat_gene.tsv {out_dir}
cp {prefix}.genome.gff3 {prefix}.gene.fasta {prefix}.protein.fasta {prefix}.genome.fasta {out_dir}
""".format(biotool=BIN,
            gmhmmp=GMHMMP_BIN,
            script=SCRIPTS,
            genome=genome,
            prefix=prefix,
           out_dir=out_dir)
    )
    gene = os.path.join(out_dir, "%s.gene.fasta" % prefix)
    protein = os.path.join(out_dir, "%s.protein.fasta" % prefix)
    gff = os.path.join(out_dir, "%s.genome.gff3" % prefix)

    return task, gene, protein, gff


def create_abricate_task(gene, prefix, database, threads, job_type,
                         work_dir, out_dir):

    task = Task(
        id="%s_%s" % (database, prefix),
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={abricate}:$PATH
abricate {gene} --threads {threads} \\
  --mincov 50 --minid 75 \\
  --db {database}> {prefix}.{database}_old.tsv
{script}/abricate_summary.py {prefix}.{database}_old.tsv \\
  -o {prefix}.{database}_summary.tsv >{prefix}.{database}.tsv
cp {prefix}.{database}.tsv {prefix}.{database}_summary.tsv {out_dir}
""".format(id=id,
        abricate=ABRICATE_BIN,
        datadir=ABRICATE_DB,
        script=SCRIPTS,
        gene=gene,
        prefix=prefix,
        database=database,
        threads=threads,
        out_dir=out_dir)
    )

    return task, os.path.join(work_dir, '%s.%s_summary.tsv' % (prefix, database))


def create_abricate_tasks(gene, prefix, database, threads,
                          job_type, work_dir, out_dir):

    option = {}
    tasks =[]
    abricates = []
    for i in database.strip().split(';'):
        option[i] = {
            "version": get_version(SOFTWARE_VERSION["abricate"]),
            "option": "abricate  --mincov 50 --minid 75  --db %s" % i
        }
        temp_work = mkdir(os.path.join(work_dir, i))
        task, abricate = create_abricate_task(
            gene=gene,
            prefix=prefix,
            database=i,
            threads=threads,
            job_type=job_type,
            work_dir=temp_work,
            out_dir=out_dir)

        tasks.append(task)
        abricates.append(abricate)

    return tasks, abricates, option


def create_function_task(protein, prefix, threads, job_type,
                         work_dir, out_dir):

    task = Task(
        id="gene_function",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
python {root}/subpipe/gene_function.py {protein} --prefix {prefix}\\
  --kingdom meta --threads {threads} --job_type {job_type} --concurrent 20 \\
  --refresh 15 --work_dir {work_dir} --out_dir {out_dir}
""".format(root=ROOT,
            protein=protein,
            prefix=prefix,
            threads=threads,
            job_type=job_type,
            work_dir=work_dir,
            out_dir=out_dir)
    )

    return task


def run_mngs_ann(genome, prefix, threads, job_type, concurrent, refresh,
                 work_dir="", out_dir=""):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    genome = check_path(genome)
    
    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    work_dict = {
        "card_vfdb": "card_vfdb",
        "function": "function",
    }
    for k, v in work_dict.items():
        mkdir(os.path.join(work_dir, v))
        mkdir(os.path.join(out_dir, v))

    dag = DAG("mngs_ann")
    gmhmmp_task, gene, protein, gff = create_gmhmmp_task(
        genome=genome,
        prefix=prefix,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(gmhmmp_task)

    tasks, abricates, option = create_abricate_tasks(
        gene=gene,
        prefix=prefix,
        database="card;vfdb",
        threads=threads,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["card_vfdb"]),
        out_dir=os.path.join(out_dir, work_dict["card_vfdb"])
    )
    options["software"].update(option)
    for i in tasks:
        dag.add_task(i)
        i.set_upstream(gmhmmp_task)

    function_task = create_function_task(
        protein=protein,
        prefix=prefix,
        threads=threads,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["function"]),
        out_dir=os.path.join(out_dir, work_dict["function"])
    )
    dag.add_task(function_task)
    function_task.set_upstream(gmhmmp_task)

    do_dag(dag, concurrent_tasks=concurrent, refresh_time=refresh)

    return 0


def mngs_ann(args):

    run_mngs_ann(genome=args.genome,
    prefix=args.prefix,
    threads=args.threads,
    job_type=args.job_type,
    concurrent=args.concurrent,
    refresh=args.refresh,
    work_dir=args.work_dir,
    out_dir=args.out_dir)

    return 0


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
name:
    mngs_ann.py : Metagenome structural and functional annotation
attention:
    mngs_ann.py genome.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_mngs_ann_args(parser).parse_args()

    mngs_ann(args)


if __name__ == "__main__":

    main()
