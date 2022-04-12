#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import shutil
import logging
import argparse

from ngsmetavirus.config import *
from ngsmetavirus.mngs_qc import run_mngs_qc
from ngsmetavirus.mngs_tax import run_mngs_tax
from ngsmetavirus.mngs_asm import run_mngs_asm
from ngsmetavirus.mngs_ann import run_mngs_ann
from ngsmetavirus.parser import add_mngs_all_args
from dagflow import DAG, Task, do_dag
from ngsmetavirus.common import check_path, mkdir


LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.1.0"


def run_report(sample, project, id, platform, out_dir, work_dir, job_type):

    dag = DAG("run_report")
    task = Task(
        id="report",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 -q all.q",
        script="""
python {script}/single_mngs_report.py --project {project} --id {id} \\
  --sequencer {platform} --sample {sample} \\
  --data {result}/01_data/{sample}.stat_qc.tsv \\
  --tax {result}/02_reads_tax/{sample}.stat_read_tax.tsv \\
  --contig {result}/03_assembly/{sample}.stat_genome.tsv \\
  --gene {result}/04_annotation/{sample}.stat_gene.tsv \\
  --function {result}/04_annotation/function/{sample}.function_summary.tsv \\
  --advanal {result}/04_annotation/function/{sample}.stat_cazy_phi.tsv \\
  --card_vfdb {result}/04_annotation/card_vfdb/{sample}.vfdb_summary.tsv {result}/04_annotation/card_vfdb/{sample}.card_summary.tsv \\
  --krona {result}/02_reads_tax/{sample}.taxonomy.html \\
  --protein {result}/04_annotation/{sample}.protein_length.png \\
  --species  {result}/04_annotation/function/{sample}.species.png  \\
  --cog  {result}/04_annotation/function/{sample}.COG.png \\
  --kegg  {result}/04_annotation/function/{sample}.KEGG.png \\
  --go  {result}/04_annotation/function/{sample}.WEGO.png \\
  --cazy  {result}/04_annotation/function/{sample}.cazy.png \\
  --phi  {result}/04_annotation/function/{sample}.phi_classify.png \\
  --html {template}/single_mngs_html \\
  --out {out_dir}
""".format(script=SCRIPTS,
            template=TEMPLATES,
            sample=sample,
            project=project,
            id=id,
            platform=platform,
            result=out_dir,
            out_dir=out_dir
        )
    )
    dag.add_task(task)
    
    md5sum_task = Task(
        id="md5sum",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 -q all.q",
        script="""
cd {outdir}
zip -r {id}_{sample}.report.zip static/ images/ report.html 
cd {outdir}
{bin}/makemd5 0* >result.md5
rm -rf static/ images/ report.html
""".format(bin=BIN,
            sample=sample,
            id=id,
            outdir=out_dir
        )
    )
    dag.add_task(md5sum_task)
    md5sum_task.set_upstream(task)

    do_dag(dag, 1, 5)

    return 0


def run_mngs_all(prefix, read1, read2, ref, nohost, dtype, trim, memory,
                 thread, job_type, concurrent, refresh, work_dir, out_dir,
                 qvalue=20, atype="metagenome", project="", pid=""):


    read1 = check_path(read1)
    read2 = check_path(read2)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    work_dict = {
        "data": "01_data",
        "tax": "02_reads_tax",
        "asm": "03_assembly",
        "ann": "04_annotation",
    }

    for k, v in work_dict.items():
        mkdir(os.path.join(work_dir, v))
        mkdir(os.path.join(out_dir, v))

    options, clean1, clean2 = run_mngs_qc(
        prefix=prefix,
        read1=read1,
        read2=read2,
        ref=ref,
        nohost=nohost, 
        dtype=dtype,
        trim=trim,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=mkdir(os.path.join(work_dir, work_dict["data"])),
        out_dir=mkdir(os.path.join(out_dir, work_dict["data"])),
        qvalue=qvalue,
        atype=atype
    )
    
    noptions = run_mngs_tax(
        prefix=prefix,
        read1=read1,
        read2=read2,
        thread=int(thread*1.5),
        job_type=job_type,
        work_dir=mkdir(os.path.join(work_dir, work_dict["tax"])),
        out_dir=mkdir(os.path.join(out_dir, work_dict["tax"])),
        concurrent=concurrent,
        refresh=refresh
    )
    options["software"].update(noptions["software"])    

    noptions, genome = run_mngs_asm(
        prefix=prefix,
        read1=clean1,
        read2=clean2,
        atype=atype,
        thread=thread,
        memory=memory,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=mkdir(os.path.join(work_dir, work_dict["asm"])),
        out_dir=mkdir(os.path.join(out_dir, work_dict["asm"])),
    )
    options["software"].update(noptions["software"])
    
    run_mngs_ann(genome=genome,
        prefix=prefix,
        threads=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=mkdir(os.path.join(work_dir, work_dict["ann"])),
        out_dir=mkdir(os.path.join(out_dir, work_dict["ann"]))
    )

    with open(os.path.join(out_dir, "ngsmeta.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    run_report(sample=prefix,
        project=project,
        id=pid,
        platform=dtype,
        out_dir=out_dir,
        work_dir=work_dir,
        job_type=job_type
    )

    return options


def mngs_all(args):

    run_mngs_all(
        prefix=args.prefix,
        read1=args.read1,
        read2=args.read2,
        ref=args.reference,
        nohost=args.nohost,
        dtype=args.dtype,
        trim=args.trim,
        memory=args.memory,
        thread=args.thread,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        qvalue=args.qvalue,
        atype=args.atype,
        project=args.project,
        pid=args.id
    )


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

    parser = add_mngs_all_args(parser)
    args = parser.parse_args()
    mngs_all(args)


if __name__ == "__main__":
    main()
