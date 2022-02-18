#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import logging
import argparse

from collections import OrderedDict
from ngsmetavirus.config import *
from ngsmetavirus.parser import add_mngs_qc_args
from dagflow import DAG, Task, do_dag
from ngsmetavirus.common import check_path, mkdir, read_tsv, read_files, get_version


LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.0.1"


def quality_control_task(read1, read2, prefix, trim, thread, job_type, work_dir, out_dir, qvalue=20):

    option = OrderedDict()
    option["fastp"] = {
        "version": get_version(SOFTWARE_VERSION["fastp"]),
        "option": "-n 0 -f {trim} -F {trim} -t {trim} -T {trim} -q {qvalue}".format(trim=trim, qvalue=qvalue)
    }

    task = Task(
        id="ngs_qc",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={fastp}:$PATH
fastp -i {read1} -I {read2} \\
    -o {prefix}.clean.r1.fq.gz -O {prefix}.clean.r2.fq.gz \\
    -w {thread} -n 0 -f {trim} -F {trim} -t {trim} -T {trim} \\
    -q {qvalue} --json {prefix}_fastp.json
export PATH={python}:$PATH
{scripts}/stat_fastp.py {prefix}_fastp.json > {prefix}.stat_qc.tsv
cp {prefix}.stat_qc.tsv {out_dir}
""".format(fastp=FASTP_BIN,
            python=PYTHON_BIN,
            scripts=SCRIPTS,
            read1=read1,
            read2=read2,
            prefix=prefix,
            thread=thread,
            trim=trim,
            qvalue=qvalue,
            out_dir=out_dir
        )
    )
    stat_qc = os.path.join(out_dir, "%s.stat_qc.tsv" % prefix)
    clean1 = os.path.join(work_dir, "%s.clean.r1.fq.gz" % prefix)
    clean2 = os.path.join(work_dir, "%s.clean.r2.fq.gz" % prefix)

    return task, option, clean1, clean2, stat_qc


def contamination_eval_task(prefix, read1, kingdom, thread, job_type, work_dir, out_dir):

    option = OrderedDict()
    option["blastn"] = {
        "version": get_version(SOFTWARE_VERSION["blastn"]),
        "option": "default"
    }

    task = Task(
        id="cont_eval",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={blast}:{python}:$PATH
python {scripts}/fq2fa.py {read1} -n 10000 >{prefix}.clean.r1.fa
blastn -query {prefix}.clean.r1.fa -db {dbase} \\
  -outfmt "6 std staxid sskingdom staxids" -max_target_seqs 5 -num_threads {thread} \\
  -out {prefix}.m6
python {scripts}/obtain_taxonomy.py -i {prefix}.m6 -t {taxonomy} -n {prefix}
python {scripts}/stat_taxonomy.py -i {prefix}.species_annotation.txt -rn 10000 -n {prefix}
python {scripts}/plot_stat_species.py {prefix}.stat_species.tsv -p {prefix} >{prefix}.top10_species.tsv
cp {prefix}.species_classify.tsv {prefix}.stat_species.tsv {out_dir}
cp {prefix}.top10_species.tsv {prefix}.top10_species.png {prefix}.top10_species.pdf {out_dir}
cp {prefix}.species_annotation.txt {out_dir}
""".format(scripts=SCRIPTS,
            blast=BLAST_BIN,
            dbase=NT_TAXON[kingdom],
            python=PYTHON_BIN,
            taxonomy=TAXONOMY,
            prefix=prefix,
            read1=read1,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, option


def create_hisat2_tasks(prefix, read1, read2, ref, nohost, thread, job_type, work_dir, out_dir):

    option = OrderedDict()
    option["hisat2"] = {
        "version": get_version(SOFTWARE_VERSION["hisat2"]),
        "option": "--un-conc-gz"
    }
    option["samtools"] = {
        "version": get_version(SOFTWARE_VERSION["samtools"]),
        "option": "default"
    }

    if nohost:
        clean1 = os.path.join(out_dir, "%s.map.r1.fq.gz" % prefix)
        clean2 = os.path.join(out_dir, "%s.map.r2.fq.gz" % prefix)
        sname = "map"
    else:
        clean1 = os.path.join(out_dir, "%s.nohost.r1.fq.gz" % prefix)
        clean2 = os.path.join(out_dir, "%s.nohost.r2.fq.gz" % prefix)
        sname = "nohost"

    task = Task(
        id="hisat2",
        work_dir=work_dir,
        type=job_type,
        option="-V -pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={hisat}:{samtools}:$PATH
export PERL5LIB="":$PERL5LIB
hisat2 --threads {thread} -x {ref} \\
    -1 {read1} -2 {read2} \\
    --al-conc-gz host --un-conc-gz nohost 2>alignment_stat.txt | samtools view --threads {thread} -bS | \\
samtools sort --threads {thread} -m 4G -o {prefix}.sort.bam
mv host.1 {prefix}.map.r1.fq.gz
mv host.2 {prefix}.map.r2.fq.gz
mv nohost.1 {prefix}.nohost.r1.fq.gz
mv nohost.2 {prefix}.nohost.r2.fq.gz
{scripts}/stat_ngs.py {prefix}.{sname}.r*.fq.gz --name {prefix}>{prefix}.stat_expect.tsv
mv {prefix}.{sname}.r*.fq.gz {out_dir}
cp {prefix}.stat_expect.tsv {out_dir}
#rm {prefix}.sort.bam
""".format(hisat=HISAT2_BIN,
            samtools=SAMTOOLS_BIN,
            scripts=SCRIPTS,
            ref=ref,
            read1=read1,
            read2=read2,
            prefix=prefix,
            sname=sname,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, option, clean1, clean2


def create_minimap2_tasks(prefix, read1, read2, ref, nohost, dtype, thread, job_type, work_dir, out_dir):

    option = OrderedDict()
    option["minimap2"] = {
        "version": get_version(SOFTWARE_VERSION["minimap2"]),
        "option": "--secondary=no -x sr"
    }
    option["samtools"] = {
        "version": get_version(SOFTWARE_VERSION["samtools"]),
        "option": "default"
    }

    if nohost:
        clean1 = os.path.join(out_dir, "%s.map.r1.fq.gz" % prefix)
        clean2 = os.path.join(out_dir, "%s.map.r2.fq.gz" % prefix)
        sname = "map"
        temp = """
minimap2 -t {thread} -x sr {ref} \\
    {read1} {read2} |paf2mapid --dtype {dtype} >map.id
""".format(ref=ref,
            read1=read1,
            read2=read2,
            dtype=dtype,
            thread=thread
        )   
    else:
        clean1 = os.path.join(out_dir, "%s.nohost.r1.fq.gz" % prefix)
        clean2 = os.path.join(out_dir, "%s.nohost.r2.fq.gz" % prefix)
        sname = "nohost"
        temp = """
minimap2 -t {thread} -ax sr {ref} \\
    {read1} {read2} |samtools view |ngs_sam2unmapid --dtype {dtype} >nohost.id
""".format(ref=ref,
            read1=read1,
            read2=read2,
            dtype=dtype,
            thread=thread
        )

    task = Task(
        id="minimap2",
        work_dir=work_dir,
        type=job_type,
        option="-V -pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={minimap2}:{samtools}:{sam2ngs}:$PATH
{temp}
export PATH={seqkit}:$PATH
seqkit grep -f {sname}.id {read1} |gzip >{prefix}.{sname}.r1.fq.gz
seqkit grep -f {sname}.id {read2} |gzip >{prefix}.{sname}.r2.fq.gz
{scripts}/stat_ngs.py {prefix}.{sname}.r*.fq.gz --name {prefix}>{prefix}.stat_expect.tsv
mv {prefix}.{sname}.r*.fq.gz {out_dir}
cp {prefix}.stat_expect.tsv {out_dir}
""".format(minimap2=MINIMAP_BIN,
            samtools=SAMTOOLS_BIN,
            sam2ngs=SAM2NGS_BIN,
            seqkit=SEQKIT_BIN,
            scripts=SCRIPTS,
            temp=temp,
            sname=sname,
            prefix=prefix,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, option, clean1, clean2


def run_mngs_qc(prefix, read1, read2, ref, nohost, dtype, trim, thread, job_type,
                concurrent, refresh, work_dir, out_dir, qvalue=20, atype="metagenome"):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    read1 = check_path(read1)
    read2 = check_path(read2)

    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }
    dag = DAG("mngs_qc")

    qc_task, option, clean1, clean2, stat_qc = quality_control_task(
        read1=read1,
        read2=read2,
        prefix=prefix,
        trim=trim,
        thread=thread,
        job_type=job_type,
        work_dir=mkdir(os.path.join(work_dir, "01_qc")),
        out_dir=out_dir,
        qvalue=qvalue
    )
    options["software"].update(option)
    dag.add_task(qc_task)

#    cont_task, option = contamination_eval_task(
#        prefix=prefix,
#        read1=read1,
#        kingdom="meta",
#        thread=int(thread*1.5),
#        job_type=job_type,
#        work_dir=mkdir(os.path.join(work_dir, "01_qc")),
#        out_dir=out_dir
#    )
#    options["software"].update(option)
#    dag.add_task(cont_task)

    if not ref:
        pass
    elif atype in ["metaviral", "rnaviral"]:
        check_path("%s.3.ht2" % ref)
        hisat_task, option, clean1, clean2 = create_hisat2_tasks(
            prefix=prefix,
            read1=clean1,
            read2=clean2,
            ref=ref,
            nohost=nohost,
            thread=thread,
            job_type=job_type,
            work_dir= mkdir(os.path.join(work_dir, "02_rm_host")),
            out_dir=out_dir
        )
        options["software"].update(option)
        dag.add_task(hisat_task)
        hisat_task.set_upstream(qc_task)
    else:
        ref = check_path(ref)
        map_task, option, clean1, clean2 = create_minimap2_tasks(
            prefix=prefix,
            read1=clean1,
            read2=clean2,
            ref=ref,
            nohost=nohost,
            dtype=dtype,
            thread=thread,
            job_type=job_type,
            work_dir= mkdir(os.path.join(work_dir, "02_rm_host")),
            out_dir=out_dir
        )
        options["software"].update(option)
        dag.add_task(map_task)
        map_task.set_upstream(qc_task)

    do_dag(dag, concurrent, refresh)

    return options, clean1, clean2


def mngs_qc(args):

    options, clean1, clean2 = run_mngs_qc(
        prefix=args.prefix,
        read1=args.read1,
        read2=args.read2,
        ref=args.reference,
        nohost=args.nohost,
        dtype=args.dtype,
        trim=args.trim,
        thread=args.thread,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        qvalue=args.qvalue,
        atype=args.atype
    )

    with open(os.path.join(args.out_dir, "mngs_qc.json"), "w") as fh:
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

    parser = add_mngs_qc_args(parser)
    args = parser.parse_args()
    mngs_qc(args)


if __name__ == "__main__":
    main()
