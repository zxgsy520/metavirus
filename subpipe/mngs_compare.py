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
from ngsmetavirus.common import check_path, check_paths, read_tsv, mkdir, get_version
from dagflow import DAG, Task, ParallelTask, do_dag


LOG = logging.getLogger(__name__)
__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []

LEFSE_BIN = "/Work/user/liyilin/software/conda_v3/miniconda3/share/lefse-1.1.1-0/"
PYTHON_BIN = "/Work/user/liyilin/software/conda_v3/miniconda3/bin/"
TXADIC = {
"p": "phylum",
"c": "class",
"o": "order",
"f": "family",
"g": "genus",
"s": "species"
}


def create_compare_group(vslist, group, work_dir):

    data = {}

    for line in read_tsv(group):
        if line[-1] not in data:
            data[line[-1]] = []
        data[line[-1]].append(line)

    r = {}

    for line in read_tsv(vslist):
        vsid = "%s--%s" % (line[0], line[1])
        ngroup = os.path.join(work_dir, "%s_group.list" % vsid)
        r[vsid] = ngroup
        fo = open(ngroup, "w")
        for i in line:
            for j in data[i]:
                fo.write("%s\t%s\n" % (j[0], j[1]))
        fo.close()

    return r


def create_merge_mpa2otu_task(mpas, job_type, work_dir="", out_dir="",
                              kingdom="Bacteria", model="uniform"):
    task = Task(
        id="merge_mpa2otu",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
python {script}/merge_mpa2otu.py {mpas} --model original >merge_mpa.tsv
python {script}/get_relative_abundance.py merge_mpa.tsv --level p --kingdom {kingdom} \\
  --model {model} >{kingdom}.abundance_phylum.xls
python {script}/get_relative_abundance.py merge_mpa.tsv --level c --kingdom {kingdom} \\
  --model {model} >{kingdom}.abundance_class.xls
python {script}/get_relative_abundance.py merge_mpa.tsv --level o --kingdom {kingdom} \\
  --model {model} >{kingdom}.abundance_order.xls
python {script}/get_relative_abundance.py merge_mpa.tsv --level f --kingdom {kingdom} \\
  --model {model} >{kingdom}.abundance_family.xls
python {script}/get_relative_abundance.py merge_mpa.tsv --level g --kingdom {kingdom} \\
  --model {model} >{kingdom}.abundance_genus.xls
python {script}/get_relative_abundance.py merge_mpa.tsv --level s --kingdom {kingdom} \\
  --model {model} >{kingdom}.abundance_species.xls
cp merge_mpa.tsv {kingdom}.abundance_*.xls {out_dir}
""".format(script=SCRIPTS,
            mpas=" ".join(mpas),
            kingdom=kingdom,
            model=model,
            out_dir=out_dir,
        )
    )

    abunds = {
    "p": os.path.join(out_dir, "%s.abundance_phylum.xls" % kingdom),
    "c": os.path.join(out_dir, "%s.abundance_class.xls" % kingdom),
    "o": os.path.join(out_dir, "%s.abundance_order.xls" % kingdom),
    "f": os.path.join(out_dir, "%s.abundance_family.xls" % kingdom),
    "g": os.path.join(out_dir, "%s.abundance_genus.xls" % kingdom),
    "s": os.path.join(out_dir, "%s.abundance_species.xls" % kingdom),
    }

    return task, abunds


def create_barplot_tasks(abunds, group, job_type, work_dir="", out_dir="", top=10):

    sorts = []
    abundans = []
    for i in abunds:
        abundans.append(abunds[i])
        sorts.append(TXADIC[i])
    tasks = ParallelTask(
        id="barplot",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
export PATH={python}:{rbin}:$PATH
python {script}/cut_relative_abundance.py {{abunds}} --top {top} >abundance_{{sorts}}_top{top}.tsv
Rscript {script}/plot_abundance_bar.R abundance_{{sorts}}_top{top}.tsv {{sorts}} abundance_top{top}
cp abundance_top{top}_{{sorts}}.p* {out_dir}
if [ -n "{group}" ]; then 
    python {script}/cut_relative_abundance.py {{abunds}} --group {group} --top {top} >group.abundance_{{sorts}}_top{top}.tsv
    Rscript {script}/plot_abundance_bar.R group.abundance_{{sorts}}_top{top}.tsv {{sorts}} group.abundance_top{top}
    cp group.abundance_top{top}_{{sorts}}.p* {out_dir} 
fi

""".format(python=PYTHON_BIN,
           rbin=R_BIN,
           script=SCRIPTS,
           top=top,
           group=group,
           out_dir=out_dir),
        abunds=abundans,
        sorts=sorts,
    )

    return tasks


def create_lefse_tasks(abunds, group, job_type, work_dir="", out_dir="", suffix=""):

    sorts = []
    abundans = []
    for i in abunds:
        abundans.append(abunds[i])
        sorts.append(TXADIC[i])

    tasks = ParallelTask(
        id="lefse%s" % suffix,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
export PATH={lefse}:{python}:$PATH
python {script}/abundance2lefse.py {{abunds}} -g {group} --cut >{{sorts}}_lefse.tsv
python {lefse}/format_input.py {{sorts}}_lefse.tsv {{sorts}}_lefse.in -c 1 -o 1000000
python {lefse}/run_lefse.py {{sorts}}_lefse.in {{sorts}}_lefse.res -l 2.0
python {lefse}/plot_res.py {{sorts}}_lefse.res {{sorts}}_lefse.png --dpi 500 --format png --width 20
python {lefse}/plot_res.py {{sorts}}_lefse.res {{sorts}}_lefse.pdf --format pdf --width 20
python {script}/filter_lefse.py {{sorts}}_lefse.res --display 1 >{{sorts}}_lefse_filter.tsv
cp {{sorts}}_lefse.png {{sorts}}_lefse.pdf {{sorts}}_lefse.res {{sorts}}_lefse_filter.tsv {out_dir}

if [ "{{sorts}}" = "species" ] ; then
  python {script}/abundance2lefse.py {{abunds}} -g {group} >{{sorts}}_lefse.tsv
  python {lefse}/format_input.py {{sorts}}_lefse.tsv {{sorts}}_lefse.in -c 1 -o 1000000
  python {lefse}/run_lefse.py {{sorts}}_lefse.in {{sorts}}_lefse.res -l 2.0
  python {script}/filter_lefse.py {{sorts}}_lefse.res >{{sorts}}.tsv
  python {lefse}/plot_cladogram.py --dpi 500 --format png {{sorts}}.tsv lefse_cladogram.png
  python {lefse}/plot_cladogram.py --dpi 500 --format pdf {{sorts}}.tsv lefse_cladogram.pdf
  cp lefse_cladogram.png lefse_cladogram.pdf {out_dir}
fi

rm {{sorts}}_lefse.in
""".format(lefse=LEFSE_BIN,
           python=PYTHON_BIN,
           script=SCRIPTS,
           group=group,
           out_dir=out_dir),
        abunds=abundans,
        sorts=sorts,
    )

    return tasks


def run_compare(mpas, group, vslist, job_type, work_dir="", out_dir="",
                kingdom="Bacteria", model="uniform", concurrent=10, refresh=30):

    mpas = check_paths(mpas)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    if group:
        group = check_path(group)
    if vslist:
        vslist = check_path(vslist)
        groups = create_compare_group(vslist, group, work_dir)
    else:
        groups = {}

    work_dict = {
        "barplot": "barplot",
        "lefse": "lefse",
    }
    for k, v in work_dict.items():
        if not group and k=="lefse":
            continue
        mkdir(os.path.join(work_dir, v))
        mkdir(os.path.join(out_dir, v))

    dag = DAG("run_compare")
    mpa2otu_task, abunds = create_merge_mpa2otu_task(
        mpas=mpas,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir,
        kingdom=kingdom,
        model=model
    )
    dag.add_task(mpa2otu_task)

    barplot_tasks = create_barplot_tasks(
        abunds=abunds,
        group=group,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["barplot"]),
        out_dir=os.path.join(out_dir, work_dict["barplot"]),
        top=10
    )
    dag.add_task(*barplot_tasks)
    mpa2otu_task.set_downstream(*barplot_tasks)

    if group:
        work = os.path.join(work_dir, work_dict["lefse"])
        out = os.path.join(out_dir, work_dict["lefse"])
        lefse_tasks = create_lefse_tasks(
            abunds=abunds,
            group=group,
            job_type=job_type,
            work_dir=work,
            out_dir=out
        )
        dag.add_task(*lefse_tasks)
        mpa2otu_task.set_downstream(*lefse_tasks)

    for i in groups:
        temp_work = mkdir(os.path.join(work, i))
        temp_out =  mkdir(os.path.join(out, i))
        lefse_tasks = create_lefse_tasks(
            abunds=abunds,
            group=groups[i],
            job_type=job_type,
            work_dir=temp_work,
            out_dir=temp_out,
            suffix=i,
        )
        dag.add_task(*lefse_tasks)
        mpa2otu_task.set_downstream(*lefse_tasks)

    do_dag(dag, concurrent, refresh)

    return 0


def mngs_compare(args):

    work_dir = mkdir(os.path.join(args.work_dir, "Bacteria"))
    out_dir = mkdir(os.path.join(args.out_dir, "Bacteria"))
    run_compare(
         mpas=args.input,
         group=args.group,
         vslist=args.vslist,
         kingdom="Bacteria",
         model=args.model,
         job_type=args.job_type,
         work_dir=work_dir,
         out_dir=out_dir,
         concurrent=args.concurrent,
         refresh=args.refresh
    )

    work_dir = mkdir(os.path.join(args.work_dir, "Viruses"))
    out_dir = mkdir(os.path.join(args.out_dir, "Viruses"))
    run_compare(
         mpas=args.input,
         group=args.group,
         vslist=args.vslist,
         kingdom="Viruses",
         model=args.model,
         job_type=args.job_type,
         work_dir=work_dir,
         out_dir=out_dir,
         concurrent=args.concurrent,
         refresh=args.refresh
    )

    return 0


def add_compare_args(parser):

    parser.add_argument("input", nargs='+', metavar="FILE", type=str,
        help="Input the abundance statistics result file of each sample, *.kreport2mpa.report")
    parser.add_argument("-g", "--group", metavar="FILE", type=str, default="",
        help="Input sample grouping table, group.list.")
    parser.add_argument("-vs", "--vslist", metavar="FILE", type=str, default="",
        help="Input the grouped files to compare, vs.list.")
    parser.add_argument("-m", "--model", metavar="STR", type=str,  default="uniform",
        choices=["mean", "uniform", "original"],
        help="Abundance conversion method. choices=[mean, uniform, original], default=uniform.")
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
    mngs_compare: Metagenome comparative analysis

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_compare_args(parser)
    args = parser.parse_args()
    mngs_compare(args)


if __name__ == "__main__":
    main()
