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
__version__ = "1.1.1"
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

    return task, abunds, os.path.join(out_dir, "merge_mpa.tsv")


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
    Rscript {script}/plot_tree_bar.R abundance_{{sorts}}_top{top}.tsv {group} {{sorts}} tree_abundance
    cp group.abundance_top{top}_{{sorts}}.p* tree_abundance_{{sorts}}.p* {out_dir}
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


def create_heatmap_tasks(abunds, group, job_type, work_dir="", out_dir="", top=50):

    sorts = []
    abundans = []
    for i in abunds:
        abundans.append(abunds[i])
        sorts.append(TXADIC[i])
    tasks = ParallelTask(
        id="heatmap",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
export PATH={python}:$PATH
#python {script}/cut_relative_abundance.py {{abunds}} --top {top} -m zscore >{{sorts}}_top{top}.tsv
python {script}/cut_relative_abundance.py {{abunds}} --top {top} -m log >{{sorts}}_top{top}.tsv
if [ -n "{group}" ]; then
    {rbin}/Rscript {script}/group_heatmap.R {{sorts}}_top{top}.tsv {group} {{sorts}}_top{top}_heatmap
else
    {rbin}/Rscript {script}/cluster_heatmap.R {{sorts}}_top{top}.tsv {{sorts}}_top{top}_heatmap
fi
cp {{sorts}}_top{top}_heatmap.p* {out_dir}
if [ -n "{group}" ]; then
    python {script}/cut_relative_abundance.py {{abunds}} --group {group} --top {top} -m log >group.{{sorts}}_top{top}.tsv
    {rbin}/Rscript {script}/cluster_heatmap.R group.{{sorts}}_top{top}.tsv group.{{sorts}}_top{top}_heatmap
    cp group.{{sorts}}_top{top}_heatmap.p* {out_dir}
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


def create_rarefaction_curve_tasks(merge_mpa, abunds, group, job_type, work_dir="",
                                  out_dir="", kingdom="Bacteria"):

    sorts = []
    abundans = []
    x = []
    for i in abunds:
        if i not in ["g", "s"]:
            continue
        abundans.append(abunds[i])
        sorts.append(TXADIC[i])
        x.append(i)
        
    tasks = ParallelTask(
        id="rarefaction_curve",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
export PATH={python}:$PATH
if [ -n "{group}" ]; then
    python {script}/plot_rarefaction_curve.py {{abunds}} --group {group} --rnumber 1000 --prefix {{sorts}}
else
    python {script}/plot_rarefaction_curve.py {{abunds}} --rnumber 1000 --prefix {{sorts}}
fi
python {script}/get_relative_abundance.py {merge_mpa} --level {{x}} --kingdom {kingdom} \
  --model original >{{sorts}}.abundance_otu.tsv
{rbin}/Rscript {script}/alpha_diversity.R {{sorts}}.abundance_otu.tsv > {{sorts}}.stat_alpha_diversity.tsv
cp {{sorts}}.rarefaction_curve.p* {{sorts}}.stat_alpha_diversity.tsv {out_dir}
""".format(python=PYTHON_BIN,
           script=SCRIPTS,
           rbin=R_BIN,
           group=group,
           merge_mpa=merge_mpa,
           kingdom=kingdom,
           out_dir=out_dir),
        abunds=abundans,
        sorts=sorts,
        x=x
    )

    return tasks


def create_venn_flower_tasks(abunds, group, job_type, work_dir="", out_dir=""):

    sorts = []
    abundans = []
    for i in abunds:
        abundans.append(abunds[i])
        sorts.append(TXADIC[i])
        
    tasks = ParallelTask(
        id="venn_flower",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
export PATH={python}:$PATH
{rbin}/Rscript {script}/venn_flower.R {{abunds}} {{sorts}}
cp {{sorts}}.venn_flower.p* {out_dir}
if [ -n "{group}" ]; then
    python {script}/cut_relative_abundance.py {{abunds}} --group {group} --top 100000000 >group.{{sorts}}.tsv
    {rbin}/Rscript {script}/venn_flower.R group.{{sorts}}.tsv group_{{sorts}}
    cp group_{{sorts}}.venn_flower.p* {out_dir}
fi
""".format(python=PYTHON_BIN,
           script=SCRIPTS,
           rbin=R_BIN,
           group=group,
           out_dir=out_dir),
        abunds=abundans,
        sorts=sorts,
    )

    return tasks


def create_pca_tasks(abunds, group, job_type, work_dir="", out_dir=""):

    sorts = []
    abundans = []
    for i in abunds:
        abundans.append(abunds[i])
        sorts.append(TXADIC[i])

    tasks = ParallelTask(
        id="pca",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
{rbin}/Rscript {script}/plot_pca.R {{abunds}} {group} {{sorts}} >{{sorts}}.stat_pca.read_tsv
cp {{sorts}}.PCA.p* {out_dir}
{rbin}/Rscript {script}/plot_pcoa.R {{abunds}} {group} {{sorts}}
cp {{sorts}}.PCoA.p* {{sorts}}.pcoa.csv {{sorts}}.compare_site.csv {out_dir}
""".format(rbin=R_BIN,
           script=SCRIPTS,
           group=group,
           out_dir=out_dir),
        abunds=abundans,
        sorts=sorts,
    )

    return tasks

def create_anosim_nmds_tasks(abunds, group, job_type, work_dir="", out_dir=""):

    sorts = []
    abundans = []
    for i in abunds:
        abundans.append(abunds[i])
        sorts.append(TXADIC[i])

    tasks = ParallelTask(
        id="anosim_nmds",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
{rbin}/Rscript {script}/plot_anosim.R {{abunds}} {group} {{sorts}} 
cp {{sorts}}*anosim.p* {{sorts}}.ANOSIM.tsv {out_dir}
{rbin}/Rscript {script}/plot_nmds.R {{abunds}} {group} {{sorts}}
cp {{sorts}}.nmds.p* {out_dir}
""".format(rbin=R_BIN,
           script=SCRIPTS,
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
rows=`cat {{sorts}}_lefse.tsv|wc -l`
if (("$rows" >= 2)); then
  python {lefse}/format_input.py {{sorts}}_lefse.tsv {{sorts}}_lefse.in -c 1 -o 1000000
  python {lefse}/run_lefse.py {{sorts}}_lefse.in {{sorts}}_lefse.res -l 2.0
  python {lefse}/plot_res.py {{sorts}}_lefse.res {{sorts}}_lefse.png --dpi 500 --format png --width 20
  python {lefse}/plot_res.py {{sorts}}_lefse.res {{sorts}}_lefse.pdf --format pdf --width 20
  python {script}/filter_lefse.py {{sorts}}_lefse.res --display 1 >{{sorts}}_lefse_filter.tsv
  if [ -s ./{{sorts}}_lefse.png ] ; then
      cp {{sorts}}_lefse.png {{sorts}}_lefse.pdf {{sorts}}_lefse.res {{sorts}}_lefse_filter.tsv {out_dir}
  fi

  if [ "{{sorts}}" = "species" ] ; then
    python {script}/abundance2lefse.py {{abunds}} -g {group} >{{sorts}}_lefse.tsv
    python {lefse}/format_input.py {{sorts}}_lefse.tsv {{sorts}}_lefse.in -c 1 -o 1000000
    python {lefse}/run_lefse.py {{sorts}}_lefse.in {{sorts}}_lefse.res -l 2.0
    python {script}/filter_lefse.py {{sorts}}_lefse.res >{{sorts}}.tsv
    python {lefse}/plot_cladogram.py --dpi 500 --format png {{sorts}}.tsv lefse_cladogram.png
    python {lefse}/plot_cladogram.py --dpi 500 --format pdf {{sorts}}.tsv lefse_cladogram.pdf
    cp lefse_cladogram.png lefse_cladogram.pdf {out_dir}
  fi
fi
rm -rf {{sorts}}_lefse.in
""".format(lefse=LEFSE_BIN,
           python=PYTHON_BIN,
           script=SCRIPTS,
           group=group,
           out_dir=out_dir),
        abunds=abundans,
        sorts=sorts,
    )

    return tasks


def create_test_tasks(abunds, group, job_type, work_dir="", out_dir="", suffix=""):

    sorts = []
    abundans = []
    for i in abunds:
        abundans.append(abunds[i])
        sorts.append(TXADIC[i])

    tasks = ParallelTask(
        id="test%s" % suffix,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
python {script}/abundance2lefse.py {{abunds}} -g {group} --cut >{{sorts}}_lefse.tsv
python {script}/difference_test.py {{sorts}}_lefse.tsv --norma 1000 >{{sorts}}.2-sample_t-test.tsv
cp {{sorts}}.2-sample_t-test.tsv {out_dir}
""".format(python=PYTHON_BIN,
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
        "heatmap": "heatmap",
        "rarefaction_curve": "rarefaction_curve",
        "venn": "venn_flower",
        "pca": "pca",
        "nmds": "anosim_nmds",
        "lefse": "lefse",
        "test": "2-sample_t-test",
    }
    for k, v in work_dict.items():
        if not group and (k in ["pca", "nmds", "lefse", "test"]):
            continue
        mkdir(os.path.join(work_dir, v))
        mkdir(os.path.join(out_dir, v))

    dag = DAG("run_compare")
    mpa2otu_task, abunds, merge_mpa = create_merge_mpa2otu_task(
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
    
    heatmap_tasks = create_heatmap_tasks(
        abunds=abunds,
        group=group,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["heatmap"]),
        out_dir=os.path.join(out_dir, work_dict["heatmap"]),
    )
    dag.add_task(*heatmap_tasks)
    mpa2otu_task.set_downstream(*heatmap_tasks)

    rare_tasks = create_rarefaction_curve_tasks(
        merge_mpa=merge_mpa,
        abunds=abunds,
        group=group,
        job_type=job_type,
        kingdom=kingdom,
        work_dir=os.path.join(work_dir, work_dict["rarefaction_curve"]),
        out_dir=os.path.join(out_dir, work_dict["rarefaction_curve"])
    )
    dag.add_task(*rare_tasks)
    mpa2otu_task.set_downstream(*rare_tasks)

    venn_tasks = create_venn_flower_tasks(
        abunds=abunds,
        group=group,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["venn"]),
        out_dir=os.path.join(out_dir, work_dict["venn"])
    )
    dag.add_task(*venn_tasks)
    mpa2otu_task.set_downstream(*venn_tasks)

    if group:
        pac_tasks = create_pca_tasks(
            abunds=abunds,
            group=group,
            job_type=job_type,
            work_dir=os.path.join(work_dir, work_dict["pca"]),
            out_dir=os.path.join(out_dir, work_dict["pca"])
        )
        dag.add_task(*pac_tasks)
        mpa2otu_task.set_downstream(*pac_tasks)
        nmds_tasks = create_anosim_nmds_tasks(
            abunds=abunds,
            group=group,
            job_type=job_type,
            work_dir=os.path.join(work_dir, work_dict["nmds"]),
            out_dir=os.path.join(out_dir, work_dict["nmds"])
        )
        dag.add_task(*nmds_tasks)
        mpa2otu_task.set_downstream(*nmds_tasks)

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
        
        temp_work = mkdir(os.path.join(os.path.join(work_dir, work_dict["test"]), i))
        temp_out =  mkdir(os.path.join(os.path.join(out_dir, work_dict["test"]), i))
        test_tasks = create_test_tasks(
            abunds=abunds,
            group=groups[i],
            job_type=job_type,
            work_dir=temp_work,
            out_dir=temp_out,
            suffix=i
        )
        dag.add_task(*test_tasks)
        mpa2otu_task.set_downstream(*test_tasks)
        

    do_dag(dag, concurrent, refresh)

    return 0


def mngs_compare(args):

    kingdoms =args.kingdom.strip()

    for i in kingdoms.split(","):
        work_dir = mkdir(os.path.join(args.work_dir, i))
        out_dir = mkdir(os.path.join(args.out_dir, i))
        run_compare(
            mpas=args.input,
            group=args.group,
            vslist=args.vslist,
            kingdom=i,
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
    parser.add_argument("--kingdom", default="Bacteria,Viruses",
       choices=["Bacteria", "Viruses", "Fungi", "Bacteria,Viruses", "Bacteria,Viruses,Fungi"],
       help="Kingdom of the sample (default: Bacteria,Viruses)"),
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
