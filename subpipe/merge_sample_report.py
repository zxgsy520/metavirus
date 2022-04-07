#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import json
import argparse
import logging
import os.path
import shutil
from datetime import datetime

from jinja2 import Template
from docxtpl import DocxTemplate

if sys.getdefaultencoding() != "utf8":
    reload(sys)
    sys.setdefaultencoding("utf8")

LOG = logging.getLogger(__name__)

__version__ = "1.2.2"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def check_path(path):

    path = os.path.abspath(path)

    if not os.path.exists(path):
        msg = "File not found '{path}'".format(**locals())
        LOG.error(msg)
        raise Exception(msg)

    return path


def check_paths(obj):
    """
    check the existence of paths
    :param obj:
    :return: abs paths
    """

    if isinstance(obj, list):
        r = []
        for path in obj:
            r.append(check_path(path))

        return r
    else:
        return check_path(obj)


def mkdir(d):
    """
    from FALCON_KIT
    :param d:
    :return:
    """
    d = os.path.abspath(d)
    if not os.path.isdir(d):
        LOG.debug('mkdir {!r}'.format(d))
        os.makedirs(d)
    else:
        LOG.debug('mkdir {!r}, {!r} exist'.format(d, d))

    return d


def read_tsv(file, sep=None):

    r = []

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        r.append(line.split(sep))

    return r


def read_tsvn(file, sep=None):

    for line in open(file):
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line:
            continue

        yield line.split(sep)


def get_prefix(file):

    file = file.split("/")[-1]

    if "." in file:
        prefix = file.split(".")[0]
    else:
        prefix = file.split("_")[0]

    return prefix


def merge_qc(files):

    n = 0
    r = []

    for file in files:
        prefix = get_prefix(file)

        for line in read_tsvn(file):
            n += 1
            if n == 1:
                line[0] = line[0].strip("#")
                print("#Sample\t%s" % "\t".join(line))
                continue
            if not line[0] or line[0].startswith("#"):
                continue
            print("%s\t%s" % (prefix, "\t".join(line)))
            r.append([prefix]+line)

    return r


def mege_read_tax(files):

    r =  merge_qc(files)

    return r


def mege_stat_genome(files):

    n = 0
    r = []

    for file in files:
        prefix = get_prefix(file)

        for line in read_tsvn(file):
            n += 1
            if n == 1:
                line[0] = line[0].strip("#")
                print("\t".join(line))
                continue
            if not line[0] or line[0].startswith("#"):
                continue
            line[0] = prefix
            print("\t".join(line))
            r.append(line)

    return r


def mege_stat_gene(files):

    r =  mege_stat_genome(files)

    return r


def read_alpha_diversity(file):

    n = 0
    r = []
    for line in read_tsv(file):
        n += 1
        if n == 1:
            continue
        r.append(line)

    return r

def read_table_annotate(file):

    return read_tsv(file, "\t")


def read_table_card_vfdb(files):

    for file in files:
        name = file.split("/")[-1].lower()
        if "card" in name:
            card = read_tsv(file, "\t")
        elif "vfdb" in name:
            vfdb = read_tsv(file, "\t")
        else:
            pass
    r = []
    maxlen = max([len(vfdb), len(card)])
    for i in range(maxlen):
        temp = []
        if i >= len(card):
            temp = ["", "", ""]
        else:
            temp = card[i][0:3]
        if i >= len(vfdb):
            temp += ["", "", ""]
        else:
            temp += vfdb[i][0:3]
        r.append(temp)

    return r


def run_report(project, id, sample, sequencer, table_datas, table_taxs, alpha_diversity,
               table_contigs, table_genes, table_function, table_advanal, table_card_vfdb,
               figures, tpl_html, out_dir):

    out_dir = mkdir(out_dir)
    now = datetime.now()

    r = {
        "author": "百易汇能",
        "reviewer": "百易汇能",
        "year": now.year,
        "month": now.month,
        "day": now.day,
        "project": project,
        "id": id,
        "sample": sample,
        "sequencer": sequencer,
        "table_data": merge_qc(table_datas),
        "table_tax": mege_read_tax(table_taxs),
        "alpha_diversity": read_alpha_diversity(alpha_diversity),
        "table_contig": mege_stat_genome(table_contigs),
        "table_gene": mege_stat_gene(table_genes),
        "table_function": read_table_annotate(table_function),
        "table_advanal": read_table_annotate(table_advanal),
        "table_card_vfdb": read_table_card_vfdb(table_card_vfdb),
        "software": {},
        "database": {}

    }

    out_dir = mkdir(out_dir)

    for i in ["images", "static"]:
        temp = os.path.join(out_dir, i)
        if os.path.exists(temp):
            shutil.rmtree(temp)
        shutil.copytree(os.path.join(tpl_html, i), temp)

    for i in check_paths(figures):
        shutil.copy(i, os.path.join(out_dir, "images/"))

    line = open(os.path.join(tpl_html, "report.html")).read()
    if isinstance(line, bytes):
        line = line.decode("utf-8")
    tpl = Template(line)

    with open(os.path.join(out_dir, "report.html"), "w") as fh:
        line = tpl.render(r)

        if isinstance(line, bytes):
            line = line.decode("utf-8")
        fh.write(line)

    return r


def report(args):

    run_report(
        project=args.project,
        id=args.id,
        sample=args.sample,
        sequencer=args.sequencer,
        table_datas=args.data,
        table_taxs=args.tax,
        alpha_diversity=args.alpha_diversity,
        table_contigs=args.contig,
        table_genes=args.gene,
        table_function=args.function,
        table_advanal=args.advanal,
        table_card_vfdb=args.card_vfdb,
        figures=args.figures,
        tpl_html=args.html,
        out_dir=args.out
    )


def add_report_args(parser):

    parser.add_argument("--sample", metavar="STR", type=str, required=True,
        help="Project sample name.")
    parser.add_argument("-p", "--project", metavar="STR", type=str, required=True,
        help="Input project name.")
    parser.add_argument("-i", "--id", metavar="STR", type=str, required=True,
        help="Input project id.")
    parser.add_argument("-s", "--sequencer", metavar="STR", type=str, default="mgi",
        help="Input sequencing platform, default=mgi")
    parser.add_argument("-o", "--out", metavar="FILE", type=str, default="./",
        help="output result path.")
    table_group = parser.add_argument_group(title="Tables", )
    table_group.add_argument("--data", nargs="+", metavar="FILE", type=str, required=True,
        help="Statistics before and after quality control.")
    table_group.add_argument("--tax", nargs="+", metavar="FILE", type=str, required=True,
        help="Statistical results of reads species annotation.")
    table_group.add_argument("-ad", "--alpha_diversity", metavar="FILE", type=str, required=True,
        help="Statistical results of reads species annotation.")
    table_group.add_argument("--contig", nargs="+", metavar="FILE", type=str, required=True,
        help="Assembled Genome Statistics")
    table_group.add_argument("--gene", nargs="+", metavar="FILE", type=str, required=True,
        help="Gene number and length statistics")
    table_group.add_argument("--function", required=True, help="")
    table_group.add_argument("--advanal", required=True, help="")
    table_group.add_argument("--card_vfdb", nargs="+", metavar='FILE', type=str, required=True,
        help="")

    figure_group = parser.add_argument_group(title="Figure", )
    figure_group.add_argument("--figures", nargs="+", required=True, help="")

    template_group = parser.add_argument_group(title="Template", )
    template_group.add_argument("--html", required=True, help="")

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

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_report_args(parser)
    args = parser.parse_args()

    report(args)


if __name__ == "__main__":
    main()
