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


def read_tsv(file):

    r = []

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        r.append(line.split("\t"))

    return r


def read_biotool_stat(file):

    line = read_tsv(file)[0]

    return line[1::]


def read_table_annotate(file):

    return read_tsv(file)


def read_table_card_vfdb(files):

    for file in files:
        name = file.split("/")[-1].lower()
        if "card" in name:
            card = read_tsv(file)
        elif "vfdb" in name:
            vfdb = read_tsv(file)
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


def run_report(project, id, sample, sequencer, table_data, table_tax, table_contig,
               table_gene, table_function, table_advanal, table_card_vfdb,
               figure_krona, figure_protein, figure_species, figure_cog,
               figure_kegg, figure_go, figure_cazy, figure_phi,
               tpl_html, out_dir):

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
        "table_data": read_tsv(table_data)[0],
        "table_tax": read_tsv(table_tax)[0],
        "table_contig": read_biotool_stat(table_contig),
        "table_gene": read_biotool_stat(table_gene),
        "table_function": read_table_annotate(table_function),
        "table_advanal": read_table_annotate(table_advanal),
        "table_card_vfdb": read_table_card_vfdb(table_card_vfdb),
        "software": {},
        "database": {}

    }


    figures = [figure_krona, figure_protein, figure_species, figure_cog, figure_kegg, figure_go, figure_cazy, figure_phi]

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
        table_data=args.data,
        table_tax=args.tax,
        table_contig=args.contig,
        table_gene=args.gene,
        table_function=args.function,
        table_advanal=args.advanal,
        table_card_vfdb=args.card_vfdb,
        figure_krona=args.krona,
        figure_protein=args.protein,
        figure_species=args.species,
        figure_cog=args.cog,
        figure_kegg=args.kegg,
        figure_go=args.go,
        figure_cazy=args.cazy,
        figure_phi=args.phi,
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
    table_group.add_argument("--data", metavar="FILE", type=str, required=True,
        help="Statistics before and after quality control.")
    table_group.add_argument("--tax", metavar="FILE", type=str, required=True,
        help="Statistical results of reads species annotation.")
    table_group.add_argument("--contig", metavar="FILE", type=str, required=True,
        help="Assembled Genome Statistics")
    table_group.add_argument("--gene", metavar="FILE", type=str, required=True,
        help="Gene number and length statistics")
    table_group.add_argument("--function", required=True, help="")
    table_group.add_argument("--advanal", required=True, help="")
    table_group.add_argument("--card_vfdb", nargs='+', metavar='FILE', type=str, required=True,
        help="")

    figure_group = parser.add_argument_group(title="Figure", )
    figure_group.add_argument("--krona", required=True, help="")
    figure_group.add_argument("--protein", required=True, help="")
    figure_group.add_argument("--species", required=True, help="")
    figure_group.add_argument("--cog", required=True, help="")
    figure_group.add_argument("--kegg", required=True, help="")
    figure_group.add_argument("--go", required=True, help="")
    figure_group.add_argument("--cazy", required=True, help="")
    figure_group.add_argument("--phi", required=True, help="")

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
