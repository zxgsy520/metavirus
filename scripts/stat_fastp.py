#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import argparse
import logging

LOG = logging.getLogger(__name__)

__version__ = "2.0.0"
__author__ = ("Junpeng Fan, Xingguo Zhang",)
__email__ = "jpfan@whu.edu.cn, 1131978210@qq.com"
__all__ = []


def json2tsv(file):

    """
    sample  total_reads  total_bases    clean_reads    clean_bases   Q20 rate (%)    Q30 rate (%)   GC

    :param json:
    :return:
    """
    j = json.load(open(file))

    return [j["summary"]["before_filtering"]["total_reads"],
            j["summary"]["before_filtering"]["total_bases"],
            j["summary"]["after_filtering"]["total_reads"],
            j["summary"]["after_filtering"]["total_bases"],
            j["summary"]["after_filtering"]["q20_rate"]*100,
            j["summary"]["after_filtering"]["q30_rate"]*100,
            j["summary"]["after_filtering"]["gc_content"]*100
            ]


def sum_data(files):

    r = []
    for file in files:
        line = json2tsv(file)
        if len(r)==0:
            r = line
            continue
        temp = []
        for i,j in zip(r, line):
            n = i+j
            temp.append(n)
        r = temp

    return r


def stat_fastp(files):

    r = sum_data(files)
    sn = len(files)
    r[-1] = r[-1]/sn
    r[-2] = r[-2]/sn
    r[-3] = r[-3]/sn
    print("""\
#Total reads\ttotal bases\tclean reads\tclean bases\tQ20 rate (%)\tQ30 rate (%)\tGC (%)
{:,}\t{:,}\t{:,}\t{:,}\t{:.2f}\t{:.2f}\t{:.2f}\
""".format(*r))


def add_args(parser):
    parser.add_argument("input", nargs='+',
       help="Input the json file obtained by fastp.")
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
stat_fastp :Consolidate statistics of fastp quality control results
stat_fastp -i *.fastp.json >stat.fastp.tsv
version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args = add_args(parser).parse_args()

    stat_fastp(args.input)


if __name__ == "__main__":
    main()
