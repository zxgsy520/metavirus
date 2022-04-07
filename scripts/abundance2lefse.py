#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.1.1"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line: #line.startswith("#")
            continue

        yield line.split(sep)

    fh.close()


def read_group(file):

    data = {}

    for line in read_tsv(file, "\t"):
        if line[0].startswith("#"):
            continue
        data[line[0]] = line[1]

    return data


def get_required_columns(header, group_data):

    cindex = [0]
    r = ["Groups"]
    n = 1

    for i in header[1::]:
        if i in group_data:
            cindex.append(n)
            r.append(group_data[i])
        n += 1

    return r, cindex


def get_abundance(abundance, cindex):

    r = []
    n = 0 
    tabund = 0

    for i in cindex:
        r.append(abundance[i])
        if n != 0:
            tabund += float(abundance[i])
        n += 1

    if tabund == 0:
        r = []
    return r


def abundance2lefse(abundance, group, cut):

    data = read_group(group)
    header = []
    n = 0

    for line in read_tsv(abundance, "\t"):
        n += 1
        if n == 1:
            header, cindex = get_required_columns(line, data)
            print("\t".join(header))
            continue
        if len(line) != len(cindex):
            line = get_abundance(line, cindex)
            if not line:
                continue
        if cut:
            line[0] = line[0].split("|")[-1]
        print("\t".join(line))

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input the species abundance table for each sample, abundance_genus.xls")
    parser.add_argument("-g", "--group", metavar="FILE", type=str, required=True,
        help="Input sample grouping table,  group.list.")
    parser.add_argument("--cut", action="store_true", default=False,
        help="Whether to trim taxid, default=False.")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
For exmple:
        abundance2lefse.py abundance_order.xls -g group.list >abundance_order_lefse.tsv
group.list:
#sample	group
C1	C_0D
C2	C_0D

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    abundance2lefse(args.input, args.group, args.cut)


if __name__ == "__main__":

    main()
