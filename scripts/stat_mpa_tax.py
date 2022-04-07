#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.0.1"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep="\t"):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)

    fh.close()


def split_tax(tax):

    r = OrderedDict()

    for i in tax.split("|"):
        level, value = i.split("__", 1)
        if (level == "k") and (level in r):
            continue
        r[level] = value

    return r


def stat_mpa_tax(file):

    data = {}

    for taxs, reads in  read_tsv(file, sep="\t"):
        taxs = split_tax(taxs)
        index = list(taxs)
        level = index[-1]
        if level not in data:
            data[level] = 0
        data[level] += int(reads)
        if level != "s":
            continue
        if "." not in taxs[level]:
             continue
        level = "sub"
        if level not in data:
            data[level] = 0
        data[level] += int(reads)

    print("#Kingdom\tPhylum\tClass\tOrder\tFamil\tGenus\tSpecies\tSub Species")
    temp = []
    for i in ["k", "p", "c", "o", "f", "g", "s", "sub"]:
        reads = 0
        if i in data:
            reads = data[i]
        temp.append(format(reads, ","))
    print("\t".join(temp))

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input the abundance statistics result file of each sample, kreport2mpa.report")

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
        stat_mpa_tax.py kreport2mpa.report >stat_tax.tsv

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    stat_mpa_tax(args.input)


if __name__ == "__main__":

    main()
