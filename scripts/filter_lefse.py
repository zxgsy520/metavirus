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

__version__ = "1.0.0"
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


def filter_lefse(file, display=100):

    n = 0
    r = []

    for line in read_tsv(file, sep="\t"):
        if ("-" != line[-1]) and ("-" != line[2]):
            print("\t".join(line))
            n += 1
        else:
            r.append(line)

    for line in r:
        print("\t".join(line))
        n += 1
        if n >= display:
            break

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input the lefse analysis table, lefse.res")
    parser.add_argument("-d", "--display", metavar="INT", type=int, default=800,
        help="Set the number of species displayed. default=80.")

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
     filter_lefse.py.py lefse.res > lefse_new.res
version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    filter_lefse(args.input, args.display)


if __name__ == "__main__":

    main()
