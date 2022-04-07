#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep="\t"):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def stat_species(file, name):

    spec_dict = {}

    for line in read_tsv(file):
        spec = line[-1].split(":")[-1]

        if spec not in spec_dict:
            spec_dict[spec] = 1
            continue
        spec_dict[spec] += 1

    n = 0
    output = open("%s.species.tsv" % name, "w")

    output.write("#Species name\tProtein number\n")
    for i in sorted(spec_dict.items(),key = lambda x:x[1],reverse = True):
        if n >=10:
            break
        output.write("%s\t%s\n" % (i[0], i[1]))
        n += 1
           
    output.close()


def add_args(parser):

    parser.add_argument('-i', '--input', metavar='FILE', type=str, required=True,
        help='Input the species annotation file.')
    parser.add_argument('-n', '--name', metavar='STR', type=str, default= 'out',
        help='Output file prefix')
    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
name:

    sure_species.py -- Confirm species by annotation results in the database

attention:
    sure_species.py -i refseq.tsv -n name
''')
    args = add_args(parser).parse_args()

    stat_species(args.input, args.name)


if __name__ == "__main__":
    main()
