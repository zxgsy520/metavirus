#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def split_fa_attr(attributes):

    r = {}
    contents = attributes.replace("[", "").replace("]", "").split()

    for content in contents:
        if not content:
            continue
        if "=" not in content:
            LOG.info("%r is not a good formated attribute: no tag!")
            continue
        tag, value = content.split("=", 1)
        r[tag] = value

    return r


def read_fa_attr(file):

    r = {}

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            seqid, line = line.strip(">").split(" ", 1)
        else:
            continue
        line = split_fa_attr(line)

        if "software" not in line:
            continue
        line = line["software"].split(":")
        r[seqid] = line

    return r


def des_fa2otu(file):

    r = read_fa_attr(file)

    soft = ["virsorter", "deepvirfinder", "vibrant", "virfinder", "virkraken"]
    print("Seq Id\t%s" % "\t".join(soft))
    for seqid in r:
        temp = [seqid]
        for i in soft:
            if i in r[seqid]:
                temp.append("1")
            else:
                temp.append("0")
        print("\t".join(temp))

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input sequence(fasta).")

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
    des_fa2otu.py Convert sequences to software otu

attention:
    des_fa2otu.py unique_virus.fasta >software_otu.tsv

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    des_fa2otu(args.input)


if __name__ == "__main__":

    main()
