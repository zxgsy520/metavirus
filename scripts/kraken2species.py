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


def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)

    if file.endswith(".gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)

    fp.close()


def read_kraken_out(file):

    r = {}

    for line in read_tsv(file, "\t"):
        if "taxid" in line[2]:
            line[2] = line[2].split("(taxid")[-1].strip(")").strip()
        if line[2] not in r:
            r[line[2]] = []
        r[line[2]].append(line[1])

    return r


def kraken2species(file, kraken_tax):

    data = read_kraken_out(file)

    for line in read_tsv(kraken_tax, "\t"):
        if line[0] not in data:
            continue
        for seqid in data[line[0]]:
            print("%s\t%s\t%s" % (seqid, line[0], line[1]))

    return 0


def add_args(parser):

    parser.add_argument('input', metavar='FILE', type=str,
        help='Input the kraken comment result.')
    parser.add_argument('-t', '--taxonomy', metavar='FILE', type=str, default= '/Work/database/taxonomy/202111/kraken.taxonomy.gz',
        help='Input the species classification file(kraken.taxonomy.gz).')

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

    kraken2species.py -- Get the classification information of reads by comparing files

attention:
    kraken2species.py kraken.out >kraken.tax
    kraken2species.py kraken.out --taxonomy /Work/database/taxonomy/202111/kraken.taxonomy.gz >kraken.tax
''')
    args = add_args(parser).parse_args()

    kraken2species(args.input, args.taxonomy)


if __name__ == "__main__":
    main()
