#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_fastq(file):
    '''Read fastq file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []
    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq.append(line.strip("@"))
            continue
        seq.append(line)

        if len(seq)==4:
            yield seq
            seq = []
    fp.close()


def stat_base(seq):

    seq = seq.upper()
    a = seq.count('A')
    t = seq.count('T')
    g = seq.count('G')
    c = seq.count('C')
    n = seq.count('N')

    return a, t, g, c, n


def stat_reads(files, name):

    nb = 0
    ab = 0
    tb = 0
    cb = 0
    gb = 0
    reads = 0

    for file in files:
        for line in read_fastq(file):
            reads += 1
            a, t, g, c, n = stat_base(line[1])
            nb += n
            ab += a
            tb += t
            cb += c
            gb += g

    base = nb+ab+tb+cb+gb

    print("#Sample\tTotal base\tReads number\tGC\tA\tT\tG\tC\tN")
    print("{0}\t{1:,}\t{2:,}\t{3:.2f}\t{4:.2f}\t{5:.2f}\t{6:.2f}\t{7:.2f}\t{8:.2f}".format(
        name, base, reads, (cb+gb)*100.0/base, ab*100.0/base, tb*100.0/base,
        gb*100.0/base, cb*100.0/base, nb*100.0/base))


def add_help_args(parser):

    parser.add_argument('input', nargs='+', metavar='FILE', type=str,
        help='Input data file, (fastq, fastq.gz)')
    parser.add_argument('-n', '--name', metavar='STR', type=str, default='out',
        help='Input sample name, default=out.')

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
    stat_ngs --Statistical sequencing data

attention:
    stat_ngs *.fastq -n name > name.stat_ngs.tsv
''')
    args = add_help_args(parser).parse_args()
    stat_reads(args.input, args.name)


if __name__ == "__main__":
    main()
