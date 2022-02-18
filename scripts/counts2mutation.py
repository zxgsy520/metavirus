#!/usr/bin/env python
#coding:utf-8

import os
import sys
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):
    """
    read tsv joined with sep
    :param file: file name
    :param sep: separator
    :return: list
    """
    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_stdin(sep=None):

    for line in sys.stdin:
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def filter_count(rd, base, a, c, g, t, mincount=3):

    if a < mincount:
        a = 0
    if c < mincount:
        c = 0
    if g < mincount:
        g = 0
    if t < mincount:
        t = 0
    tall = a+c+g+t
    if tall == 0:
        tall = 1
    dl = [a, c, g, t]
    if base == "A":
        rd = a
        dl = [c, g, t]
    elif base == "C":
        rd = c
        dl = [a, g, t]
    elif base == "G":
        rd = g
        dl = [a, c, t]
    elif base == "T":
        rd = t
        dl = [a, c, g]
    else:
        rd = 0
        dl = [a, c, g, t]
    mut = max(dl)*1000.0/tall

    return rd, a, c, g, t, tall, mut


def counts2mutation(file, mincount=3):

    if file:
        fh = read_tsv(file, '\t')
    else:
        fh = read_stdin('\t')
    snp = 0

    print("""#Reference id\tPosition\tReference base\tReference depth\tDepth\
\tA\tC\tG\tT\tMaxmutation rate(‰)""")
    for line in fh:
        if line[0] == "chr":
            continue
        base = line[3]
        rd, a, c, g, t, tall, mut = filter_count(int(line[2]), base, int(line[6]),
                                                 int(line[7]), int(line[8]),
                                                 int(line[9]), mincount)
        if mut > 0:
            snp += 1

        print("""{ref}\t{pos}\t{base}\t{refd}\t{dp}\t{a}\t{c}\t{g}\t{t}\
\t{mut:.2f}""".format(
            ref=line[0],
            pos=line[1],
            base=base,
            refd=rd,
            dp=tall,
            a=a,
            c=c,
            g=g,
            t=t,
            mut=mut)
        )
    LOG.info("The number of mutation sites is: %s" % snp)

    return 0


def add_hlep_args(parser):

    parser.add_argument('input', metavar='FILE', type=str,
        help='Input the bam file after genome alignment.')
    parser.add_argument('-mc', '--mincount', metavar='INT', type=int, default=3,
        help='Filter base support number.')

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
    counts2mutation.py Count the mutation frequency of each point of the virus genome.

attention:
    counts2mutation.py mpileup_counts.tsv >mutation.xls
    samtools mpileup -aa -f genome.fa rmdup.bam  |mpileup2readcounts |counts2mutation.py >mutation.xls
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    counts2mutation(args.input, args.mincount)


if __name__ == "__main__":

    main()
