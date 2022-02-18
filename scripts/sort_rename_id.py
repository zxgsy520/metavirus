#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import gzip
import logging
import argparse

import collections

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []



def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''
    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')

        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            if seq:
                yield seq.split('\n')
            seq = "%s\n" % line
        else:
            seq += line
    if seq:
        yield seq.split('\n')
    fp.close()


def sort_raname_id(file, prefix):

    data = {}
    r = {}
    for seqid, seq in read_fasta(file):
        data[seqid] = seq
        r[seqid] = len(seq)

    n = 1
    for seqid, seqlen in sorted(r.items(), key=lambda x:x[1], reverse=True):
        newid = '%s_%s' % (prefix, n)
        print('>%s\n%s' % (newid, data[seqid]))
        n += 1
    return 0


def add_args(parser):

    parser.add_argument('fasta', help='')
    parser.add_argument('-p', '--prefix', metavar='STR', type=str, default='otu',
        help='Input the sequence prefix.')

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
description:
    sort_raname_id : Sort and rename the genome sequence.

usage:
    sort_raname_id.py final.contigs.fa > genome.fasta
version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_args(parser)
    args = parser.parse_args()

    sort_raname_id(args.fasta, args.prefix)


if __name__ == "__main__":
    main()
