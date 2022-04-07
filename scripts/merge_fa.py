#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_fasta(file):

    '''Read fasta file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    seq = []
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            line = line.strip(">")
            if len(seq) == 2:
                yield seq
            seq = []
            seq.append(line)
            continue
        if len(seq) == 2:
            seq[1] += line
        else:
            seq.append(line)

    if len(seq) == 2:
        yield seq
    fp.close()


def format_seq(seq, length=100):

    r = re.findall('.{'+str(length)+'}', seq)
    r.append(seq[(len(r)*length):])

    return "\n".join(r)


def get_prefix(file):

    file = file.split('/')[-1]

    if "." in file:
        prefix = file.split(".")
    else:
        prefix = file.split("_")

    return prefix[0]


def merge_fa(files):

    for file in files:
        prefix = get_prefix(file)

        for seqid, seq in read_fasta(file):
            print(">%s|%s\n%s" % (prefix, seqid, format_seq(seq)))

    return 0



def add_hlep_args(parser):

    parser.add_argument("input", nargs="+", metavar="FILE", type=str,
        help="Input genome sequence file, (*.fasta)")

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
    merge_fa.py: Merge genome sequences from each sample
attention:
    merge_fa.py *.fasta >all.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    merge_fa(args.input)


if __name__ == "__main__":

    main()
