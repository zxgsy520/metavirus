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
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            line = line.strip(">")
            if len(seq) == 2:
                yield seq
            seq = []
            seq.append(line.split()[0])
            continue
        if len(seq) == 2:
            seq[1] += line
        else:
            seq.append(line)

    if len(seq) == 2:
        yield seq
    fp.close()


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


def get_seq(files, genome):


    seqids = set()

    for file in files:
        for line in read_tsv(file, "\t"):
            seqids.add(line[0].split("|")[0])

    for seqid, seq in read_fasta(genome):
        if seqid not in seqids:
            continue
        print(">%s\n%s" % (seqid, seq))

    return 0


def add_hlep_args(parser):

    parser.add_argument('genome', metavar='FILE', type=str,
        help='Input genome file.')
    parser.add_argument("-i", "--input", nargs='+', metavar="FILE", type=str, required=True,
        help="Input the id information of the virus sequence.")

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
    get_seq.py: Extract sequence based on sequence id
attention:
    get_seq.py genome.fasta -i dvfpred.txt >virus.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    get_seq(args.input, args.genome)


if __name__ == "__main__":

    main()
