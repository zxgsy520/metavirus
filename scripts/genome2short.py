#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import random
import logging
import argparse


LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def convert_size(string):

    size = string.lower()

    if size.endswith("g") or size.endswith("gb"):
        size, unit = size.split("g", 1)
        size = float(size)*1e9
    elif size.endswith("m") or size.endswith("mb"):
        size, unit = size.split("m", 1)
        size = float(size)*1e6
    elif size.endswith("k") or size.endswith("kb"):
        size, unit = size.split("k", 1)
        size = float(size)*1e3
    else:
        try:
            size = float(size)
        except:
            raise Exception("Size %s input error" % string)

    return size


def read_fasta(file):

    '''Read fastq file'''
    if file.endswith(".gz"):
        fa = gzip.open(file)
    else:
        fa = open(file)
    seq = ''
    for line in fa:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            seq = seq.split('\n')
            if len(seq)==2:
                yield seq[0], seq[1]
            seq = ''
            line = line.strip(">").split()[0]
            seq += "%s\n" % line
            continue
        seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield seq[0], seq[1]
    fa.close()


def cut_seq(seq, minlen=250, maxlen=500):

    seqlen = len(seq)

    if seqlen > maxlen:
        lenseq = random.randint(minlen, maxlen)
        end = random.randint(lenseq, seqlen)
        start = end-lenseq
        seq = seq[start:end]

    return seq


def genome2short(file, window="5kb", minlen=300, maxlen=600):

    window = int(convert_size(window))

    fo = open("contig2sequence.tsv", "w")

    fo.write("#Contig id\tSequence number\n")
    for seqid, seq in read_fasta(file):
        n = 0
        for i in range(window, len(seq), window):
            n += 1
            nseqid = "%s_sub%s" % (seqid, n)
            nseq = seq[i-window:i]
            nseq = cut_seq(nseq)
            print('>%s\n%s' % (nseqid,nseq))
        fo.write("%s\t%s\n" % (seqid, n))
    fo.close()

    return 0


def add_help_args(parser):

    parser.add_argument("genome", metavar="FILE", type=str,
        help='Input genome file(fasta).')
    parser.add_argument("-w", "--window", metavar="STR", type=str, default="5kb",
        help='cut window size, default=5k')
    parser.add_argument("--minlen", metavar="INT", type=int, default=300,
        help="Minimum length of split sequence,  default=300.")
    parser.add_argument("--maxlen", metavar="INT", type=int, default=600,
        help="Maximum length of the split sequence,  default=600.")
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
    genome2short.py:Randomly segmented genomes for blast alignment of annotated species

attention:
    genome2short.py genome.fasta >short.fasta
    genome2short.py genome.fasta --window 10000 >short.fasta
    genome2short.py genome.fasta --window 10kb --minlen 500 --maxlen 1000 >short.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))
    args = add_help_args(parser).parse_args()

    genome2short(args.genome, args.window, args.minlen, args.maxlen)


if __name__ == "__main__":
    main()
