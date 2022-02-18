#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fa = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fa = open(file)
    else:
        raise Exception("%r file format error" % file)

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


def read_fastq(file):
    '''Read fastq file'''

    if file.endswith(".gz"):
        fp = gzip.open(file, 'r')
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
        if line.startswith("@") and (len(seq)==0 or len(seq)>=5):
            seq = []
            seq.append(line.strip("@").split()[0])
            continue
        if line.startswith("@") and len(seq)==4:
            yield seq[0], seq[1]
            seq = []
            seq.append(line.strip("@").split()[0])
            continue
        seq.append(line)

    if len(seq)==4:
        yield seq[0], seq[1]
    fp.close()


def sort_genome(genome):

    genome_dict = {}

    if genome.endswith(".fastq") or genome.endswith(".fq") or genome.endswith(".fastq.gz") or genome.endswith(".fq.gz"):
        fh = read_fastq(genome)
    elif genome.endswith(".fasta") or genome.endswith(".fa") or genome.endswith(".fasta.gz") or genome.endswith(".fa.gz"):
        fh = read_fasta(genome)
    else:
        raise Exception("%r file format error" % genome)

    for seqid, seq in fh:
        genome_dict[seqid] = seq

    sn = 0
    cn = 0
    for line in sorted(genome_dict.items(),key = lambda x:len(x[1]),reverse = True):
        seq = line[1].upper()

        if seq.count('N')==0:
            cn += 1
            seqid = 'contig{0:06d}'.format(cn)
        else:
            sn += 1
            seqid = 'scaffold{0:06d}'.format(sn)

        print('>%s\n%s' % (seqid, seq))


def add_help(parser):

    parser.add_argument('genome',
        help='Input genome file.')

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
    sort_genome.py  Sort the genome and format it

attention:
    sort_genome.py  genome.fa >genome_sort.fa

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help(parser).parse_args()

    sort_genome(args.genome)


if __name__ == "__main__":

    main()
