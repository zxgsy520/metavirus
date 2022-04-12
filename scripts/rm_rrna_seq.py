#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.2.0"
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
            if seq:
                yield seq
            seq = [line, ""]
            continue
        seq[1] += line
    if seq:
        yield seq

    fp.close()


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
        if isinstance(line, bytes):
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


def read_tsv(file, step="\t"):

    if file.endswith(".gz"):
        ft = gzip.open(file)
    else:
        ft = open(file)

    for line in ft:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line or line.startswith("#"):
            continue
        yield line.split(step)

    ft.close()


def rm_rrna_seq(file, gff):

    seqs = set()
    data = {}

    for line in read_tsv(gff):
        if line[2] == "rRNA":
            seqs.add(line[0])
            continue
        if line[0] not in data:
            data[line[0]] = set()
        data[line[0]].add(line[2])
        
    LOG.info("The number of deleted sequences is %s" % len(seqs))

    if file.endswith(".fastq.gz") or file.endswith(".fq.gz") or file.endswith(".fastq") or file.endswith(".fq"):
        fp = read_fastq(file)
    elif file.endswith(".fasta.gz") or file.endswith(".fa.gz") or file.endswith(".fasta") or file.endswith(".fa"):
        fp = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)

    for record in fp:
        seqid = record[0].split()[0]
        if seqid in seqs:
            continue
        if seqid in data:
            record[0] += " [structure=%s]" % ";".join(list(data[seqid]))
        print(">%s\n%s" % (record[0], record[1]))
    return 0


def add_help_args(parser):

    parser.add_argument("input",  metavar="FILE", type=str,
        help="Input sequence file(fasta, fastq).")
    parser.add_argument("-g", "--gff", metavar="FILE", type=str, required=True,
        help="gff file of ncRNA annotation results(gff, gff3).")

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
    rm_rrna_seq.py --Delete sequences containing rRNA

attention:
    rm_rrna_seq.py input.fasta --gff input.ncRNA.gff >output.fasta
''')
    args = add_help_args(parser).parse_args()

    rm_rrna_seq(args.input, args.gff)


if __name__ == "__main__":
    main()
