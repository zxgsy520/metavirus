#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import sys
import pysam
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
            seq[1] += line.upper()
        else:
            seq.append(line)

    if len(seq) == 2:
        yield seq
    fp.close()


def read_bam_coverage(file):

    if file.endswith(".bam"):
        fh = pysam.AlignmentFile(file, "rb", check_sq=False)
    elif file.endswith(".sam"):
        fh = pysam.AlignmentFile(file, 'r')
    else:
        raise Exception("%r file format error" % file)

    for i in fh.references:
        yield i, fh.count_coverage(contig=i)

    fh.close()


def stat_mutation(genome, bam):

    seqs = {}
    for seqid, seq in read_fasta(genome):
        seqs[seqid] = seq

#   RB = {0:"A", 1:"C", 2:"G", 3:"T"}
    print("""#Reference id\tPosition\tReference base\tReference depth\tDepth\
\tA\tC\tG\tT\tMaxmutation rate(%)""")
    for seqid, covs in read_bam_coverage(bam):
        for i in range(len(seqs[seqid])):
            base = seqs[seqid][i]
            a = covs[0][i]
            c = covs[1][i]
            g = covs[3][i]
            t = covs[3][i]
            tall = a+c+g+t
            if tall == 0:
               tall = 1
            rd = 0
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
            mut = max(dl)*100.0/tall

            print("""{ref}\t{pos}\t{base}\t{refd}\t{dp}\t{a}\t{c}\t{g}\t{t}\
\t{mut:.2f}""".format(
                ref=seqid,
                pos=i+1,
                base=base,
                refd=rd,
                dp=tall,
                a=a,
                c=c,
                g=g,
                t=t,
                mut=mut)
            )

    return 0


def add_hlep_args(parser):

    parser.add_argument('bam', metavar='FILE', type=str,
        help='Input the bam file after genome alignment.')
    parser.add_argument('-g', '--genome', metavar='FILE', type=str, required=True,
        help='Input genome file, format(fasta).')

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
    bam2mutation.py Count the mutation frequency of each point of the virus genome.

attention:
    bam2mutation.py input.bam -g genome.fasta >mutation.xls

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    stat_mutation(args.genome, args.bam)


if __name__ == "__main__":

    main()
