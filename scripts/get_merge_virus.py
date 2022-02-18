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


def read_deepvirfinder(file):

    data = {}

    for line in read_tsv(file, "\t"):
        if "name" == line[0]:
            continue

        data[line[0]] = line[2]
    return data


def read_virsorter(file):

    data = {}

    for line in read_tsv(file, "\t"):
        if "seqname " == line[0]:
            continue
        seqid = line[0].split("||")[0]
        data[seqid] = line[7]

    return data


def read_blst(file):

    data = {}

    for line in read_tsv(file, "\t"):
        if "Viruses" != line[1]:
            continue
        data[line[0]] = line[5]

    return data


def get_merge_virus(genome, deepvirfinder, virsorter, blst):

    dpf_dict = read_deepvirfinder(deepvirfinder)
    vsort_dict = read_virsorter(virsorter)
    blst_dict = read_blst(blst)

    virus = list(dpf_dict.keys()) + list(vsort_dict.keys()) + list(blst_dict.keys())
    virus = sorted(list(set(virus)))

    LOG.info("#Seq Id\tDeepvirfinder\tVirsorter\tBlstn")
    data = {}
    for i in virus:
        temp = [i]
        newid = i
        if i in dpf_dict:
            temp.append("Y")
        else:
            temp.append("N")
        if i in vsort_dict:
            temp.append("Y")
            newid += " [moltype=%s]" % vsort_dict[i]
        else:
            temp.append("N")
        if i in blst_dict:
            temp.append("Y")
            newid += " [organism=%s]" % blst_dict[i]
        else:
            temp.append("N")
        data[i] = newid
        LOG.info("\t".join(temp))

    for seqid, seq in read_fasta(genome):
        if seqid not in data:
            continue
        print(">%s\n%s" % (data[seqid], seq))

    return 0


def add_help_args(parser):

    parser.add_argument("genome", metavar="FILE", type=str,
        help='Input genome file(fasta).')
    parser.add_argument("-dvf", "--deepvirfinder", metavar="FILE", type=str, required=True,
        help='Input deepvirfinder virus prediction results.')
    parser.add_argument("-vsort", "--virsorter", metavar="FILE", type=str, required=True,
        help="Input virsorter virus prediction results.")
    parser.add_argument("-blst", "--blstn", metavar="FILE", type=str, required=True,
        help="Input blstn virus map results.")

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
    get_merge_virus.py: Pick up viral sequences

attention:
    get_merge_virus.py genome.fasta --deepvirfinder dvfpred.txt --virsorter viral-score.tsv --blstn stat_contig_taxonomy.tsv >WL15.raw_virus.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))
    args = add_help_args(parser).parse_args()

    get_merge_virus(args.genome, args.deepvirfinder, args.virsorter, args.blstn)


if __name__ == "__main__":
    main()
