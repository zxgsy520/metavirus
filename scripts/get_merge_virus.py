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


def read_virsorter(file, minlen="10kb"):

    data = {}
    minlen = convert_size(minlen)
    for line in read_tsv(file, "\t"):
        if "seqname" in line[0]:
            continue
        if int(line[8]) < minlen:
            continue

        seqid = line[0].split("||")[0]
        data[seqid] = line[7]

    return data


def read_deepvirfinder(file, minlen="10kb", score=0.9, pvalue=0.01):

    r = []
    minlen = convert_size(minlen)
    for line in read_tsv(file, "\t"):
        if "name" == line[0]:
            continue
        if int(line[1]) < minlen:
            continue
        if float(line[2]) < score or float(line[3]) > pvalue:
            continue

        r.append(line[0])
    return r


def read_vibrant(file):

    data = {}

    for line in read_tsv(file, "\t"):
        data[line[0]] = [line[1], line[2]]

    return data


def read_virfinder(file, minlen="10kb", score=0.9, pvalue=0.01):

    r = []
    minlen = convert_size(minlen)
    for line in read_tsv(file, "\t"):
        if int(line[1]) <minlen:
            continue
        if float(line[2]) < score or float(line[3]) > pvalue:
            continue
        r.append(line[0])

    return r


def read_virkraken(file):

    r = []

    for line in read_tsv(file, "\t"):
        if "ContigID" in line[0]:
            continue
        r.append(line[0])

    return r


def get_merge_virus(genome, virsorter, deepvirfinder, vibrant,virfinder, virkraken,
                    minlen="10kb", score=0.9, pvalue=0.01):

    virs_dict = read_virsorter(virsorter, minlen)
    depf_list = read_deepvirfinder(deepvirfinder, minlen, score, pvalue)
    vibr_dict = read_vibrant(vibrant)
    virf_list = read_virfinder(virfinder, minlen, score, pvalue)

    virus = list(virs_dict.keys()) + depf_list + list(vibr_dict.keys()) + virf_list
    if virkraken:
        virk_list = read_virkraken(virkraken)
        virus = virus + virk_list
    else:
        virk_list = []

    virus = sorted(list(set(virus)))
    data = {}
    for i in virus:
        if i not in data:
            data[i] = ["virus", "", "", ""]
        if i in virs_dict:
            data[i][1] = virs_dict[i]
            data[i][3] = "virsorter"
        if i in depf_list:
            data[i][3] += ";deepvirfinder"
        if i in vibr_dict:
            data[i][0] = vibr_dict[i][0]
            data[i][2] = vibr_dict[i][1]
            data[i][3] += ";vibrant"
        if i in virf_list:
            data[i][3] += ";virfinder"
        if i in virk_list:
            data[i][3] += ";virkraken"

    LOG.info("#Seq Id\tOrganism\tSeq type\tPhage type\tSupported software")
    minlen = convert_size(minlen)
    for seqid, seq in read_fasta(genome):
        if len(seq) < minlen:
            continue
        if seqid not in data:
            continue
        data[seqid][3] = data[seqid][3].strip(";")
        LOG.info("%s\t%s" % (seqid, "\t".join(data[seqid])))

        organism, vtype, ptype, software = data[seqid]
        if organism:
            seqid = "%s [organism=%s]" % (seqid, organism)
        if vtype:
            seqid = "%s [moltype=%s]" % (seqid, vtype)
        if ptype:
            seqid = "%s [describe=%s]" % (seqid, ptype)
        if software:
            seqid = "%s [software=%s]" % (seqid, software)

        print(">%s\n%s" % (seqid, seq))

    return 0


def add_help_args(parser):

    parser.add_argument("genome", metavar="FILE", type=str,
        help='Input genome file(fasta).')
    parser.add_argument("-vsort", "--virsorter", metavar="FILE", type=str, required=True,
        help="Input virsorter virus prediction results.")
    parser.add_argument("-devf", "--deepvirfinder", metavar="FILE", type=str, required=True,
        help='Input deepvirfinder virus prediction results.')
    parser.add_argument("-vibr", "--vibrant", metavar="FILE", type=str, required=True,
        help="Input vibrant phages prediction results.")
    parser.add_argument("-virf", "--virfinder", metavar="FILE", type=str, required=True,
        help="Input virfinder virus prediction results.")
    parser.add_argument("-virk", "--virkraken", metavar="FILE", type=str, default="",
        help="Input virkraken virus prediction results.")
    parser.add_argument("-ml", "--minlen", metavar="STR", type=str, default="5kb",
        help="Set the minimum length of sequence filtering, default=5kb")
    parser.add_argument("--score", metavar="FLOAT", type=float, default=0.9,
        help="set score, default=0.9")
    parser.add_argument("--pvalue", metavar="FLOAT", type=float, default=0.01,
        help="Set pvalue, default=0.01")

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
    get_merge_virus.py genome.fasta --virsorter viral-score.tsv --deepvirfinder dvfpred.txt
      --vibrant vibrant.tsv --virfinder virfinder.txt > virus.fasta
   get_merge_virus.py genome.fasta --virsorter viral-score.tsv --deepvirfinder dvfpred.txt
    --vibrant vibrant.tsv --virfinder virfinder.txt --virkraken virus.csv --minlen 10kb --score 0.9 --pvalue 0.01 > virus.fasta

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))
    args = add_help_args(parser).parse_args()

    get_merge_virus(args.genome, args.virsorter, args.deepvirfinder, args.vibrant,
        args.virfinder, args.virkraken, args.minlen, args.score, args.pvalue)


if __name__ == "__main__":
    main()
