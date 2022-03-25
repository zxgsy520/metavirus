#!/usr/bin/env python
#coding:utf-8

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
    else:
        fp = open(file)

    r = ""
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            if r:
                yield r.split("\n", 1)
            r = "%s\n" % line.strip(">")
            continue
        r += line.upper()

    if r:
        yield r.split("\n", 1)
    fp.close()


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


def read_vibrant_fasta(file, minlen="10kb"):

    fname = file.split("/")[-1]
    ptypes = fname.split(".")[-2]
    phage, ptypes = ptypes.split("_", 1)

    minlen = convert_size(minlen)

    for seqid, seq in read_fasta(file):
        if len(seq) < minlen:
            continue
        seqid = seqid.split()[0]
        print("%s\t%s\t%s" % (seqid, phage, ptypes))

    return 0


def read_vibrant_phage(files, minlen):

    print("#Seqid\tClassify\tTypes")
    for file in files:
        read_vibrant_fasta(file, minlen)

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", nargs="+", metavar="FILE", type=str,
        help='Input vibrant sequence, format(fna).')
    parser.add_argument("-ml", "--minlen", metavar="STR", type=str, default="5kb",
        help="Set the minimum length of sequence filtering, default=5kb")

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
    read_vibrant_phage.py: Export phages from vibrant prediction results.

attention:
    read_vibrant_phage.py phages_*.fna >stat_vibrant.tsv
    read_vibrant_phage.py phages_*.fna -ml 10kb >stat_vibrant.tsv
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    read_vibrant_phage(args.input, args.minlen)


if __name__ == "__main__":

    main()
