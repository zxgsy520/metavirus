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


def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)
    if file.endswith("gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)
    fp.close()


def read_abricate(file):

    print("#qseqid\tsseqid\tdbxref\tname\tdesc\tqstart\tqend\tstartend\tcoverage\tidentity")
    r = {}
    for line in read_tsv(file, "\t"):
        print("{0}\t{1}\t{2}:{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}".format(
            line[1], line[12], line[11], line[12], line[5], line[13],
            line[2], line[3], line[6], line[9], line[10])
        )
        if line[5] not in r:
            r[line[5]] = [[], line[13]]
        r[line[5]][0].append(line[1])

    return r


def abricate_summary(file, out):

    data = read_abricate(file)
    fh = open(out, "w")

    fh.write("#Gene_name\tSeqid\tGenes\tFunction\n")
    for i in data:
        fh.write("{0}\t{1}\t{2}\t{3}\n".format(
            i, ",".join(data[i][0]), len(data[i][0]), data[i][1])
        )
    fh.close()

    return 0


def add_hlep_args(parser):

    parser.add_argument("abricate", metavar='FILE', type=str,
        help="Input the result of abricate comparison.")
    parser.add_argument("-o", "--out", metavar='STR', type=str, default="out",
        help="Output result file.")

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
    abricate_summary.py Summarize the abricate analysis results.
attention:
    abricate_summary.py vfdb.tsv -o vfdb.summary.tsv >vfdb.stat.tsv
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    abricate_summary(args.abricate, args.out)


if __name__ == "__main__":

    main()
