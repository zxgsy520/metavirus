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
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''
    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()

        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            if len(seq)>=1:
                yield seq.split('\n')
            seq = "%s\n" % line
        else:
            seq += line
    yield seq.split('\n')
    fp.close()


def stat_ann(file):

    n = 0

    for line in read_tsv(file, '\t'):
        n += 1

    return n


def stat_advanal_ann(protein, cazy, phi):

    total = 0

    for seqid, seq in read_fasta(protein):
        if len(seq)<=0:
            continue
        total += 1

    cazy = stat_ann(cazy)
    phi = stat_ann(phi)

    print('#Database\tAnnotated number\t% all')
    print('{0}\t{1:,}\t{2:,.2f}\n{3}\t{4:,}\t{5:,.2f}'.format('cazy', cazy, cazy*100.0/total, 'phi', phi, phi*100.0/total))


def add_hlep_args(parser):

    parser.add_argument('protein',
        help='Input protein sequence, format(fasta).')
    parser.add_argument('-c', '--cazy', metavar='FILE', type=str, required=True,
        help='Input cazy data to annotate results(cazy.tsv).')
    parser.add_argument('-p', '--phi', metavar='FILE', type=str, required=True,
        help='Input cazy data to annotate results(phi.tsv).')

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
    stat_advanal_ann.py Statistics and advanal analysis output results.

attention:
    stat_advanal_ann.py protein.fasta --cazy cazy.tsv --phi phi.tsv
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    stat_advanal_ann(args.protein, args.cazy, args.phi)


if __name__ == "__main__":

    main()
