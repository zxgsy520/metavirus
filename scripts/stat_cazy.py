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


def read_tsv(file, sep=None):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_cazy(file):

    data = {}

    for line in read_tsv(file, "\t"):
        for i in line[2].split(";"):
            if i not in data:
                data[i] = [0, []]
            data[i][0] += 1
            data[i][1].append(line[0])
    return data


def read_activ(file):

    data = {}

    for line in read_tsv(file, "\t"):
        if len(line)<=1:
            line.append("")
        data[line[0]] = line[1].strip()

    return data


def stat_cazy(cazy, activ):

    cazy_dict = read_cazy(cazy)
    activ_dict = read_activ(activ)

    print("#Class\tNumber\tDesc\tGene")
    for gene, line in sorted(cazy_dict.items(), key=lambda x:x[1][0], reverse = True):
        if gene in activ_dict:
            desc = activ_dict[gene]
        else:
            desc = ""
        print("%s\t%s\t%s\t%s" % (gene, line[0], desc, ";".join(line[1])))

    return 0


def add_help(parser):

    parser.add_argument('input', metavar='FILE', type=str,
        help='Input the annotation file of CAZy, CAZy.tsv')
    parser.add_argument('--activ', metavar='FILE', type=str,  required=True,
        help='Input CAZy types of documentation, CAZy.activities.txt.')
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
    stat_cazy.py CAZy Statistical Classification

attention:
    stat_cazy.py CAZy.cazy.tsv --activ CAZy.activities.txt >stat_cazy.tsv

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help(parser).parse_args()

    stat_cazy(args.input, args.activ)


if __name__ == "__main__":

    main()
