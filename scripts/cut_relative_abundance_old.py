#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.0.1"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep="\t"):

    for line in open(file):
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue

        yield line.split(sep)


def addlist(list1, list2):

    r = []
    n = 0
    for i in list2:
        r.append(list1[n]+i)
        n += 1

    return r


def str_list2float(otulist):

    r = []

    for i in otulist:
        r.append(float(i))

    return r


def read_abundance(file, top=20):

    data = OrderedDict()
    samples = []
    total_abunds = []

    n = 0
    for line in read_tsv(file, sep="\t"):
        n += 1
        if n ==1:
            samples = line[1::]
            continue
        taxid = line[0].split("|")[-1]
        abunds = str_list2float(line[1::])
        if len(total_abunds) == 0:
            total_abunds = abunds
        else:
            total_abunds = addlist(total_abunds, abunds)
        if taxid not in data:
            data[taxid] = abunds
        else:
            data[taxid] = addlist(data[taxid], abunds)
        if n >= top+1:
            break

    return samples, total_abunds, data


def process_uniform_abundance(abunds, otulist):

    r = []
    n = 0

    for i in otulist:
        r.append("{:.6f}".format(i*100.0/abunds[n]))
        n += 1

    return r


def cut_relative_abundance(file, top=20):

    samples, abunds, data =  read_abundance(file, top)

    print("Taxid\t%s" % "\t".join(samples))
    for taxid, line in data.items():
        line = process_uniform_abundance(abunds, line)
        print("%s\t%s" % (taxid, "\t".join(line)))

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input merge-sorted abundance file, abundance_species.xls")
    parser.add_argument("-t", "--top", metavar='INT', type=int, default=20,
        help="Species showing top X in abundance, default=20.")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
For exmple:
        cut_relative_abundance.py abundance_species.xls -t 20 > abundance_species_top20.xls

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    cut_relative_abundance(args.input, args.top)


if __name__ == "__main__":

    main()
