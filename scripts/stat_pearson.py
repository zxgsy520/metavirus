#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import logging
import argparse

import numpy as np
from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []



def read_tsv(file, sep=None):

    for line in open(file):
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue

        yield line.split(sep)


def str_list2float(otulist):

    r = []

    for i in otulist:
        r.append(float(i))

    return r


def read_bundance(file):

    n = 0
    samples = []
    data = OrderedDict()

    for line in read_tsv(file, "\t"):
        n += 1
        if n == 1:
            samples = line[1::]
            continue
        for i in range(len(samples)):
            if samples[i] not in data:
                data[samples[i]] = []
            data[samples[i]].append(float(line[i+1]))

        if min(str_list2float(line[1::])) <= 0:
            break

    return data, samples


def process_zscore(otulist):

    otulist = np.array(otulist)
    meanotu = np.mean(otulist)
    std = np.std(otulist)

    return (otulist-meanotu)/std


def stat_pearson(file):

    data, samples = read_bundance(file)

    print("#Sample\t%s" % "\t".join(samples))
    maxbd = max(data[samples[0]])

    for i in data:
        temp = []
        pil = data[i]
        if maxbd >= 100:
            pil = process_zscore(data[i])
        for j in data:
            if i == j:
                temp.append("1")
                continue
            pjl = data[j]
            if maxbd >= 100:
                pjl = process_zscore(data[j])
            result = np.corrcoef(pil, pjl)
            temp.append("{:.6f}".format(result[0][1]))

        print("%s\t%s" % (i, "\t".join(temp)))

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input abundance file, abundance_species.xls")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
stat_pearson.py: Calculate the Pearson correlation coefficient.

For exmple:
        stat_pearson.py abundance_species.xls > stat_pearson.xls

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    stat_pearson(args.input)


if __name__ == "__main__":

    main()
