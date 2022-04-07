#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

from scipy import stats
from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line: #line.startswith("#")
            continue

        yield line.split(sep)

    fh.close()


def get_index(header):

    r = OrderedDict()
    n = 1

    for i in header[1::]:
        if i not in r:
            r[i] = []
        r[i].append(n)
        n += 1

    return r


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


def read_abundance(file):

    data = OrderedDict()
    groups = []
    total_abunds = []

    n = 0
    for line in read_tsv(file, "\t"):
        n += 1
        if n ==1:
            groups = get_index(line)
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

    return groups, total_abunds, data


def process_normalization_abundance(abunds, otulist, norma=10):

    r = []
    n = 0

    for i in otulist:
        if abunds[n] == 0:
            r.append("0")
        else:
            r.append(i*100.0/abunds[n]*norma)
        n += 1

    return r


def get_group_abundance(groups, abunds):

    r = OrderedDict()
    group = []

    for i in groups:
        if i not in group:
            group.append(i)
            r[i] = []
        for j in groups[i]:
            r[i].append(abunds[j-1])

    return r[group[0]], r[group[1]], group[0], group[1]


def difference_test(file, norma=1):

    groups, abunds, data =  read_abundance(file)

    print("#Tax_id\tGroup1\tGroup\tT-value\tP-value")
    for taxid, line in data.items():
        if norma:
            line = process_normalization_abundance(abunds, line, norma)

        a, b, aname, bname =  get_group_abundance(groups, line)
        try:
            r = stats.ttest_ind(a, b, alternative="two-sided")
        except:
            ds = stats.levene(a, b)
            if ds.pvalue > 0.05:
                r = stats.ttest_ind(a, b, equal_var=False)
            else:
                r = stats.ttest_ind(a, b, equal_var=False)

        if r.pvalue > 0.05:
            continue
        print("%s\t%s\t%s\t%s\t%s" % (taxid, aname, bname, r.statistic, r.pvalue))
    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input abundance files for 2 groups, lefse.tsv")
    parser.add_argument("-n", "--norma", metavar="INT", type=int, default=1,
        help="set the normalization value (default 0 meaning no normalization), default=1")

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
        difference_test.py lefse.tsv > difference_test.tsv

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    difference_test(args.input, args.norma)


if __name__ == "__main__":

    main()
