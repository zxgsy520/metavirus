#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

import math
from matplotlib import pyplot as plt

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


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


def read_mutation(file):

    r = {}

    for line in read_tsv(file, "\t"):
        r[int(line[1])] = float(line[-1])

    return r


def get_sample(file):

    sample = file.split("/")[-1]
    sample = sample.split(".")[0]

    return sample.split("-")[-1]


def stat_mutation(file1, file2):

    r1 = read_mutation(file1)
    r2 = read_mutation(file2)

    m1 = 0
    m2 = 0
    mall = 0
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    
    for i in r1:
        if r1[i] > 0 and r2[i] >0:
            mall += 1
            x1.append(i)
            x2.append(i)
            y1.append(math.log(r1[i]))
            y2.append(math.log(r2[i]))
            continue    
        if r1[i] > 0:
            m1 += 1
            x1.append(i)
            y1.append(math.log(r1[i]))
        if r2[i] >0:
            m2 += 1
            x2.append(i)
            y2.append(math.log(r2[i]))
    label1 = get_sample(file1)
    label2 = get_sample(file2)
    print("#Unique and public mutation statistics.")
    print("#%s\t%s\tAll" % (label1, label2))
    print("{0:,}\t{1:,}\t{2:,}".format(m1, m2, mall))

    return x1, y1, x2, y2, label1, label2, len(r1)


def add_legend(ax, x, y, location="upper left"):

    ax.legend(loc=location,
        frameon=False,
        fontsize="large",
        labelspacing=0.3,
        handlelength=0.6,
        handleheight=0.6,
        handletextpad=0.5,
        bbox_to_anchor=(x, y)
    )

    return 0


def plot_mutation(x1, y1, x2, y2, label1, label2, xmax):

    reload(sys)
    sys.setdefaultencoding('utf-8')
    fig = plt.figure(figsize=(12, 4.5))
    ax = fig.add_subplot(111)

    plt.grid(True, which='minor', axis='both', lw=1, color='#E5C700', alpha=0.3)
    plt.grid(True, which='major', axis='both', lw=1.5, color='#E2BDD5', alpha=0.3)
    ax.tick_params(axis='both', which='both', color='#ffffff', length=5, width=1.5)
    font = {'family': 'Times New Roman', 'weight': 'normal', 'color': '#212121', 'size': 16}
    ax.set_xlabel('Genomic Position(bp)', font)
    ax.set_ylabel(u'Log10 Mutation Rate(‰)', font)
    ax.plot(x1, y1, "o", color='#212121', markersize=4, label=label1)
    ax.plot(x2, y2, "o", color='#FFA000', markersize=4, label=label2)

    ymax = max(y1+y2)
    plt.ylim(0, ymax)
    plt.xlim(0, xmax)

    add_legend(ax, 1, 1, 'upper left')
    fig.tight_layout()
    plt.savefig('mutation.png', dpi=700)
    plt.savefig("mutation.pdf")

    return 0


def add_hlep_args(parser):

    parser.add_argument('-f1', '--file1', metavar='FILE', type=str,
        help='Input mutation rate file1.')
    parser.add_argument('-f2', '--file2', metavar='FILE', type=str,
        help='Input mutation rate file2.')

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
URL: https://github.com/zxgsy520/nextmeta
name:
    plot_mutation.py: Plot a mutation rate comparison chart.

attention:
    plot_mutation.py -f1 prefix1.mutation.xls -f2 prefix2.mutation.xls >mutation.tsv

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    x1, y1, x2, y2, label1, label2, xmax = stat_mutation(args.file1, args.file2)

    plot_mutation(x1, y1, x2, y2, label1, label2, xmax)


if __name__ == "__main__":

    main()
