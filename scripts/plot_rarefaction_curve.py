#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

import math
import random
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

LOG = logging.getLogger(__name__)

__version__ = "v1.1.0"
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


def get_repeat_species(string, repeats):

    if repeats <=1:
        repeats = 1

    r = [string] * int(repeats)

    return r


def read_bundance(file):

    n = 0
    samples = []
    data = {}

    for line in read_tsv(file, "\t"):
        n += 1
        if n == 1:
            samples = line[1::]
            continue
        for i in range(len(samples)):
            if samples[i] not in data:
                data[samples[i]] = []

            abund = float(line[i+1])
            if abund <=0 :
                continue
            data[samples[i]].append(abund)

    return data


def get_diversity(file, rnumber=10e3):

    data = read_bundance(file)
    minab = 1000

    for i in data:
        mintp = min(data[i])
        if mintp >= minab:
            continue
        minab = mintp

    rnum = int(1.0/minab+1)
    if rnum >= rnumber:
        LOG.info("The parameter rnumber is %s, which is seriously too small. The minimum setting is %s." % (rnumber, rnum))
        rnumber = rnum

    for i in data:
        temp = []
        n = 0
        for j in data[i]:
            n += 1
            temp += get_repeat_species("A%s" % n, j*rnumber)
        data[i] = temp

    return data


def stat_saturation(datalist, step=5):

    datalen = len(datalist)
    step = int(5*datalen/100)
    x = []
    y = []
    yerr = []

    for i in range(step, datalen+step, step):
        if i >= datalen:
            i = datalen
        x.append(i*100.0/datalen)
        temp = []
        for j in range(0, 3, 1):
            temp.append(len(set(random.sample(datalist, i))))
        y.append(np.mean(temp))
        yerr.append(np.std(temp, ddof=1))

    return x, y, yerr


def read_group(file):

    data = {}

    for line in read_tsv(file, "\t"):
        if line[0].startswith("#") or line[0]=="sample":
            continue
        data[line[0]] = line[1]

    return data


def plot_rarefaction_curve(file, group, prefix, rnumber=10e3, step=5):

    data = get_diversity(file, rnumber)

    if group:
        group = read_group(group)
  
    font = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 14,
         }
    colors = ["#0e72cc", "#6ca30f", "#f59311", "#fa4343", "#16afcc", "#85c021",
              "#d12a6a", "#0e72cc", "#6ca30f", "#f59311", "#fa4343", "#16afcc"]

    plt.switch_backend("agg")
    #plt.style.use('ggplot')
    fig = plt.figure(figsize=[12, 6.5])
    ax = plt.axes([0.08, 0.08, 0.75, 0.88])
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    #ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_visible(False) #去除边框
    ax.spines['right'].set_visible(False)

    r = {}
    n = 0
    for i in data:
        x, y, yerr = stat_saturation(data[i], step)
        if group:
            gid = group[i]
            if gid not in r:
                r[gid] = colors[n]
                n += 1    
                ax.errorbar(x, y, yerr=yerr, label=gid, color=r[gid],  ecolor=r[gid])
            else:
                ax.errorbar(x, y, yerr=yerr, color=r[gid], ecolor=r[gid])
        else:
            ax.errorbar(x, y, yerr=yerr, label=i)

    ax.set_xlabel("Rarefraction Percentage(%)", font)
    ax.set_ylabel("Richness (Observed OTU)", font)

    ax.legend(loc="center left",
        frameon=False,
        labelspacing=0.5,
        handletextpad=0.3,
        columnspacing=3.0,
        bbox_to_anchor=(1.02, 0.5))
    plt.xlim([0, 100])
    plt.savefig('%s.rarefaction_curve.pdf' % prefix)
    plt.savefig('%s.rarefaction_curve.png' % prefix, dpi=700)

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input abundance file, abundance_species.xls")
    parser.add_argument("-g", "--group", metavar="FILE", type=str, default=None,
        help="Input sample grouping table,  group.list.")
    parser.add_argument("-p", "--prefix", metavar="STR", type=str, default="OUT",
        help="Output file prefix, default=OUT")
    parser.add_argument("-r", "--rnumber", metavar="INT", type=int, default=10e3,
        help="Select the number of reads to simulate, default=10e3")
    parser.add_argument("-s", "--step", metavar="INT", type=int, default=5,
        help="Set the simulation step size (interval percentage), default=5")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
plot_rarefaction_curve.py: Calculate the Pearson correlation coefficient.

For exmple:
        plot_rarefaction_curve.py abundance_species.xls > stat_pearson.xls

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    plot_rarefaction_curve(args.input, args.group, args.prefix, args.rnumber, args.step)


if __name__ == "__main__":

    main()
