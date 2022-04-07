#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
from collections import OrderedDict

__author__, __email__, __version__ = ("fanjunpeng",), "jpfan@whu.edu.cn", "0.1.0"

LOG = logging.getLogger(__name__)

__all__ = []


def read_tbl(file):
    """
    read table
    :param file:
    :return:
    """

    for line in open(file):
        line = line.strip()

        if line.startswith("#") or not line:
            continue

        yield line.split("\t")


def cluster_protein(file):

    """
    cluster protein by pathway and ko
    :param file: kegg annotation result consist protein id, ko and pathways joined with "\t"
    :return: dict {pathway: {ko: [proteins]}}
    """
    path_dict = {}
    LOG.info("reading kegg result from '%r'" % file)

    for protein, ko, pathway in read_tbl(file):
        paths = pathway.split(";")

        if ko.startswith("ko:"):
            ko = ko[3:]

        for path in paths:

            if path.startswith("path:"):
                path = path[5:]

            if path not in path_dict:
                path_dict[path] = {}

            if ko not in path_dict[path]:
                path_dict[path][ko] = []

            path_dict[path][ko].append(protein)

    return path_dict


def output_keg(keg, path_dict, output):
    """
    output .keg by kegg annotation result
    :param keg: ko00001.keg
    :param path_dict: see function cluster_protein
    :param output: output file
    :return: 0
    """
    path_id = ""

    LOG.info("output kegg map to '%r'" % output)
    fh = open(output, "w")

    for line in open(keg):
        line = line.strip()

        if not line:
            continue

        tag = line[0]

        if tag == "C":
            path_id = "ko" + line.split()[1]
            fh.write("%s\n" % line)
            continue
        elif tag == "D":

            if path_id not in path_dict:
                continue

            mess = line.split()
            ko = mess[1]
            name = " ".join(mess[2:])

            if ko not in path_dict[path_id]:
                continue

            for p in path_dict[path_id][ko]:

                if ko == "-":
                    fh.write("D      %s\t\n" % p)
                else:
                    fh.write("D      %s\t%s %s\n" % (p, ko, name))
        else:
            fh.write("%s\n" % line)

    return output


def stat_keg(keg):
    """
    read pathway information from .keg
    :param keg: .keg file
    :return: a dict {pathway_A: {pathway_B: [proteins]}}
    """

    r = OrderedDict()
    LOG.info("reading kegg map from %r" % keg)

    path1 = ""
    path2 = ""

    for line in open(keg):
        line = line.strip()

        if not line:
            continue

        tag = line[0]

        if tag == "A":
            path1 = " ".join(line.split()[1:])
            r[path1] = OrderedDict()
            continue

        if tag == "B":
            if len(line) == 1:
                continue
            path2 = " ".join(line.split()[2:])
            r[path1][path2] = []
            continue

        if tag == "D":

            r[path1][path2].append(line.split()[1])

    return r


def plot_keg(keg_dict, out):
    """
    plot function
    :param keg_dict: see stat_keg
    :param out: output filename
    :return: 0
    """
    x = []
    y = []
    n = 1

    for path1 in keg_dict:
        _x = [n]
        _y = []
        _genes = set()
        n += 1

        for path2 in keg_dict[path1]:
            tmp = set(keg_dict[path1][path2])
            num = len(tmp)

            if not num:
                continue

            _x.append(n)
            _y.append(num)
            _genes |= tmp
            n += 1

        x += _x
        y += [len(_genes)] + _y

    y_max = max(y) * 1.1

    colors = []
    color = ["", "blue", "green", "red", "purple", "skyblue", "orange", "pink", "salmon", "gray"]
    lv = 0
    n = 1

    LOG.info("plot KEGG annotation result to %r" % out)

    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    #plt.style.use('seaborn')
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, )

    for path1 in keg_dict:

        lv += 1
        colors.append(color[lv])
        ax.text(y[n-1], n, "{0:,}".format(y[n-1]), fontsize=8, verticalalignment='center', family="Arial", )
        ax.text(y_max / -1.7, n, path1, fontsize=10, verticalalignment='center', horizontalalignment='left', family="Arial",
                color=color[lv])
        n += 1

        for path2 in keg_dict[path1]:

            tmp = set(keg_dict[path1][path2])
            num = len(tmp)

            if not num:
                continue

            ax.text(num, n, num, fontsize=8, verticalalignment='center', family="Arial",)
            ax.text(y_max / -1.8, n, path2, fontsize=8, verticalalignment='center', horizontalalignment='left', family="Arial",
                    color=color[lv])
            colors.append(color[lv])
            n += 1

    ax.barh(x, y, color=colors, alpha=0.5)

    ax.set_xlim([-y_max/100, y_max])
    ax.set_ylim([0, len(x)+1])
    plt.xticks(fontsize=8, family="Arial",)

    ax.set_yticks([])
    plt.subplots_adjust(top=0.97, left=0.35, right=0.95, bottom=0.07)

    plt.xlabel("Number of Genes", fontsize=10, family="Arial", weight="bold")

    ax.invert_yaxis()
    plt.savefig(out+".pdf")
    plt.savefig(out+".png", dpi=900)

    return out


def set_args():

    args = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description="""
create .keg file from kegg annotation result

version: %s
contact: %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args.add_argument("--keg", metavar="FILE", required=True,
                      help="KO file downloaded from KEGG, usually named 'ko00001.keg'")
    args.add_argument("--in", metavar="FILE", dest="input", required=True,
                      help="KEGG annotation result consist protein id, KO, pathway joined with '\t'")
    args.add_argument("--out", metavar="STR", default="out", help="output prefix (default: out)")
    args.add_argument("--plot", action="store_true", help="plot graph")

    return args.parse_args()


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    args = set_args()

    path_dict = cluster_protein(args.input)
    keg = output_keg(args.keg, path_dict, args.out+".keg")

    if args.plot:
        plot_keg(stat_keg(keg), args.out)


if __name__ == "__main__":
    main()

