#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
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


def plot_pep(lengths, window, prefix):

    x = [(i+0.5)*window for i in range(int(max(lengths)/window)+1)]
    y = [0 for _ in x]

    for i in lengths:
        y[int(i/window)] += 1

    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    #plt.style.use('seaborn')
    fig, ax = plt.subplots(figsize=(6.5, 4.5), )
    plt.subplots_adjust(top=0.95, left=0.10, right=0.95, bottom=0.10)
    ax.bar(x, y, width=window, linewidth=0.5, edgecolor="white", )

    miny = window*-0.1
    if miny < 0:
        miny = 0

    plt.ylim([miny, plt.ylim()[1]])
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.xlabel("Protein length (aa)", fontsize=10, weight="bold")
    plt.ylabel("Number", weight="bold", fontsize=10,)
    plt.savefig("%s.protein_length.pdf" % prefix)
    plt.savefig("%s.protein_length.png" % prefix, dpi=900)

    return 0

def plot_cds2pep(file, prefix, window=100):

    lengths = []

    for seqid, seq in read_fasta(file):
        if len(seq) % 3 != 0:
            continue
        lengths.append(len(seq)/3)

    plot_pep(lengths, window, prefix)

    return 0


def add_help_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input sequence file, format(fasta, fasta.gz)")
    parser.add_argument("-p", "--prefix", metavar="FILE", type=str, required=True,
        help="output file prefix.")

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
    plot_cds2pep.py -- Count the length of genes and draw a graph of protein length.
attention:
    plot_cds2pep.py 1.gene.fasta -p txt
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()
    plot_cds2pep(args.input, args.prefix)


if __name__ == "__main__":
    main()
