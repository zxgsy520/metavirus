#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse
import collections

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


HEAD = ["##gff-version", "##source-version", "##date:", "# Sequence file name:",
        "# Model file name:", "# RBS: false"]


def check_head(line):

    r = False

    for i in HEAD:
        if i in line:
            r = True
            break
    return r


def read_mgm(file):

    if file.endswith(".gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode("utf-8")

        line = line.strip()
        if not line or check_head(line):
            continue
        yield line
    fp.close()


def split_attr(attributes):

    r = collections.OrderedDict()
    contents = attributes.split(";")

    for content in contents:
        if not content:
            continue
        if "=" not in content:
            print("%r is not a good formated attribute: no tag!")
            continue
        tag, value = content.split("=", 1)
        r[tag] = value

    return r


def to_string(attributes):

    attr = []

    for key, value in attributes.items():
        if key in "ID":
            attr.insert(0, '%s=%s' % (key, value))
        elif key in "Name":
            attr.insert(1, '%s=%s' % (key, value))
        elif key in "Parent":
            attr.insert(2, '%s=%s' % (key, value))
        else:
            attr.append('%s=%s' % (key, value))

    return ';'.join(attr)


def process_mgm(file, prefix="MNGS"):

    r = {}
    protein = []
    gene = []
    seq = ""
    n = 0

    print("##gff-version 3")
    for line in read_mgm(file):
        if "# Model information:" in line:
            n = 0
            continue
        if not line.startswith("#"):
            n += 1
            line = line.split('\t')
            attr = split_attr(line[-1])
            #print(line)
            gengid = attr["gene_id"]
            newid = "%s.%s.%s" % (prefix, line[0], n)
            r[gengid] = newid
            attr["Parent"] = newid
            attr["gene_id"] = newid
            line[8] = to_string(attr)
            print("\t".join(line))
            continue
        elif line.startswith("##Protein") or line.startswith("##DNA "):
            seq = ""
            seqid = line.split()[-1]
            continue
        elif line.startswith("##end-DNA"):
            gene.append(">%s\n%s" % (r[seqid], seq))
            continue
        elif line.startswith("##end-Protein"):
            protein.append(">%s\n%s" % (r[seqid], seq))
            continue
        else:
            seq += line.strip().strip("##")

    return protein, gene


def out_seq(protein, gene, prefix="MNGS"):

    fp = open("%s.protein.fasta" % prefix, "w")
    for line in protein:
        fp.write("%s\n" % line)
    fp.close()

    fg = open("%s.gene.fasta" % prefix, "w")
    for line in gene:
        fg.write("%s\n" % line)
    fg.close()


def add_hlep_args(parser):

    parser.add_argument("input", metavar='FILE', type=str,
        help="Input GeneMark annotation results,(.mgm, gff)")
    parser.add_argument("-p", "--prefix", metavar='FILE', type=str, default="MNGS",
        help="Set the prefix of the output file, the sample name, default=NGS")

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
    mgm2gff.py Extract gff files, genome sequences and protein sequences from the prediction results of GeneMark.

attention:
    mgm2gff.py name.mgm -p name.gff >name.gff

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()
    protein, gene = process_mgm(args.input, args.prefix)

    out_seq(protein, gene, args.prefix)


if __name__ == "__main__":

    main()
