#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging

LOG = logging.getLogger(__name__)

__version__ = "1.1.1"
__author__ = ("Junpeng Fan, Xingguo Zhang",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def read_tsv(file, sep=None):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def refseqproc(file):

    print("#qseqid\tsseqid\tdbxref\tclass\tname\tdesc\tqstart\tqend\tevalue\tscore\tnote")
    for record in read_tsv(file, sep="\t"):
        if len(record[4].split(">")) <=1:
            annotate = record[4].strip()
        else:
            for i in record[4].split(">"):
                annotate = i.strip()
                if "hypothetical protein" in i:
                    continue
                annotate = i.strip()
                break
                
        desc = "[".join(annotate.split("[")[:-1])
        taxon = annotate.split("[")[-1].rstrip("]")
        gene = "-"
        if "protein" in desc:
            gene = desc.split("protein")[-1]
            if len(gene) >= 2:
               gene = gene.split()[0]
            if len(gene) <= 2:
                gene = "-"
            if gene[0].isupper():
                pass
            else:
                gene = "-"
        print("\t".join([record[0], record[1], "Refseq:%s" % record[1], "-", gene,
                         desc, record[2], record[3], record[5], record[6], "Taxon:%s" % taxon]))


def set_args():
    args = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description="""
description:

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args.add_argument("file", help="")

    return args.parse_args()


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    args = set_args()

    refseqproc(args.file)



if __name__ == "__main__":
    main()
