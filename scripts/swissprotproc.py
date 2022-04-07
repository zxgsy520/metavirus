#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging

import collections
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


def split_attr(attributes):

    r = collections.OrderedDict()
    contents = attributes.split()
    desc = ""
    temp = ""

    for content in contents:
        if "=" not in content:
            if not temp:
                temp = content
            else:
                temp += " %s" % content
        else:
            if "=" not in temp:
                desc = temp
                temp = content
            else:
                tag, value = temp.split("=", 1)
                r[tag] = value
                temp = content

    if "=" in temp:
        tag, value = temp.split("=", 1)
        r[tag] = value

    return r, desc


def swissprotproc(file):

    print("#qseqid\tsseqid\tdbxref\tclass\tname\tdesc\tqstart\tqend\tevalue\tscore\tnote")
    for record in read_tsv(file, sep="\t"):
        data, desc = split_attr(record[5].split(" ", 1)[-1])

        if "GN" in data:
            try:
                k = float(data["GN"])
                data["GN"] = "-"
            except:
                pass
        else:
            data["GN"] = "-"

        print("\t".join([record[0], record[1], "SwissProt:%s" % record[1], "-", data["GN"],
                         desc, record[2], record[3], record[6], record[7], "Taxon:%s" % data["OS"]]))

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input swissprot annotation results, swissprot.out")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
Exmple:
        python swissprotproc.py swissprot.out >swissprot.tsv

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    swissprotproc(args.input)


if __name__ == "__main__":

    main()
