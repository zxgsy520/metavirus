#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Junpeng Fan",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def read_tsv(file):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split("\t")


def read_obo(file):

    term = False
    r = []
    for line in open(file):
        line = line.strip()

        if not line:
            term = False
            continue

        if line == "[Term]":
            term = True
            if r:
                yield r
                r = []

            continue

        if term is True:
            r.append(line)

    if r:
        yield r


def get_term(file):
    r = {}

    for record in read_obo(file):
        term_id = record[0][4:]
        name = record[1][6:]
        name_space = record[2][11:]

        r[term_id] = [name, name_space]

    return r


def go_proc(anno, obo):

    term_dict = get_term(obo)

    print("#qseqid\tbiological_process\tcellular_component\tmolecular_function")

    for record in sorted(read_tsv(anno)):
        name = record[0]
        go = record[1:]

        _tmp = {
            "biological_process": [],
            "molecular_function": [],
            "cellular_component": []
        }

        for g in go:

            if g in term_dict:
                _name, _name_space = term_dict[g]
                _tmp[_name_space].append("%s - %s" % (g, _name))

        _mess = []
        for k in ["biological_process", "cellular_component", "molecular_function"]:

            if not _tmp[k]:
                _mess.append("-")
            else:
                _mess.append(";".join(_tmp[k]))

        print("%s\t%s" % (name, "\t".join(_mess)))


def add_args(parser):

    parser.add_argument("--go", required=True, help="")
    parser.add_argument("--go_obo", required=True, help="")

    return parser


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""


version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_args(parser)
    args = parser.parse_args()

    go_proc(args.go, args.go_obo)


if __name__ == "__main__":
    main()

