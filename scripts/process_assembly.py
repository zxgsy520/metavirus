#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
from Bio import SeqIO

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Junpeng Fan, Xingguo Zhang",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def process_assembly(fasta, _format, min_length, organism, strain, gcode):

    n = 1

    if _format:
        records = sorted(SeqIO.parse(fasta, "fasta"), key=lambda k: len(k), reverse=True)
    else:
        records = SeqIO.parse(fasta, "fasta")

    for record in records:

        if len(record.seq) < min_length:
            continue
        record.seq = record.seq.upper()
        record.description = record.description.split(None, 1)
        
        if len(record.description)>=2:
            record.description = record.description[1]
        else:
            record.description = ""

        if _format:
            record.id = _format % n
        if organism:
            record.description += " [organism=%s]" % organism
        if strain:
            record.description += " [strain=%s]" % strain
        if gcode:
            record.description += " [gcode=%s]" % gcode

        print(record.format("fasta"))
        n += 1


def add_args(parser):

    parser.add_argument("fasta", help="")
    parser.add_argument("--format", help="contig name format")
    parser.add_argument("--min_length", type=int, default=500, help="")
    parser.add_argument("--organism", help="organism name")
    parser.add_argument("--strain", help="strain")
    parser.add_argument("--gcode", type=int, help="genetic code")

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
    process_assembly(args.fasta, args.format, args.min_length, args.organism, args.strain, args.gcode)


if __name__ == "__main__":
    main()

