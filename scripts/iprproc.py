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

        if not line:
            continue

        yield line.split("\t")


def read_ipr(file):

    for record in read_tsv(file):
        qseqid = record[0]
        _type, _id, desc, start, end, _evalue = record[3:9]

        if len(record) == 14:
            go = record[13]
        else:
            go = "-"

        yield qseqid, _type, _id, desc, start, end, _evalue, go


def _join(p):

    r = []

    for record in sorted(p, key=lambda d: float(d[6])):
        if not r:
            r = record
            continue

        r[4] += ";"+record[4]
        r[5] += ";" + record[5]
        r[6] += ";" + record[6]

    return r


def _combine(file):

    p = []

    for record in sorted(read_ipr(file), key=lambda d: (d[0], d[1], d[2])):

        if not p:
            p = [record]
            continue

        if p[0][0] == record[0] and p[0][1] == record[1] and p[0][2] == record[2]:
            p.append(record)
        else:
            yield _join(p)
            p = [record]

    if p:
        yield _join(p)


def read_tigrfam(file):

    r = {}

    for record in read_tsv(file):
        r[record[0]] = record[1:]

    return r


def iprproc(ipr, tigrfam, prefix):

    tigrfam_dict = read_tigrfam(tigrfam)

    tigrfams_out = open("%s.TIGRFAMs.tsv" % prefix, "w")
    tigrfams_out.write("#qseqid\tsseqid\tdbxref\tclass\tname\tdesc\tqstart\tqend\tevalue\tscore\tnote\n")
    pfam_out = open("%s.Pfam.tsv" % prefix, "w")
    pfam_out.write("#qseqid\tsseqid\tdbxref\tclass\tname\tdesc\tqstart\tqend\tevalue\tscore\tnote\n")
    go_out = open("%s.WEGO.txt" % prefix, "w")

    go_dict = {}

    print("#qseqid\tsseqid\tdbxref\tclass\tname\tdesc\tqstart\tqend\tevalue\tscore\tnote")

    for qseqid, _type, _id, desc, start, end, _evalue, go in sorted(read_ipr(ipr), key=lambda d: (d[0], d[1], float(d[6].split(";")[0]))):

        if _type == "TIGRFAM":
            dbxref = "TIGRFAMs:%s" % _id
            _class = "-"
            desc, name, ec, go = tigrfam_dict[_id]

            if ec == "-":
                note = "-"
            else:
                note = "EC:%s" % ec
            tigrfams_out.write("{qseqid}\t{_id}\t{dbxref}\t{_class}\t{name}\t"
                               "{desc}\t{start}\t{end}\t{_evalue}\t-\t{note}\n".format(**locals()))
        elif _type == "Pfam":
            dbxref = "Pfam:%s" % _id
            _class = note = name = "-"
            pfam_out.write(
                "{qseqid}\t{_id}\t{dbxref}\t{_class}\t{name}\t{desc}\t{start}\t{end}\t{_evalue}\t-\t{note}\n".format(
                    **locals()))
        else:
            dbxref = "%s:%s" % (_type, _id)
            _class = note = name = "-"
            print("{qseqid}\t{_id}\t{dbxref}\t{_class}\t{name}\t{desc}\t{start}\t{end}\t{_evalue}\t-\t{note}".format(**locals()))

        if go != "-":
            if qseqid not in go_dict:
                go_dict[qseqid] = []

            go_dict[qseqid] += go.replace(";", "|").split("|")

    for k, v in go_dict.items():
        if v:
            go_out.write("%s\t%s\n" % (k, "\t".join(set(v))))

    tigrfams_out.close()
    pfam_out.close()
    go_out.close()


def filter_ipr_result(file):

    pqseqid = ptype = ""

    for record in sorted(read_ipr(file), key=lambda d: (d[0], d[3], d[8])):

        qseqid, _type, _id, _name, _evalue = record[0], record[3], record[4], record[5], record[8]

        if len(record) == 13:
            go = record[13]
        else:
            go = ""

        if not pqseqid:
            pqseqid = qseqid
            ptype = _type
            continue

        if qseqid == pqseqid and _type == ptype:
            continue
        else:
            yield qseqid, _type, _id, _evalue, _name, go
            pqseqid = qseqid
            ptype = _type

    if not pqseqid:
        yield qseqid, _type, _id,  _evalue, _name, go


def set_args():

    args = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description="""
description:

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args.add_argument("ipr", help="")
    args.add_argument("--tigrfams", help="")
    args.add_argument("--prefix", help="")

    return args.parse_args()


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    args = set_args()

    iprproc(args.ipr, args.tigrfams, args.prefix)


if __name__ == "__main__":
    main()

