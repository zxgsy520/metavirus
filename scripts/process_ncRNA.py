#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
from collections import OrderedDict


LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Junpeng Fan, Xingguo Zhang",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def read_tsv(file, sep=None):

    for line in open(file):
        line = line.strip()

        if not line:
            continue

        yield line.split(sep)


def read_trna(file):

    r = []

    for record in read_tsv(file, "\t"):

        if record[0].startswith("#VERSION"):
            version = record[0].split()[1]
            continue
        elif record[0].startswith("#ARGS"):
            args = record[0].split()[1]
            continue
        elif record[0].split()[0] in "Sequence|Name|--------":
            continue
        else:
            pass

        seq, start, end, _type, _anticoden, intro_start, intro_end, score = record[0], record[2], record[3], record[4], record[5], record[6], record[7], record[8]

        if _type == "Undet":
            continue

        start, end, intro_start, intro_end = map(int, (start, end, intro_start, intro_end))

        if start > end:
            start, end = end, start
            strand = "-"
            anticoden = "pos:complement(%s..%s),aa:%s,seq:%s" % (start, end, _type, _anticoden)
        else:
            strand = "+"
            anticoden = "pos:%s..%s,aa:%s,seq:%s" % (start, end, _type, _anticoden)

        attributes = OrderedDict([
            ("inference", ["COORDINATES: profile:%s" % version, ]),
            ("product", ["tRNA-%s" % _type]),
            ("anticodon", [anticoden])]
        )

        feature = [seq.strip(), version.strip(), "tRNA", start, end, score, strand, ".", attributes]
        r.append(feature)

    return r


def read_rrna(file):

    r = []

    for record in read_tsv(file, "\t"):

        if record[0].startswith("#VERSION"):
            version = record[0].split()[1]
            continue
        elif record[0].startswith("#ARGS"):
            args = record[0].split()[1]
            continue
        elif record[0].startswith("#"):
            continue
        else:
            pass

        seq, _type, start, end, score, strand, _attr = record[0], record[2], record[3], record[4], record[5], record[6], record[8]

        if "barrnap" in record[1]:
            version = record[1]

        start, end = int(start), int(end)

        if "5s_rRNA" in _attr:
            product = "5S ribosomal RNA"
        elif "16s_rRNA" in _attr:
            product = "16S ribosomal RNA"
        elif "23s_rRNA" in _attr:
            product = "23S ribosomal RNA"
        elif "8s_rRNA" in _attr:
            product = "5S ribosomal RNA"
        elif "18s_rRNA" in _attr:
            product = "18S ribosomal RNA"
        elif "28s_rRNA" in _attr:
            product = "28S ribosomal RNA"
        else:
            product = _attr

        attributes = OrderedDict([
            ("inference", ["COORDINATES: profile:%s" % version, ]),
            ("product", [product])]
        )

        feature = [seq, version, "rRNA", start, end, score, strand, ".", attributes]
        r.append(feature)

    return r


def discover_rna(start, end, product):

    start, end = int(start), int(end)
    inference = product

    if start > end:
        start, end = int(end), int(start)
    rna_len = end-start+1

    if "ribosomal" in product and "subunit" in product:
        if "Bacterial" in product or "Archaeal" in product:
            if rna_len<=500:
                product = "5S ribosomal RNA"
            elif rna_len>=2500:
                product = "23S ribosomal RNA"
            else:
                product = "16S ribosomal RNA"
        elif "Eukaryotic" in product:
            if rna_len<=140:
                product = "5S ribosomal RNA"
            elif rna_len<=500:
                product = "5.8S ribosomal RNA"
            elif rna_len>=3000:
                product = "28S ribosomal RNA"
            else:
                product = "18S ribosomal RNA"
    else:
        pass
    return product, inference


def read_rfam(file, family):

    r = []

    family_dict = {}

    for record in read_tsv(family, sep="\t"):
        family_dict[record[0]] = record[1:]

    for record in read_tsv(file):

        if record[0].startswith("#VERSION"):
            version = record[1]
            continue
        elif record[0].startswith("#ARGS"):
            args = record[1]
            continue
        elif record[0].startswith("#"):
            continue
        else:
            pass

        name, _id, seq, start, end, strand, score, evalue, olp = record[1], record[2], record[3], record[9], record[10], record[11], record[16], record[17], record[19]
        start, end = int(start), int(end)
        if start > end:
            start, end = end, start

        if olp == "=":
            continue

        product, _type, clen = family_dict[_id][1:]

        if end-start+1 < int(clen)*0.8:
            sys.stderr.write("RNA in %s:%s-%s len: %s match %s %s region < 0.8*%s\n" % (seq, start, end, end-start+1,_id, name, clen))
            continue

        _class = ""

        if _type.startswith("Gene; rRNA"):
            _type = "rRNA"
            product, inference = discover_rna(start, end, product)
#            if product!=inference:
#                product = "%s;inference=Rfam:%s" % (product, inference)
        elif _type.startswith("Gene; tRNA"):
            _type = "tRNA"
        elif _type.startswith("Gene;"):
            _type = "ncRNA"
            _class = ["other"]
        elif _type.startswith("Cis-reg;"):
            _type = "regulatory"
        else:
            continue

        if _class:
            attributes = OrderedDict([
                ("ncRNA_class", _class),
                ("inference", ["COORDINATES: profile:%s" % version, "COORDINATES: nucleotide motif:Rfam:14.0:%s" % _id]),
                ("product", [product]),
            ])
        else:
            attributes = OrderedDict([
                ("inference", ["COORDINATES: profile:%s" % version, "COORDINATES: nucleotide motif:Rfam:14.0:%s" % _id]),
                ("product", [product]),
            ])

        feature = [seq, version, _type, start, end, score, strand, ".", attributes]
        r.append(feature)

    return r


def read_crispr(file):

    r = []

    for record in read_tsv(file, sep="\t"):
        r.append(record)

    return r


def _overlap(q, s):

    s1, e1 = q[3], q[4]
    s2, e2 = s[3], s[4]

    assert s2 >= s1

    if (e1-s2+1) > min((e1-s1+1, e2-s2+1))*0.5:
        return 1
    else:
        return 0


def _deduplicate(trna, rrna, rfam, rfamily):

    r = []
    p = []

    for g in sorted(read_trna(trna)+read_rrna(rrna)+read_rfam(rfam, rfamily), key=lambda i: (i[0], i[3])):

        if not p:
            p = g
            continue

        if g[0] != p[0]:
            r.append(p)
            p = g
            continue

        if _overlap(p, g):
            if p[2] == g[2]:
                if g[1].startswith("RNAmmer") or g[1].startswith("tRNAscan"):
                    sys.stderr.write("%r overlap %r \n" % (p, g))
                    p = g

                continue
            else:
                if g[4] - g[3] > p[4] - p[3]:
                    r.append(p)
                    p = g
                else:
                    r.append(g)

                continue

        r.append(p)
        p = g

    if p:
        r.append(p)

    return r


def process_ncrna(trna, rrna, rfam, rfamily):

    print("##gff-version 3")

    for g in sorted(_deduplicate(trna, rrna, rfam, rfamily), key=lambda d: (d[0], int(d[3]))):

        if isinstance(g[-1], OrderedDict):
            print("%s\t%s" % ("\t".join(map(str, g[:-1])), ";".join(["%s=%s" % (k, ",".join(v)) for k, v in g[-1].items()])))
        else:
            print("\t".join(g))

    return 0


def set_args():
    args = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description="""
description:

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args.add_argument("--tRNA", help="")
    args.add_argument("--rRNA", help="")
    args.add_argument("--rfam", help="")
    args.add_argument("--rfamily", help="")

    return args.parse_args()


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    args = set_args()
    process_ncrna(args.tRNA, args.rRNA, args.rfam, args.rfamily)

if __name__ == "__main__":
    main()
