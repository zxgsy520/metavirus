#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

if sys.getdefaultencoding() != "utf8":
    reload(sys)
    sys.setdefaultencoding("utf8")

LOG = logging.getLogger(__name__)

__version__ = "2.2.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith('#'):
            continue

        yield line.split('\t')


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fa = gzip.open(file)
    else:
        fa = open(file)

    seq = ""
    for line in fa:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            if seq:
                yield seq.split("\n", 1)
            seq = "%s\n" % line.strip(">").split()[0]
            continue
        seq += line
    yield seq.split("\n", 1)
    fa.close()


def _combime(nr, SwissProt, KEGG, COG, TIGRFAMs, Pfam, GO):

    stat_dict = {
        "nr":{},
        "Pfam":{},
        "TIGRFAMs":{},
        "COG":{},
        "KEGG":{},
        "GO":{},
        "SwissProt":{},
        "name":{}
        }

    for record in sorted(read_tsv(nr)):
        stat_dict["nr"][record[0]] = record[5].replace(",", "%2C").split(";")[0]

    for record in sorted(read_tsv(SwissProt)):
        product = record[5]
        stat_dict["SwissProt"][record[0]] = product
        if "rotein" in product:
            gene = product.split("rotein")[-1].strip()
        else:
            gene = product.split()[-1]
        if len(gene) <=2:
            continue
        if record[0] not in stat_dict["name"]:
            stat_dict["name"][record[0]] = gene
        else:
            if stat_dict["name"][record[0]] == "-":
                stat_dict["name"][record[0]] = gene
            else:
                stat_dict["name"][record[0]] += ",%s" % gene

    for record in sorted(read_tsv(Pfam)):
        stat_dict["Pfam"][record[0]] = [record[2], record[5].replace(",", "%2C").split(";")[0]]

    for record in sorted(read_tsv(TIGRFAMs)):
        stat_dict["TIGRFAMs"][record[0]] = [record[2], record[5].replace(",", "%2C").split(";")[0]]
        if record[0] not in stat_dict["name"]:
            stat_dict["name"][record[0]] = record[4]
        else:
            if stat_dict["name"][record[0]] == "-":
                stat_dict["name"][record[0]] = record[4]
            else:
                stat_dict["name"][record[0]] += ",%s" % record[4]

    for record in sorted(read_tsv(COG)):
        stat_dict["COG"][record[0]] = [record[2], record[5].replace(",", "%2C").split(";")[0], record[3]]
        if record[0] not in stat_dict["name"]:
            stat_dict["name"][record[0]] = record[4]
        else:
            if stat_dict["name"][record[0]] == "-":
                stat_dict["name"][record[0]] = record[4]
            else:
                stat_dict["name"][record[0]] += ",%s" % record[4]

    for record in sorted(read_tsv(KEGG)):
        stat_dict["KEGG"][record[0]] = [record[2], record[5].replace(",", "%2C").split(";")[0], record[3]]
        if record[0] not in stat_dict["name"]:
            stat_dict["name"][record[0]] = record[4]
        else:
            if stat_dict["name"][record[0]] == "-":
                stat_dict["name"][record[0]] = record[4]
            else:
                stat_dict["name"][record[0]] += ",%s" % record[4]

    for record in sorted(read_tsv(GO)):

        goid = ""
        godes = ""

        for i in record[1:]:
            if i == "-":
                continue
            for j in i.split(";"):
                n, d = j.split(" - ")
                goid += "%s;" % n
                godes += "%s;" % d.replace(",", "%2C").replace("(", " - ").replace(")", "")

        stat_dict["GO"][record[0]] = [goid.strip(";"), godes.strip(";")]

    return stat_dict


def merge_anno(file, nr, SwissProt, KEGG, COG, TIGRFAMs, Pfam, GO):

    print("#Gene_Id\tGene_Length(bp)\tGene_Name\tSwissProt_Description\t"\
          "nr_Description\tPfam_Id\tPfam_Description\t"\
          "TIGRFAMs_Id\tTIGRFAMs_Description\tCOG_Id\tCOG_Description\tCOG_Type\t"\
          "KO_Id\tKO_Description\tPathway\tGO_Id\tGO_Description")

    stat_dict = _combime(nr, SwissProt, KEGG, COG, TIGRFAMs, Pfam, GO)
    pep_num = 0
    for seqid, seq in read_fasta(file):
        temp = [seqid, str(len(seq))]
        pep_num += 1

        if seqid in stat_dict["name"]:
            temp.append(stat_dict["name"][seqid])
        else:
            temp.append("-")

        if seqid in stat_dict["SwissProt"]:
            temp.append(stat_dict["SwissProt"][seqid])
        else:
            temp.append("-")

        if seqid in stat_dict["nr"]:
            temp.append(stat_dict["nr"][seqid])
        else:
            temp.append("-")

        if seqid in stat_dict["Pfam"]:
            temp += stat_dict["Pfam"][seqid]
        else:
            temp += ["-", "-"]

        if seqid in stat_dict["TIGRFAMs"]:
            temp += stat_dict["TIGRFAMs"][seqid]
        else:
            temp += ["-", "-"]

        if seqid in stat_dict["COG"]:
            temp += stat_dict["COG"][seqid]
        else:
            temp += ["-", "-", "-"]

        if seqid in stat_dict["KEGG"]:
            temp += stat_dict["KEGG"][seqid]
        else:
            temp += ["-", "-", "-"]

        if seqid in stat_dict["GO"]:
            temp += stat_dict["GO"][seqid]
        else:
            temp += ["-", "-"]

        print("\t".join(temp))

    fh = open("function_summary.tsv", "w")
    _all = set()
    _one = set()

    for k, v in stat_dict.items():
        v = set(v.keys())
        if k == "name":
            continue
        fh.write("%s\t%s\t%.2f\n" % (k, len(v), len(v)*100.0/pep_num))

        if not _all:
            _all = v
        _all = _all & v

        _one = _one | v

    fh.write("""\
all databases\t%s\t%.2f
at least one databases\t%s\t%.2f
""" % (len(_all), len(_all)*100.0/pep_num, len(_one), len(_one)*100.0/pep_num))

    fh.close()

    return 0


def add_help_args(parser):

    parser.add_argument("protein", metavar='FILE', type=str,
        help="Protein file for input gene annotations.")
    parser.add_argument("-sp", "--swissprot", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the SwissProt database.")
    parser.add_argument("-nr", "--nr", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the nr(refseq) database.")
    parser.add_argument("-kg", "--KEGG", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the kegg database.")
    parser.add_argument("-cog", "--COG", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the kegg database.")
    parser.add_argument("-tf", "--TIGRFAMs", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the TIGRFAMs database.")
    parser.add_argument("-pf", "--Pfam", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the Pfam database.")
    parser.add_argument("-go", "--GO", metavar='FILE', type=str, required=True,
        help="Input the annotation results of the GO database.")

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
    merge_protein_anno.py  Statistical summary of each database.

attention:
    merge_protein_anno.py protein.fasta -nr nr.tsv(refseq.tsv) -kg KEGG.tsv -cog COG.tsv >stat_anno.tsv

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    merge_anno(args.protein, args.nr, args.swissprot, args.KEGG, args.COG, args.TIGRFAMs, args.Pfam, args.GO)


if __name__ == "__main__":

    main()
