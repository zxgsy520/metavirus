#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.0.1"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


GET_NEXT_LEVEL = {"k": ["p__", 1],
    "p": ["c__", 2],
    "c": ["o__", 3],
    "o": ["f__", 4],
    "f": ["g__", 5],
    "g": ["s__", 6],
    "s": ["ITF", 7]
    }


def read_tsv(file, sep="\t"):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue

        yield line.split(sep)

    fh.close()


def split_tax(tax):

    r = OrderedDict()

    for i in tax.split("|"):
        level, value = i.split("__", 1)
        if (level == "k") and (level in r):
            continue
        if "s" == level:
            value = value.split(".")[0]
        r[level] = value

    return r


def tax2str(taxdict):

    r = []

    for key, value in taxdict.items():
        r.append('%s__%s' % (key, value))

    return '|'.join(r)


def str_list2float(otulist):

    r = []

    for i in otulist:
        r.append(float(i))

    return r


def addlist(list1, list2):

    r = []
    n = 0
    for i in list2:
        r.append(list1[n]+i)
        n += 1

    return r


def read_mpa(file, kingdom="Bacteria", level="g"):

    nlevel = ""
    if level in GET_NEXT_LEVEL:
        nlevel, rank = GET_NEXT_LEVEL[level]
    level = "%s__" % level

    data = {}
    samples = []
    total_abunds = []

    for line in read_tsv(file, sep="\t"):
        if line[0].startswith("#"):
            samples = line[1::]
            continue
        if kingdom=="Bacteria":
            line[0] = line[0].replace("Archaea", "Arch%s" % kingdom, 1)
        if kingdom not in line[0]:
            continue
        if level not in line[0]:
            continue
        if kingdom=="Bacteria":
            line[0] = line[0].replace("Arch%s" % kingdom, "Archaea", 1)
        tax = split_tax(line[0])
        if (nlevel in line[0]) or len(tax) > rank:
            continue

        taxid = tax2str(tax)
        #taxid = taxid.split(nlevel)[0].strip("|")
        abunds = str_list2float(line[1::])
        if len(total_abunds) == 0:
            total_abunds = abunds
        else:
            total_abunds = addlist(total_abunds, abunds)
        if taxid not in data:
            data[taxid] = abunds
        else:
            data[taxid] = addlist(data[taxid], abunds)

    return samples, total_abunds, data


def process_mean_abundance(abunds, otulist):

    r = []
    meanabund = sum(abunds)*1.0/len(abunds)

    n = 0
    for i in otulist:
        r.append("{:.2f}".format(i*1.0*meanabund/abunds[n]))
        n += 1

    return r


def process_uniform_abundance(abunds, otulist):

    r = []
    n = 0

    for i in otulist:
        r.append("{:.6f}".format(i*100.0/abunds[n]))
        n += 1

    return r


def tlist_str(otulist):

    r = []
    for i in otulist:
        r.append(str(i))

    return r


def get_relative_abundance(file, kingdom="Bacteria", level="g", model="mean", cut=False):

    samples, abunds, data =  read_mpa(file, kingdom, level)

    print("Tax Id\t%s" % "\t".join(samples))
    for taxid, line in sorted(data.items(), key=lambda x: sum(x[1])/len(x[1]), reverse=True):
        if model=="mean":
            line = process_mean_abundance(abunds, line)
        elif model=="uniform":
            line = process_uniform_abundance(abunds, line)
        else:
            line = tlist_str(line)
        if cut:
            taxid = taxid.split("|")[-1]
        print("%s\t%s" % (taxid, "\t".join(line)))

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input the merged abundance table, merge_mpa2otu.tsv")
    parser.add_argument("-l", "--level", metavar='STR', type=str, default="",
        choices=["p", "c", "o", "f", "g", "s", ""],
        help="Input the displayed species level. choices=[p,c,o,f,g,s], default=all.")
    parser.add_argument("-k", "--kingdom", metavar='STR', type=str, default="",
        choices=["Bacteria", "Eukaryota", "Viruses"],
        help="Choose the object of the study. choices=[Bacteria, Eukaryota, Viruses], default=all.")
    parser.add_argument("-m", "--model", metavar="STR", type=str,  default="uniform",
        choices=["mean", "uniform", "original"],
        help="Abundance conversion method. choices=[mean, uniform, original], default=uniform.")
    parser.add_argument("--cut", action="store_true", default=False,
        help="Whether to trim taxid, default=False.")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
For exmple:
        get_relative_abundance.py merge_mpa2otu.tsv > relative_abundance.tsv

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    get_relative_abundance(args.input, args.kingdom, args.level, args.model, args.cut)


if __name__ == "__main__":

    main()
