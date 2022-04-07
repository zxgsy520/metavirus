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



def read_tsv(file, sep="\t"):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)

    fh.close()


def split_tax(tax):

    r = OrderedDict()

    for i in tax.split("|"):
        level, value = i.split("__", 1)
        if (level == "k") and (level in r):
            continue
        r[level] = value

    return r


def tax2str(taxdict):

    r = []

    for key, value in taxdict.items():
        r.append('%s__%s' % (key, value))

    return '|'.join(r)


def generate_list(length):

    r = []

    for i in range(length):
        r.append(0)

    return r


def read_stat_mpa(files):

    samples = []
    data = OrderedDict()
    abunds = []
    n = 0

    for file in files:
        n += 1
        sample = file.split("/")[-1].split(".")[0]
        samples.append(sample)

        abund = 0
        for line in read_tsv(file):
            if "k__Fungi" in line[0]:
                line[0] = line[0].split("|", 1)[-1]
            if "__" in line[0]:
                tax_dict = split_tax(line[0])
            else:
                tax_dict = {"": line[0]}
            line[0] = tax2str(tax_dict)
            if len(tax_dict) == 1:
                abund += int(line[1])
            if line[0] not in data:
                data[line[0]] = generate_list(len(files))
            data[line[0]][n-1] = int(line[1])
        abunds.append(abund)

    return samples, abunds, data


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


def merge_mpa2otu(files, model="mean"):

    samples, abunds, data = read_stat_mpa(files)

    print("#ID\t%s" % "\t".join(samples))
    for otuid, line in data.items(): #sorted(data.items(), key=lambda x: (x[0], -min(x[1]))): #reverse=True
        if model=="mean":
            line = process_mean_abundance(abunds, line)
        elif model=="uniform":
            line = process_uniform_abundance(abunds, line)
        else:
            line = tlist_str(line)
        print("%s\t%s" % (otuid, "\t".join(line)))

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", nargs='+', metavar="FILE", type=str,
        help="Input the abundance statistics result file of each sample, *.kreport2mpa.report")
    parser.add_argument("-m", "--model", metavar="STR", type=str, choices=["mean", "uniform", "original"], default="original",
        help="Abundance conversion method,choices=[mean, uniform, original].default=original.")

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
        merge_mpa2otu.py *.kreport2mpa.report >all.otu.tsv

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    merge_mpa2otu(args.input, args.model)


if __name__ == "__main__":

    main()
