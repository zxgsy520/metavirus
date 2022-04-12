#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

import math
from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.2.3"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    for line in open(file):
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue

        yield line.split(sep)


def addlist(list1, list2):

    r = []
    n = 0
    for i in list2:
        r.append(list1[n]+i)
        n += 1

    return r


def str_list2float(otulist):

    r = []

    for i in otulist:
        r.append(float(i))

    return r


def read_abundance(file, top=20):

    data = OrderedDict()
    samples = []
    total_abunds = []

    n = 0
    for line in read_tsv(file, "\t"):
        n += 1
        if n ==1:
            samples = line[1::]
            continue
        taxid = line[0].split("|")[-1]
        abunds = str_list2float(line[1::])
        if sum(abunds) <= 0:
            continue
        if len(total_abunds) == 0:
            total_abunds = abunds
        else:
            total_abunds = addlist(total_abunds, abunds)
        if taxid not in data:
            data[taxid] = abunds
        else:
            data[taxid] = addlist(data[taxid], abunds)
        if n >= top+1:
            break

    return samples, total_abunds, data


def process_uniform_abundance(abunds, otulist):

    r = []
    n = 0

    for i in otulist:
        if abunds[n] == 0:
            r.append("0")
        else:
            r.append("{:.6f}".format(i*100.0/abunds[n]))
        n += 1

    return r


def read_group(file):

    data = {}

    for line in read_tsv(file, "\t"):
        if line[0].startswith("#") or line[0]=="sample":
            continue
        data[line[0]] = line[1]

    return data


def process_zscore(otulist):

    total = sum(otulist)
    meanotu = total*1.0/len(otulist)
    sd = 0

    for i in otulist:
        sd += math.pow((i-meanotu), 2)
    sd = math.sqrt(sd/len(otulist))

    r = []
    for i in otulist:
        r.append("{:.6f}".format((i-meanotu)/sd))

    return r


def process_log(otulist):

    r = []

    for i in otulist:
        if i <= 0:
            i = 1e-6
        else:
            i = i*1e4
        r.append("{:.6f}".format(math.log10(i)))

    return r


def cut_relative_abundance(file, group, top=20, model="zscore"):

    samples, abunds, data =  read_abundance(file, top)
    if group:
        group_dict = read_group(group)
        result = OrderedDict()
        groups = list(set(list(group_dict.values())))
        print("Taxid\t%s" % "\t".join(groups))
    else:
        group_dict = {}
        result = {}
        print("Taxid\t%s" % "\t".join(samples))

    for taxid, line in data.items():
        line = process_uniform_abundance(abunds, line)
        if group:
            result[taxid] = {}
            for i in range(len(samples)):
                sample = samples[i]
                group = group_dict[str(sample)]
                if group not in result[taxid]:
                    result[taxid][group] = []
                result[taxid][group].append(float(line[i]))
        else:
            if "zscore" == model:
                line = process_zscore(str_list2float(line))
            elif "log" == model:
                line = process_log(str_list2float(line))
            else:
                pass
            print("%s\t%s" % (taxid, "\t".join(line)))

    for taxid in result:
        temp = result[taxid]
        abunds = []
        for i in groups:
            abunds.append("{:.6f}".format(sum(temp[i])/len(temp[i])))
        if "zscore" == model:
            abunds = process_zscore(str_list2float(abunds))
        elif "log" == model:
            abunds = process_log(str_list2float(abunds))
        else:
            pass
        print("%s\t%s" % (taxid, "\t".join(abunds)))

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input merge-sorted abundance file, abundance_species.xls")
    parser.add_argument("-g", "--group", metavar="FILE", type=str, default=None,
        help="Input sample grouping table,  group.list.")
    parser.add_argument("-t", "--top", metavar='INT', type=int, default=20,
        help="Species showing top X in abundance, default=20.")
    parser.add_argument("-m", "--model",  metavar="STR", type=str, default=None,
        choices=["zscore", "log", ""],
        help="zscore: zscore normalization of the data; log:, Logarithmic processing of data; choices=[zscore, log, None],  default=None.")

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
        cut_relative_abundance.py abundance_species.xls -t 20 > abundance_species_top20.xls
        #提取丰度前20的物种，计算相对丰度。
        cut_relative_abundance.py abundance_species.xls --group group.list -t 20 > group.abundance_species_top20.xls
        #根据组提取丰度前20的物种，计算相对丰度。
        cut_relative_abundance.py abundance_species.xls -t 20 -m zscore > abundance_zscore_top20.xls
        #提取丰度前20的物种，计算相对丰度并进行zscore标准化处理。
        cut_relative_abundance.py abundance_species.xls -t 20 -m log > abundance_zscore_top20.xls
        #提取丰度前20的物种，计算相对丰度并进行对数处理标准化处理。

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    cut_relative_abundance(args.input, args.group, args.top, args.model)


if __name__ == "__main__":

    main()
