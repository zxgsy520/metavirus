#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

import numpy as np
from matplotlib import pyplot as plt

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


COLOR = ['#080B74', '#303285', '#A63800', '#BF6030', '#91A200', '#ACBB2F', '#FFDB73', '#FFCE40',
     '#FFBE00', '#DFFA00', '#FF5600', '#1A1EB2', '#7375D8', '#4E51D8', '#FFA273',
     '#FF8040', '#EDFC72', '#E8FC3F', '#A67B00', '#BF9B30']


def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def analyze_phi(string):

    data = {}
    string = string.strip().split('#')
    data['uniprot_id'] = string[0]
    data['phi_id'] = ','.join(string[1].split('__'))
    data['gene_name'] = string[2]
    data['tax_id'] = string[3]
    data['pathogen_species'] = string[4].replace('_', ' ')
    data['effector'] = []
    data['chemistry_target'] = []
    data['description'] = []

    for i in string[-1].split('__'):
        if i.startswith("effector"):
            data['effector'].append(i.split('_(')[-1].strip(')').replace('_', ' '))
            continue
        if i.startswith("chemistry"):
            data['chemistry_target'].append(i.split(':_')[-1].replace('_', ' '))
            continue
        data['description'].append(i.replace('_(', '(').replace('_', ' '))

    return data


def list2str(lists):

    string = ''

    for i in lists:
        if not i:
            continue
        string += '%s;' % i
    string = string.strip(';')

    if not string:
        string = '.'

    return string


def stat_phi_m8(file, prefix):

    desc_dict = {}
    phi_dict = {}
    fp = open('%s.phi.tsv' % prefix, 'w')

    fp.write('#Gene_id\tPhi_id\tUniprot_id\tGene_name\tTax_id\tPathogen_species\tEffector\tDescription\tChemistry_target\tQstart\tQend\tEvalue\tBitscore\n')
    for line in read_tsv(file, '\t'):
        descrip = analyze_phi(line[1])

        fp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (line[0],
            descrip['phi_id'], descrip['uniprot_id'], descrip['gene_name'],
            descrip['tax_id'], descrip['pathogen_species'], list2str(descrip['effector']),
            list2str(descrip['description']), list2str(descrip['chemistry_target']), line[2], line[3], line[5], line[6]))

        for i in descrip['phi_id'].split(','):
            if i not in phi_dict:
                phi_dict[i] = [[] ,[], [], 0]
                phi_dict[i][0] += descrip['effector']
                phi_dict[i][2] += descrip['chemistry_target']
                phi_dict[i][1] += descrip['description']
            phi_dict[i][3] += 1

        for i in descrip['description']:
            if not i:
                continue
            if i not in desc_dict:
                desc_dict[i] = 0
            desc_dict[i] += 1

        for i in descrip['effector']:
            i = 'effector(%s)' % i
            if i not in desc_dict:
                desc_dict[i] = 0
            desc_dict[i] += 1

        for i in descrip['chemistry_target']:
            i = 'chemistry target(%s)' % i
            if i not in desc_dict:
                desc_dict[i] = 0
            desc_dict[i] += 1

    fp.close()

    fc = open('%s.phi_class.tsv' % prefix, 'w')
    fc.write('#Phi_id\tEffector\tDescription\tChemistry_target\tGene_number\n')
    for i in phi_dict:
        effes = list2str(set(phi_dict[i][0]))
        descs = list2str(set(phi_dict[i][1]))
        chem_target = list2str(set(phi_dict[i][2]))
        fc.write('%s\t%s\t%s\t%s\t%s\n' % (i, effes, descs, chem_target, phi_dict[i][3]))
    fc.close()

    print('#Gene_type\tNumber')
    for i in desc_dict:
        print('%s\t%s' % (i, desc_dict[i]))

    return desc_dict


def draw_pie(data_dict, prefix, sangle=120):

    data = []
    label = []
    sdata =  sum(data_dict.values()) 

    for i in data_dict:
        data.append(data_dict[i])
        label.append('{0}:{1:.2f}%'.format(i, data_dict[i]*100.0/sdata))

    #fig, ax = plt.subplots(figsize=(10, 5), subplot_kw=dict(aspect="equal"))
    plt.figure(figsize=(10, 5))
    ax =  plt.axes([0, 0.025, 0.50, 0.95])

    wedges, texts = ax.pie(data, wedgeprops= {'linewidth':0.5,'edgecolor':'white'}, startangle=sangle, colors=COLOR[0:len(data)])

    ax.legend(wedges, label,
        loc="upper left",
        frameon=False,
        fontsize="large",
        labelspacing=0.4,
        handlelength=0.7,
        handleheight=0.7,
        handletextpad=0.4,
        bbox_to_anchor=(0.90, 0, 0.70, 1))

    plt.savefig('%s.phi_classify.png' % prefix, dpi=700)
    plt.savefig("%s.phi_classify.pdf" % prefix)


def add_hlep_args(parser):

    parser.add_argument("-i", "--input", metavar='FILE', type=str, required=True,
        help="Inpu phi annotation result file.")
    parser.add_argument("-p", "--prefix", metavar='STR', type=str, default='ont',
        help="Prefix of input and output results; default=out.")

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
    stat_phi.py Subsequent processing of phi annotation results

attention:
    stat_phi.py -i *.phi.out -p name

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()
    desc_dict = stat_phi_m8(args.input, args.prefix)

    draw_pie(desc_dict, args.prefix, sangle=120)


if __name__ == "__main__":

    main()
