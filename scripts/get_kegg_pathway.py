#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def list_to_str(lists, start, end):
    
    out_list = []
    n = 0

    for i in lists:
        n += 1
        if n<=start:
            out_list.append(i)
        else:
            break
    out_list.append(" ".join(lists[start:end]))

    return out_list


def read_kegg_keg(input_files):
    
    keg_dict = {}

    for line in open(input_files, 'r'):
        line = line.strip()

        if not line or line.startswith('#'):
            continue
        
        if line.startswith('C'):
            line = line.split()
            line = list_to_str(line, 2, -1)
        else:
            continue

        keg_dict['ko'+line[1]] = line[2]
    
    return keg_dict


def read_pathway(input_filse):
    
    for line in open(input_filse, 'r'):
        line = line.strip()
        
        if not line or line.startswith('#'):
            continue

        line = line.replace(':', '\t').split('\t')
        yield line


def generate_smon(pep2ko, ko2pathway):
    
    ko_dict = {}
    pathway_dict = {}

    for line in read_pathway(pep2ko):
        if line[1] not in ko_dict:
            ko_dict[line[1]] = []
            ko_dict[line[1]].append(line[3])
        else:
            ko_dict[line[1]].append(line[3])

    for line in read_pathway(ko2pathway):
       if line[3].startswith('map'):
           continue
       
       if line[1] not in pathway_dict:
           pathway_dict[line[1]] = []
           pathway_dict[line[1]].append(line[3])
       else:
           pathway_dict[line[1]].append(line[3])

    return ko_dict, pathway_dict


def read_m6(input_filse):
    
    for line in open(input_filse, 'r'):
        line = line.strip()

        if not line or line.startswith('#'):
            continue

        line = line.replace(':', '\t').split('\t')
        yield line[0:3]


def construct_map(input_m6, pep2ko, ko2pathway):
    
    map_dict = {}
    ko_dict, pathway_dict = generate_smon(pep2ko, ko2pathway)

    for line in read_m6(input_m6):
        if line[2] not in ko_dict:
            continue        
        for i in ko_dict[line[2]]:
            if i not in pathway_dict:
                continue
            for j in pathway_dict[i]:
                if j in map_dict:
                    map_dict[j][0].add(i)
                    map_dict[j][1].add(line[0])                    
                else:
                    map_dict[j] = [set(),set()]
                    map_dict[j][0].add(i)
                    map_dict[j][1].add(line[0])
    ko_dict = {}
    pathway_dict = {}
    return map_dict


def out_pathway(map_dict, keg_files, out_files):

     keg_dict = read_kegg_keg(keg_files)
     output = open(out_files, 'w')
     output.write('#Pathway\tDefinition\tInternet site\tGene\tNumber\n')

     for name in map_dict:
         n = 0
         if name in keg_dict:
             cmd = '%s\t%s\thttps://www.kegg.jp/kegg-bin/show_pathway?%s' % (name,keg_dict[name],name)
         else:
             cmd = '%s\t%s\thttps://www.kegg.jp/kegg-bin/show_pathway?%s' % (name,'',name)
         for i in map_dict[name][0]:
             cmd += '+'+i
         cmd += '\t'
         for i in map_dict[name][1]:
             cmd += i+','
             n += 1
         output.write(cmd+'\t'+str(n)+'\n')


def get_pathway(parser):

    parser.add_argument('-i', '--input', metavar='STR', type=str, required=True,
        help='Enter the protein and kegg library for comparison and filtered result files.')
    parser.add_argument('-pd', '--pepid', metavar='STR', type=str,
        default='/export/pipeline/SGenome/NPGAP/databases/KEGG/87.0-r20180701/prokaryotes.kegg.pep2ko',
        help='Enter the protein id and ko correspondence file.')
    parser.add_argument('-p', '--pathway',  metavar='STR', type=str,
        default='/export/pipeline/SGenome/NPGAP/databases/KEGG/87.0-r20180701/ko2pathway.tsv',
        help='Enter the correspondence file between ko and map.')
    parser.add_argument('-k', '--keg', metavar='STR', type=str,
        default='/export/pipeline/SGenome/NPGAP/databases/KEGG/87.0-r20180701/ko00001.keg',
        help='Enter the comment file for ko and map.')
    parser.add_argument('-o', '--out', metavar='STR', type=str, default='kegg.http.tsv',
        help='The name of the output file')
    
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
     get_kegg_pathway.py Search for all possible metabolic pathways of this species.

attention:
     get_kegg_pathway.py -i ZFPS3.KEGG.out
     get_kegg_pathway.py -i ZFPS3.KEGG.out -pd prokaryotes.kegg.pep2ko -p ko2pathway.tsv -k ko00001.keg -o kegg.http.tsv
''')
    args = get_pathway(parser).parse_args()
    map_dict = construct_map(args.input, args.pepid, args.pathway)
    out_pathway(map_dict, args.keg, args.out)


if __name__ == "__main__":

    main()
