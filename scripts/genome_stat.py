#!/usr/bin/env python

import re
import os
import sys
import logging
import argparse
from collections import defaultdict

from Bio import SeqIO


def n50(lengths):

    sum_length = sum(lengths)
    accu_length = 0

    for i in sorted(lengths, reverse=True):
        accu_length += i

        if accu_length >= sum_length*0.5:
            return i


def stat_genome(genome, out):

    lengths = []
    base_dict = {}

    for record in SeqIO.parse(genome, "fasta"):
        seq = record.seq.upper()
        a_number = seq.count('A')
        t_number = seq.count('T')
        g_number = seq.count('G')
        c_number = seq.count('C')
        n_number = seq.count('N')
    
        lengths.append(len(seq))
        base_dict[record.id] = [a_number, t_number, g_number, c_number, n_number]

    with open(out, "w") as fh:
        fh.write("""\
#Total bases\tContig number\tContig N50\tLongest Contig\tShortest Contig
{0:,}\t{1:,}\t{2:,}\t{3:,}\t{4:,}
""".format(sum(lengths), len(lengths), n50(lengths), max(lengths), min(lengths) ))
    
    return base_dict


def stat_gc(gc_dict, out):

    sum_a = sum_t = sum_g = sum_c = sum_n = 0

    for a, t, g, c, n in gc_dict.values():
        sum_a += a
        sum_t += t
        sum_g += g
        sum_c += c
        sum_n += n

    sum_all = sum_a + sum_t + sum_g + sum_c + sum_n

    with open(out, "w") as fh:
        fh.write("""\
#Base\tNumber\t% of total
A\t{0:,}\t{1:.2f}
T\t{2:,}\t{3:.2f}
G\t{4:,}\t{5:.2f}
C\t{6:,}\t{7:.2f}
N\t{8:,}\t{9:.2f}
GC\t{10:,}\t{11:.2f}
Total\t{12:,}\t100.00
""".format(
            sum_a, sum_a*100.0/sum_all,
            sum_t, sum_t*100.0/sum_all,
            sum_g, sum_g*100.0/sum_all,
            sum_c, sum_c*100.0/sum_all,
            sum_n, sum_n*100.0/sum_all,
            sum_g+sum_c, (sum_g+sum_c)*100.0/sum_all,
            sum_all
    ))
    

def genome_gc_args(parser):
    parser.add_argument('--fasta', metavar='FILE', required=True,
        help='enter the genome file (fasta).')
    parser.add_argument('--contig', metavar='STR', default='genome_stat.tsv',
        help='out the genome stat.')
    parser.add_argument('--base', metavar='STR', default='genoe_base.tsv',
        help='out the genome base.')
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
     stat_genome_gc.py -- Statistical genome length and base content

attention:
    stat_genome_gc.py --fasta genome.fasta
''')
    args = genome_gc_args(parser).parse_args()
    base_dict = stat_genome(args.fasta,args.contig)
    stat_gc(base_dict,args.base)


if __name__ == "__main__":
    main()
