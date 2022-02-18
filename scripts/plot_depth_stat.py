#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

LOG = logging.getLogger(__name__)
__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep=None):
    """
    read tsv joined with sep
    :param file: file name
    :param sep: separator
    :return: list
    """
    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


# Dictionary output is in gff format
def output_dictionary(gff_dict, out_gff):

    output_gff = open(out_gff+".window_depth.tsv", 'w')

    output_gff.write("#Contig\tStart\tEnd\tDepth(X)\n")
    for line in gff_dict:
        scf_id, start, end ,depth = gff_dict[line]
        output_gff.write("{0}\t{1}\t{2}\t{3:.2f}\n".format(scf_id, start, end ,depth))
    output_gff.close()


# Draw a depth overlay map
def draw_depth(draw_list, out):

    height=len(draw_list)//2+len(draw_list)%2
    plt.style.use('ggplot')

    if len(draw_list)>1:
        plt.figure(figsize=(12, height*4.5))
    else:
        plt.figure(figsize=(12, 4.5))
    matplotlib.rcParams['ytick.direction'] = 'out'
    matplotlib.rcParams['xtick.direction'] = 'out'

    i = 0

    for scf in draw_list:
        i=i+1
        if len(draw_list)>1:
            ax = plt.subplot(height, 1, i)
        else:
            ax = plt.subplot(1,1,1)
#        matplotlib.rcParams['xtick.direction'] = 'out'
        ax.spines['top'].set_visible(False) #去掉上边框
        ax.spines['bottom'].set_visible(False) #去掉下边框
        ax.spines['left'].set_visible(False) #去掉左边框
        ax.spines['right'].set_visible(False) #去掉右边框
        plt.plot(scf['end'], scf['depth'], label=scf['scf'], linewidth ='1',color="#cb416b")
#        plt.legend(loc='upper right',frameon=False)
#        plt.legend("boxoff")
        plt.tight_layout(pad=2.5, w_pad=2.5, h_pad=1.0)
        font = {'weight': 'bold','size': 12,}

        ymax = max(scf['depth'])
        ymin = 0
        xmax = max(scf['end'])

        plt.ylim(ymin,ymax*1.2)
        plt.xlabel("Position", font)
        plt.ylabel("Depth", font)
#        plt.spines['right'].set_color('none')
#        plt.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
    plt.savefig("%s.depth.png" % out, dpi=700)
    plt.savefig("%s.depth.pdf" % out)


###
def extract_depth(window,depth_dict,out_stat):
#window=1000

    output_coverage = open(out_stat+".seq_depth.tsv",'w')
    output_coverage.write("#Contig\tSequencing depth (X)\n")

    gff_dict = {}
    draw_list = []
    k = 0

    for line  in depth_dict:
        scf_lenth = len(depth_dict[line])
        if scf_lenth >=window:
            site = range(0,scf_lenth,window)

            output_coverage.write("{0}\t{1:.2f}\n".format(line, sum(depth_dict[line])*1.00/scf_lenth))

            end = []
            depth = []
            j = 0

            for i in site:
                j = j+1
                k = k+1

                if j<len(site):
                    average_depth = sum(depth_dict[line][i:site[j]])*1.00/window
                    end.append(site[j])
                    depth.append(average_depth)

                    gff_dict[k] = (line, i+1, site[j], average_depth)
                else:
                    average_depth = sum(depth_dict[line][i:scf_lenth])*1.00/(scf_lenth-i)
                    end.append(scf_lenth)
                    depth.append(average_depth)

                    gff_dict[k] = (line, i+1, scf_lenth, average_depth)

            draw_dict = {'scf':line, 'end':end, 'depth':depth}
            draw_list.append(draw_dict)

    return gff_dict, draw_list


def read_depth(input_file):

    scf_name = []
    depth_dict = {}

    #for line in sorted(read_txt(input_file), key=lambda k: (int(k[1]), k[0])):
    for line in read_tsv(input_file, '\t'):

        scf = line[0]
        depth = int(line[2])
        if scf in scf_name:
            depth_dict[scf].append(depth)
        else:
            scf_name.append(scf)

            depth_dict[scf] = [depth]

    return depth_dict


def add_help_args(parser):

    parser.add_argument('-i', '--input', metavar='FILE', required=True,
                        help='output file of samtools depth.')
    parser.add_argument('-w', '--windows', type=int, default=500,
                        help='Set the override window size(bp).')
    parser.add_argument('-o', '--out', type=str, default='out',
                        help='Output prefix.')

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
        plot_depth_stat.py -- Draw sequencing coverage map
attention:
        plot_depth_stat.py -i <inputfile> -w <window_size>
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()
    depth_dict=read_depth(args.input)
    gff_dict, draw_list=extract_depth(args.windows, depth_dict, args.out)

    #output_dictionary(gff_dict, args.out)
    draw_depth(draw_list, args.out)


if __name__ == "__main__":
    main()
