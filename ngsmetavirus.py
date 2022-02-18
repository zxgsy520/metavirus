#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from ngsmetavirus.parser import *
#from ngsmeta.obs_data import obs_data
from ngsmetavirus.mngs_qc import mngs_qc
from ngsmetavirus.mngs_asm import mngs_asm
#from ngsmeta.annotation import annotation
from ngsmetavirus.mngs_all import mngs_all
from ngsmetavirus.mngs_multi import mngs_multi, add_mngs_multi_args
from ngsmetavirus import __version__, __email__, __author__

def add_ngsmeta_parser(parser):

    subparsers = parser.add_subparsers(
        title='command',
        dest='commands')
    subparsers.required = True

    mngs_multi_parser = subparsers.add_parser('multi', help="Multiple samples in parallel")
    mngs_multi_parser = add_mngs_multi_args(mngs_multi_parser)
    mngs_multi_parser.set_defaults(func=mngs_multi)

    mngs_all_parser = subparsers.add_parser("all", help="Metagenomics and macrotranscriptome analysis process")
    mngs_all_parser = add_mngs_all_args(mngs_all_parser)
    mngs_all_parser.set_defaults(func=mngs_all)

    mngs_qc_parser = subparsers.add_parser('mngs_qc', help="Perform data quality control and remove the host")
    mngs_qc_parser = add_mngs_qc_args(mngs_qc_parser)
    mngs_qc_parser.set_defaults(func=mngs_qc)

    mngs_asm_parser = subparsers.add_parser('mngs_asm', help="Assembly of metagenomics or Macrotranscriptome")
    mngs_asm_parser = add_mngs_asm_args(mngs_asm_parser)
    mngs_asm_parser.set_defaults(func=mngs_asm)
    
#    annotation_parser = subparsers.add_parser('annotation', help="Annotation of genome and transgroup")
#    annotation_parser = add_annotation_args(annotation_parser)
#    annotation_parser.set_defaults(func=annotation)
    
    return parser


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Next-generation sequencing metagenomic analysis process

version: %s
contact:  %s <%s>\
        """ % (__version__, " ".join(__author__), __email__))

    parser = add_ngsmeta_parser(parser)
    args = parser.parse_args()

    args.func(args)

    return parser.parse_args()


if __name__ == "__main__":
    main()
