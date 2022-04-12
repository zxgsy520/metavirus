#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import shutil
import logging
import argparse

from ngsmetavirus.config import *
from ngsmetavirus.mngs_qc import run_mngs_qc
from ngsmetavirus.mngs_tax import run_mngs_tax
from ngsmetavirus.mngs_asm import run_mngs_asm
from ngsmetavirus.parser import add_mngs_all_args
from dagflow import DAG, Task, do_dag
from ngsmetavirus.common import check_path, mkdir


LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.0.0"


def run_mngs_all(prefix, read1, read2, ref, nohost, dtype, trim, memory, thread, job_type,
                concurrent, refresh, work_dir, out_dir, qvalue=20, atype="metagenome"):


    read1 = check_path(read1)
    read2 = check_path(read2)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    work_dict = {
        "data": "01_data",
        "rtax": "02_reads_tax",
        "asm": "03_assembly",
    }

    for k, v in work_dict.items():
        mkdir(os.path.join(work_dir, v))
        mkdir(os.path.join(out_dir, v))

    options, clean1, clean2 = run_mngs_qc(
        prefix=prefix,
        read1=read1,
        read2=read2,
        ref=ref,
        nohost=nohost, 
        dtype=dtype,
        trim=trim,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=mkdir(os.path.join(work_dir, work_dict["data"])),
        out_dir=mkdir(os.path.join(out_dir, work_dict["data"])),
        qvalue=qvalue,
        atype=atype
    )
    
    noptions = run_mngs_tax(
        prefix=prefix,
        read1=read1,
        read2=read2,
        thread=int(thread*1.5),
        job_type=job_type,
        work_dir=mkdir(os.path.join(work_dir, work_dict["rtax"])),
        out_dir=mkdir(os.path.join(out_dir, work_dict["rtax"])),
        concurrent=concurrent,
        refresh=refresh
    )
    options["software"].update(noptions["software"])    

    noptions, genome = run_mngs_asm(
        prefix=prefix,
        read1=clean1,
        read2=clean2,
        atype=atype,
        thread=thread,
        memory=memory,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=mkdir(os.path.join(work_dir, work_dict["asm"])),
        out_dir=mkdir(os.path.join(out_dir, work_dict["asm"])),
    )
    options["software"].update(noptions["software"])

    with open(os.path.join(out_dir, "ngsmeta.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    return options


def mngs_all(args):

    run_mngs_all(
        prefix=args.prefix,
        read1=args.read1,
        read2=args.read2,
        ref=args.reference,
        nohost=args.nohost,
        dtype=args.dtype,
        trim=args.trim,
        memory=args.memory,
        thread=args.thread,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        qvalue=args.qvalue,
        atype=args.atype
    )


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""


version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_mngs_all_args(parser)
    args = parser.parse_args()
    mngs_all(args)


if __name__ == "__main__":
    main()
