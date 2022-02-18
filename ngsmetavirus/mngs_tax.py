#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import logging
import argparse

from ngsmetavirus.config import *
from ngsmetavirus.common import check_path, mkdir, get_version
from dagflow import DAG, Task, ParallelTask, do_dag
from ngsmetavirus.parser import add_mngs_tax_args
from ngsmetavirus import __author__, __email__, __version__

LOG = logging.getLogger(__name__)


def create_kraken_task(prefix, read1, read2, thread, job_type, work_dir, out_dir):

    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    options["software"]["kraken2"] = {
        "version": get_version(SOFTWARE_VERSION["kraken2"]),
        "option": "default"
    }

    if read1.endswith(".gz"):
        types = "--gzip-compressed"
    elif read1.endswith(".bz2"):
        types = "--bzip2-compressed"
    else:
        types = ""

    memory = thread*4
    if memory <= 50:
        memory = 50


    task = Task(
        id="kraken_%s" % prefix,
        work_dir=work_dir,
        type=job_type,
        option="-l vf=%dG -pe smp %s %s" % (memory, thread, QUEUE),
        script="""
export PATH={kraken2}:$PATH
kraken2 --db {database} --threads {thread} --use-names {types} \\
  --report {prefix}.kraken.report --output {prefix}.kraken.out \\
  --paired {read1} {read2}
python {script}/kraken2species.py {prefix}.kraken.out --taxonomy {taxonomy} >{prefix}.kraken.tax
export PATH={bracken}:$PATH
bracken -d {database} -i {prefix}.kraken.report -o {prefix}.bracken.out -w {prefix}.bracken.report -r 150 -l S
{krakentools}/kreport2mpa.py  -r {prefix}.bracken.report -o {prefix}.report
{krakentools}/kreport2krona.py -r  {prefix}.bracken.report -o {prefix}.krona_report
{krona_bin}/ktImportText {prefix}.krona_report -o {prefix}.taxonomy.html
cp {prefix}.bracken.out {prefix}.report {prefix}.taxonomy.html {out_dir}
""".format(kraken2=KRAKEN2_BIN,
            bracken=BRACKEN_BIN,
            krakentools=KrakenTools,
            krona_bin=BRONA_BIN,
            script=SCRIPTS,
            database=KRAKEN_DB,
            taxonomy=TAXONOMY,
            read1=read1,
            read2=read2,
            prefix=prefix,
            thread=thread,
            types=types,
            out_dir=out_dir
        )
    )

    return task, options, os.path.join(work_dir, "%s.taxonomy.html" % prefix)


def run_mngs_tax(prefix, read1, read2, thread, job_type, work_dir, out_dir,
                 concurrent=10, refresh=30):

    read1 = check_path(read1)
    read2 = check_path(read2)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    dag = DAG("mngs_tax")

    task, options, taxhtml = create_kraken_task(
        prefix=prefix,
        read1=read1,
        read2=read2,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(task)    

    do_dag(dag, concurrent, refresh)

    return options


def mngs_tax(args):

    options = run_mngs_tax(
        prefix=args.prefix,
        read1=args.read1,
        read2=args.read2,
        thread=args.thread,
        job_type=args.job_type,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        concurrent=args.concurrent,
        refresh=args.refresh
    )

    with open(os.path.join(args.out_dir, "mngs_tax.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    return 0


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

    parser = add_mngs_tax_args(parser)
    args = parser.parse_args()
    mngs_tax(args)


if __name__ == "__main__":
    main()
