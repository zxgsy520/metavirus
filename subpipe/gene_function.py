#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
*.nr.m6      blast results of nr database
*.nr.tsv     nr annotation file
*.KEGG.m6        blast results of KEGG database
*.KEGG.tsv       KEGG annotation file
*.KOG.out        blast results of KOG database
*.KOG.tsv        KOG annotation file
*.ipr.out        interproscan annotation file
*.ipr.tsv        interproscan annotation file

"""
import os
import sys
import argparse
import logging

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../"))
from dagflow import Task, ParallelTask, DAG, do_dag
from ngsmetavirus import __author__, __email__, __version__
from ngsmetavirus.config import *
from ngsmetavirus.common import check_paths, cd, mkdir, get_version
from seqkit.split import seq_split


LOG = logging.getLogger(__name__)
__version__ = "1.2.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


BLAST_BIN = "/Work/pipeline/software/Base/blast+/bin/"
DIAMOND_BIN = "/Work/pipeline/software/Base/diamond/v2.0.3/"
INTERPROSCAN_BIN = "/Work/pipeline/software/RNAseq/interproscan-5.51-85.0/"
RSCRIPT = "Rscript"
GS_BIN = "/Work/pipeline/software/Base/ghostscript/v9.53.3/"

SOFTWARE_VERSION = {
    "blastp": {
        "GETVER": "%s/blastp -version 2>&1| grep 'blastp'" % BLAST_BIN,
        "REGEXP": "\d+\.\d+\.\d+\+",
        "MINVER": "2.7.1+"
    },
    "rpsblast": {
        "GETVER": "%s/rpsblast -version 2>&1| grep 'rpsblast'" % BLAST_BIN,
        "REGEXP": "\d+\.\d+\.\d+\+",
        "MINVER": "2.7.1+"
    },
    "interproscan": {
        "GETVER": "%s/interproscan.sh -version 2>&1| grep 'version'" % INTERPROSCAN_BIN,
        "REGEXP": "\d+\.\d+\-\d+\.\d+",
        "MINVER": "5.30-69.0"
    },
}


NR = "/Work/database/nr/202104"
NR_TAXON = {
    "fungi": os.path.join(NR, "fungi.dmnd"),
    "animals": os.path.join(NR, "metazoa.dmnd"),
    "plants": os.path.join(NR, "viridiplantae.dmnd"),
    "meta": os.path.join(NR, "microbe.dmnd"),
}

KEGG = "/Work/database/kegg/2021/"
KEGG_KEG = os.path.join(KEGG, "ko00001.keg")
KEGG_PATH = os.path.join(KEGG, "ko2pathway.tsv")

#KEGG_KO = os.path.join(KEGG, "eukaryotes.kegg.pep2ko.txt")
KEGG_KO = os.path.join(KEGG, "microbe/microbe.kegg.pep2ko.txt")
KEGG_TAXON = {
    "fungi": os.path.join(KEGG, "fungi/fungi.dmnd"),
    "animals": os.path.join(KEGG, "animals/animals.dmnd"),
    "plants": os.path.join(KEGG, "plants/plants.dmnd"),
    "meta": os.path.join(KEGG, "microbe/microbe.dmnd")
}
COG = "/Work/database/COG/2020/"
COG_TAXON = {
    "meta": os.path.join(COG, "cog.dmnd"),
}
COG_NAME = os.path.join(COG, "cog-20.def.tab")
COG_FUNC = os.path.join(COG, "fun-20.tab")

KOG = "/Work/database/KOG/2003/"
KOG_TAXON = {
    "fungi": os.path.join(KOG, "kog.dmnd"),
    "animals": os.path.join(KOG, "kog.dmnd"),
    "plants": os.path.join(KOG, "kog.dmnd"),
}
KOG_NAME = os.path.join(KOG, "kog.def.tab")
KOG_FUNC = os.path.join(KOG, "fun-20.tab")

GO = "/Work/database/GO/20180828/"
TIGRFAMS = "/Work/database/TIGRFAMs/"
PFAM = "/Work/database/Pfam/v34.0/"
SWISSPROT_DB = "/Work/database/SwissProt/202111/swissprot.dmnd"

CAZY = "/Work/database/CAZy/v10/"
CAZY_DB = os.path.join(CAZY, "CAZy.dmnd")
CAZY_ACTIV = os.path.join(CAZY, "CAZy.activities.txt")
CAZY_SUBFAM = os.path.join(CAZY, "CAZy.subfam.txt")
PHI_DB = "/Work/database/phi/v4-12/phi"
PUL_DB = "/Work/database/dbCAN-PUL/v202010/PUL"

def create_nr_task(proteins, prefix, db, evalue, coverage, threads,
                       job_type, work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "nr"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option= "-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={diamond}:$PATH
time diamond blastp --query {{protein}} --db {db} \\
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \\
--max-target-seqs 5 --evalue 1e-05 --threads {threads} --out {{prefixs}}.nr.m6
""".format(diamond=DIAMOND_BIN,
           db=db,
           evalue=evalue,
           script=SCRIPTS,
           threads=threads),
        prefixs=prefixs,
        protein=proteins,
    )

    join_task = Task(
        id="merge_nr",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.nr.m6 > {prefix}.nr.m6
time {script}/blast_filter.py {prefix}.nr.m6 \\
--outfmt std qlen slen stitle --out qseqid sseqid qstart qend stitle evalue bitscore \\
--min_qcov {coverage} --min_scov 0 --evalue {evalue} --best > {prefix}.nr.out
{script}/refseqproc.py {prefix}.nr.out > {prefix}.nr.tsv
{script}/sure_species.py -i {prefix}.nr.tsv -n {prefix}
Rscript {script}/plot_species.R {prefix}.species.tsv {prefix}.species
cp {prefix}.nr.tsv {prefix}.nr.m6 {prefix}.species.* {out_dir}
""".format(id=id,
           prefix=prefix,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def create_interproscan_task(proteins, prefix, job_type, work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id="ipr"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp 6 -q all.q",
        script="""
export PATH={interproscan}:$PATH
time interproscan.sh -i {{protein}} -appl Pfam,TIGRFAM,SMART -iprlookup -goterms \\
  --cpu 6 -t p -f TSV -o {{prefixs}}.ipr.out
""".format(interproscan=INTERPROSCAN_BIN,
           script=SCRIPTS,
           ),
        protein=proteins,
        prefixs=prefixs

    )

    join_task = Task(
        id="merge_ipr",
        work_dir=work_dir,
        type="local",
        option="-pe smp 1 %s" % QUEUE,
        script="""
export PATH={script}:$PATH
cat {id}*/*.ipr.out > {prefix}.ipr.out
time iprproc.py {prefix}.ipr.out --tigrfams {tigrfams}/TIGRFAMS.link --prefix {prefix}  > {prefix}.other.tsv

goproc.py --go {prefix}.WEGO.txt --go_obo {go}/go.obo > {prefix}.GO.tsv
new.get_GO.classify.py {prefix}.WEGO.txt {go}/GO_level4.deal.txt > {prefix}.go.classify.xls
{rscript} {script}/get.level2_3.counts.R {prefix}.go.classify.xls {prefix}.go2.xls {prefix}.go3.xls
{rscript} {script}/GO.level.stat.R {prefix}.go2.xls {prefix}.WEGO.pdf
{gs}/gs -dNOSAFER -r600 -dBATCH -sDEVICE=pngalpha -dNOPAUSE -dEPSCrop -sOutputFile={prefix}.WEGO.png {prefix}.WEGO.pdf
cp *.pdf *.png {out_dir}
cp {prefix}.ipr.out {prefix}.*.tsv {out_dir}
""".format(gs=GS_BIN,
            id=id,
            prefix=prefix,
            tigrfams=TIGRFAMS,
            go=GO,
            rscript=RSCRIPT,
            script=SCRIPTS,
            out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def create_kegg_task(proteins, prefix, db, evalue, coverage, threads, job_type,
                     work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "KEGG"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={diamond}:$PATH
time diamond blastp --query {{protein}} --db {db} \\
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \\
--max-target-seqs 5 --evalue {evalue} --threads {threads} --out {{prefixs}}.KEGG.m6
""".format(diamond=DIAMOND_BIN,
           db=db,
           evalue=evalue,
           threads=threads),
        protein=proteins,
        prefixs=prefixs,
    )

    join_task = Task(
        id="merge_KEGG",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.KEGG.m6 > {prefix}.KEGG.m6
time {script}/blast_filter.py {prefix}.KEGG.m6 \\
--outfmt std qlen slen stitle --out qseqid sseqid qstart qend evalue bitscore \\
--min_qcov {coverage} --min_scov 0 --evalue {evalue} --best > {prefix}.KEGG.out

{script}/keggproc.py {prefix}.KEGG.out --ko {pep2ko} --pathway {pathway} --keg {keg} > {prefix}.KEGG.tsv
{script}/get_kegg_pathway.py --input {prefix}.KEGG.out -pd {pep2ko} \\
--pathway {pathway} --keg {keg} --out {prefix}.kegg_pathway.tsv
cut -f 1,3,4 {prefix}.KEGG.tsv > {prefix}.KEGG.KO
{script}/make_keg.py --keg {keg} --in {prefix}.KEGG.KO --out {prefix}.KEGG --plot

cp {prefix}.KEGG.png {prefix}.KEGG.pdf {out_dir}
cp {prefix}.kegg_pathway.tsv {prefix}.KEGG.KO {prefix}.KEGG.keg {prefix}.KEGG.tsv {prefix}.KEGG.m6 {out_dir}
""".format(id=id,
           pep2ko=KEGG_KO,
           pathway=KEGG_PATH,
           keg=KEGG_KEG,
           prefix=prefix,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def create_kog_task(proteins, prefix, db, evalue, coverage, threads, job_type,
                    work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "KOG"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={diamond}:$PATH
time diamond blastp --query {{protein}} --db {db} \\
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \\
--max-target-seqs 5 --evalue {evalue} --threads {threads} --out {{prefixs}}.KOG.m6
""".format(diamond=DIAMOND_BIN,
           db=db,
           evalue=evalue,
           threads=threads),
        protein=proteins,
        prefixs=prefixs
    )

    join_task = Task(
        id="merge_KOG",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.KOG.m6 > {prefix}.KOG.m6
time {script}/blast_filter.py {prefix}.KOG.m6 \\
--outfmt std qlen slen stitle --out qseqid sseqid qstart qend slen stitle evalue bitscore \\
--min_qcov {coverage} --min_scov 0 --evalue {evalue} --best > {prefix}.KOG.out
time {script}/cogproc.py {prefix}.KOG.out --name {kog} > {prefix}.KOG.tsv
{script}/plot_cog.py {prefix}.KOG.tsv --func {func} --name {kog} --out {prefix}.KOG
#cp {prefix}.KOG.pdf {prefix}.KOG.png {out_dir}
cp {prefix}.KOG.m6 {prefix}.KOG.tsv {out_dir}
""".format(id=id,
           prefix=prefix,
           kog=KOG_NAME,
           func=KOG_FUNC,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def create_cog_task(proteins, prefix, db, evalue, coverage, threads, job_type,
                    work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "COG"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={diamond}:$PATH
time diamond blastp --query {{protein}} --db {db} \\
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \\
--max-target-seqs 5 --evalue {evalue} --threads {threads} --out {{prefixs}}.COG.m6
""".format(diamond=DIAMOND_BIN,
           db=db,
           evalue=evalue,
           threads=threads),
        protein=proteins,
        prefixs=prefixs
    )

    join_task = Task(
        id="merge_COG",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.COG.m6 > {prefix}.COG.m6
time {script}/blast_filter.py {prefix}.COG.m6 \\
--outfmt std qlen slen stitle --out qseqid sseqid qstart qend slen stitle evalue bitscore \\
--min_qcov {coverage} --min_scov 0 --evalue {evalue} --best > {prefix}.COG.out
time {script}/cogproc.py {prefix}.COG.out --name {cog} > {prefix}.COG.tsv
{script}/plot_cog.py {prefix}.COG.tsv --func {func} --name {cog} --out {prefix}.COG
cp {prefix}.COG.pdf {prefix}.COG.png {out_dir}
cp {prefix}.COG.m6 {prefix}.COG.tsv {out_dir}
""".format(id=id,
           prefix=prefix,
           cog=COG_NAME,
           func=COG_FUNC,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def create_swissprot_task(proteins, prefix, db, evalue, coverage, threads, job_type,
                    work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "SwissProt"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={diamond}:$PATH
time diamond blastp --query {{protein}} --db {db} \\
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \\
--max-target-seqs 5 --evalue {evalue} --threads {threads} --out {{prefixs}}.SwissProt.m6
""".format(diamond=DIAMOND_BIN,
           db=db,
           evalue=evalue,
           threads=threads),
        protein=proteins,
        prefixs=prefixs
    )

    join_task = Task(
        id="merge_SwissProt",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.SwissProt.m6 > {prefix}.SwissProt.m6
time {script}/blast_filter.py {prefix}.SwissProt.m6 \\
--outfmt std qlen slen stitle --out qseqid sseqid qstart qend slen stitle evalue bitscore \\
--min_qcov {coverage} --min_scov 0 --evalue {evalue} --best > {prefix}.SwissProt.out
time {script}/swissprotproc.py {prefix}.SwissProt.out > {prefix}.SwissProt.tsv
cp {prefix}.SwissProt.m6 {prefix}.SwissProt.tsv {out_dir}
""".format(id=id,
           prefix=prefix,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def create_cazy_task(proteins, prefix, evalue, coverage, threads, job_type,
                     work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "cazy"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={diamond}:$PATH
time diamond blastp --query {{proteins}} --db {db} \\
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \\
--max-target-seqs 5 --evalue 1e-05 --threads {threads} --out {{prefixs}}.CAZy.m6
""".format(diamond=DIAMOND_BIN,
           db=CAZY_DB,
           evalue=evalue,
           coverage=coverage,
           threads=threads),
        proteins=proteins,
        prefixs=prefixs,
    )

    join_task = Task(
        id="merge_CAZy",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.CAZy.m6 > {prefix}.CAZy.m6
time {script}/blast_filter.py {prefix}.CAZy.m6 \\
  --outfmt std qlen slen stitle --out qseqid sseqid qstart qend evalue bitscore stitle \\
  --min_qcov {coverage} --min_scov 0 --evalue {evalue} --best >{prefix}.CAZy.out
{script}/cazyproc.py {prefix}.CAZy.out --activ {activ} \\
  --subfam {subfam} -o {prefix}.cazy_classify.tsv >{prefix}.cazy.tsv
{script}/stat_cazy.py {prefix}.cazy.tsv --activ {activ} >{prefix}.stat_cazy.tsv
{script}/plot_cazy.py {prefix}.cazy_classify.tsv -p {prefix}
cp {prefix}.CAZy.m6 {prefix}.CAZy.out {prefix}.cazy.tsv {out_dir}
cp {prefix}.cazy_classify.tsv {prefix}.cazy.png {prefix}.cazy.pdf {out_dir}
""".format(id=id,
           prefix=prefix,
           activ=CAZY_ACTIV,
           subfam=CAZY_SUBFAM,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task, os.path.join(work_dir, "%s.cazy.tsv" % prefix)


def create_phi_task(proteins, prefix,  evalue, coverage, threads, job_type,
                    work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "phi"
    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={blast}:$PATH
time blastp -query {{proteins}} -db {db} \\
-outfmt '6 std qlen slen stitle' \\
-max_target_seqs 5 -evalue 1e-05 -num_threads {threads} -out {{prefixs}}.phi.m6 \
""".format(blast=BLAST_BIN,
           db=PHI_DB,
           evalue=evalue,
           coverage=coverage,
           threads=threads),
        proteins=proteins,
        prefixs=prefixs,
    )

    join_task = Task(
        id="merge_phi",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.phi.m6 > {prefix}.phi.m6
time {script}/blast_filter.py {prefix}.phi.m6 \\
  --outfmt std qlen slen stitle --out qseqid sseqid qstart qend stitle evalue bitscore \\
  --min_qcov {coverage} --min_scov 0 --evalue {evalue} --best >{prefix}.phi.out
{script}/stat_phi.py -i {prefix}.phi.out \\
  -p {prefix} >{prefix}.phi_species.tsv
cp {prefix}.phi.m6 {prefix}.phi.out {prefix}.phi.tsv {prefix}.phi_class.tsv {out_dir}
cp {prefix}.phi_species.tsv {prefix}.phi_classify.png {prefix}.phi_classify.pdf {out_dir}
""".format(id=id,
           prefix=prefix,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task, os.path.join(work_dir, "%s.phi_class.tsv" % prefix)


def create_stat_function_task(protein, prefix, fun_dir, job_type,
                    work_dir="", out_dir=""):

    task = Task(
        id="stat_function",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
python {script}/merge_protein_anno.py {protein} --nr {fun_dir}/{prefix}.nr.tsv\\
  --KEGG {fun_dir}/{prefix}.KEGG.tsv --COG {fun_dir}/{prefix}.COG.tsv \\
  --TIGRFAMs {fun_dir}/{prefix}.TIGRFAMs.tsv --Pfam {fun_dir}/{prefix}.Pfam.tsv \\
  --GO {fun_dir}/{prefix}.GO.tsv --swissprot {fun_dir}/{prefix}.SwissProt.tsv > {prefix}.merge.annotate.xls
python {script}/stat_advanal_ann.py {protein} --cazy {fun_dir}/{prefix}.cazy.tsv \\
  --phi {fun_dir}/{prefix}.phi_class.tsv >{prefix}.stat_cazy_phi.tsv
cp function_summary.tsv {out_dir}/{prefix}.function_summary.tsv
cp {prefix}.merge.annotate.xls {prefix}.stat_cazy_phi.tsv {out_dir}
""".format(script=SCRIPTS,
            protein=protein,
            prefix=prefix,
            fun_dir=fun_dir,
            out_dir=out_dir
        )
    )

    return task


def create_cds_annotation_dag(protein, prefix, kingdom, evalue, coverage,
                              threads, job_type, work_dir, out_dir):

    kingdom = kingdom.lower()
    protein = check_paths(protein)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    work_dict = {
        "split": "00_data",
        "nr": "01_NR",
        "ipr": "02_ipr",
        "kegg": "03_KEGG",
        "kog": "04_KOG",
        "swissprot": "05_SwissProt",
        "cazy": "06_CAZy",
        "phi": "07_PHI",
    }
    if kingdom == "meta":
        work_dict["kog"] = "04_COG"

    for k, v in work_dict.items():
        work_dict[k] = mkdir(os.path.join(work_dir, v))

    proteins = seq_split([protein], mode="length", num=5000000, output_dir=work_dict["split"])

    _options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    _options["software"]["blastp"] = {
        "version": get_version(SOFTWARE_VERSION["blastp"]),
        "option": "-evalue %s -outfmt '6 std qlen slen stitle' -max_target_seqs 5" % evalue
    }
    _options["database"]["nr"] = {
        "version": os.path.basename(NR),
    }

    dag = DAG("CDS")
    nr_tasks, nr_join = create_nr_task(
        proteins=proteins,
        prefix=prefix,
        db=NR_TAXON[kingdom],
        evalue=evalue,
        coverage=coverage,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["nr"],
        out_dir=out_dir
    )
    dag.add_task(*nr_tasks)
    dag.add_task(nr_join)

    _options["software"]["Interproscan"] = {
        "version": get_version(SOFTWARE_VERSION["interproscan"]),
        "option": "-appl Pfam,TIGRFAM,SMART -iprlookup -goterms -t p -f TSV"
    }
    _options["database"]["Pfam"] = {
        "version": os.path.basename(PFAM),
    }
    _options["database"]["TIGRFAMs"] = {
        "version": os.path.basename(TIGRFAMS),
    }
    _options["database"]["GO"] = {
        "version": os.path.basename(GO),
    }
    ipr_tasks, ipr_join = create_interproscan_task(
        proteins=proteins,
        prefix=prefix,
        job_type=job_type,
        work_dir=work_dict["ipr"],
        out_dir=out_dir
    )
    dag.add_task(*ipr_tasks)
    dag.add_task(ipr_join)


    _options["database"]["KEGG"] = {
        "version": os.path.basename(KEGG),
    }
    kegg_tasks, kegg_join = create_kegg_task(
        proteins=proteins,
        prefix=prefix,
        db=KEGG_TAXON[kingdom],
        evalue=evalue,
        coverage=coverage,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["kegg"],
        out_dir=out_dir
    )
    dag.add_task(*kegg_tasks)
    dag.add_task(kegg_join)


    _options["software"]["rpsblast"] = {
        "version": get_version(SOFTWARE_VERSION["rpsblast"]),
        "option": "-evalue 0.01 -seg no -outfmt 5"
    }
    _options["database"]["KOG"] = {
        "version": os.path.basename(KOG),
    }
    if kingdom == "meta":
        cog_tasks, cog_join = create_cog_task(
            proteins=proteins,
            prefix=prefix,
            db=COG_TAXON[kingdom],
            evalue=evalue,
            coverage=coverage,
            threads=threads,
            job_type=job_type,
            work_dir=work_dict["kog"],
            out_dir=out_dir
        )
        dag.add_task(*cog_tasks)
        dag.add_task(cog_join)
    else:
        kog_tasks, kog_join = create_kog_task(
            proteins=proteins,
            prefix=prefix,
            db=KOG_TAXON[kingdom],
            evalue=evalue,
            coverage=coverage,
            threads=threads,
            job_type=job_type,
            work_dir=work_dict["kog"],
            out_dir=out_dir
        )
        dag.add_task(*kog_tasks)
        dag.add_task(kog_join)

    swissprot_tasks, swissprot_join = create_swissprot_task(
        proteins=proteins,
        prefix=prefix,
        db=SWISSPROT_DB,
        evalue=evalue,
        coverage=coverage,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["swissprot"],
        out_dir=out_dir
    )
    dag.add_task(*swissprot_tasks)
    dag.add_task(swissprot_join)

    cazy_tasks, cazy_join, cazy_stat = create_cazy_task(
        proteins=proteins,
        prefix=prefix,
        evalue=evalue,
        coverage=coverage,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["cazy"],
        out_dir=out_dir
    )
    dag.add_task(*cazy_tasks)
    dag.add_task(cazy_join)

    phi_tasks, phi_join, phi_stat = create_phi_task(
        proteins=proteins,
        prefix=prefix,
        evalue=evalue,
        coverage=coverage,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["phi"],
        out_dir=out_dir)
    dag.add_task(*phi_tasks)
    dag.add_task(phi_join)
    
    stat_task = create_stat_function_task(
        protein=protein,
        prefix=prefix,
        fun_dir=out_dir,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(stat_task)
    stat_task.set_upstream(nr_join)
    stat_task.set_upstream(ipr_join)
    stat_task.set_upstream(cazy_join)

    return dag, _options


def run_cds_annotation(proteins, prefix, kingdom, evalue, coverage, threads, job_type, concurrent, refresh, work_dir, out_dir):

    dag, options = create_cds_annotation_dag(proteins, prefix, kingdom, evalue, coverage, threads, job_type, work_dir, out_dir)
    do_dag(dag, concurrent_tasks=concurrent, refresh_time=refresh)


def add_cds_args(parser):
    """
    CDS annotation arguments
    :param parser:
    :return: parser
    """

    parser.add_argument("protein", metavar='FILE', type=str,
        help="Proteins in fasta file")
    parser.add_argument("-p", "--prefix", metavar="STR", required=True,
        help="Name used in pipeline")
    parser.add_argument("-k", "--kingdom", choices=["fungi", "animals", "plants", "meta"], default="meta",
        help="Kingdom of the sample[fungi, animals, plants, meta] (default: meta)")
    parser.add_argument("-t", "--threads", metavar="INT", type=int, default=1,
        help="Threads used to run blastp (default: 1)")
    parser.add_argument("-e", "--evalue", metavar="NUM", type=float, default=1e-05,
        help="Evalue cutoff of blastp for Refseq and KEGG (default: 1e-05)")
    parser.add_argument("-c", "--coverage", metavar="NUM", type=float, default=30,
        help="Coverage cutoff of blastp for Refseq and KEGG (default: 30)")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    parser.add_argument("--work_dir", metavar="DIR", default="NPGAP.work",
        help="Work directory (default: current directory)")
    parser.add_argument("--out_dir", metavar="DIR", default="NPGAP.out",
        help="Output directory (default: current directory)")

    return parser



def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
CDS annotation pipeline

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_cds_args(parser)
    args = parser.parse_args()

    run_cds_annotation(
        proteins=args.protein,
        prefix=args.prefix,
        kingdom=args.kingdom,
        evalue=args.evalue,
        coverage=args.coverage,
        threads=args.threads,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir
    )


if __name__ == "__main__":
    main()
