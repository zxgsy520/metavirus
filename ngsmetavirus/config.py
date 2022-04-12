
import os.path
from collections import OrderedDict

ROOT = "/Work/user/zhangxg/pipeline/metavirus/v1.1.1/"
BIN = os.path.join(ROOT, "bin")
TEMPLATES = os.path.join(ROOT, "template")
SCRIPTS = os.path.join(ROOT, "scripts")
DATABASE = os.path.join(ROOT, "database") 
TAXONOMY = os.path.join(DATABASE, "species.taxonomy.gz")
QUEUE = '-q all.q,s01'


OBSUTIL_BIN = ""
PYTHON_BIN = "/Work/pipeline/software/Base/miniconda3/bin/"
R_BIN = "/Work/pipeline/software/Base/miniconda/v4.10.3/bin/"

FASTP_BIN = "/Work/pipeline/software/Base/fastp/lastest/"
FASTQC_BIN = "/Work/pipeline/software/Base/fastqc/lastest/"
BLAST_BIN = "/Work/pipeline/software/Base/blast+/bin/"
HISAT2_BIN = "/Work/pipeline/software/RNAseq/hisat2-2.1.0/"
MINIMAP_BIN = "/Work/pipeline/software/Base/minimap2/"
SAMTOOLS_BIN = "/Work/pipeline/software/Base/samtools/samtools/"
SAM2NGS_BIN = "/Work/pipeline/software/Base/sam2ngs/v1.1.0/"
SEQKIT_BIN = "/Work/pipeline/software/Base/seqkit/"
PIGZ_BIN = "/Work/pipeline/software/Base/pigz/lastest/"
##taxonomy
KRAKEN2_BIN = "/Work/pipeline/software/meta/kraken2/v2.1.2/bin/"
BRACKEN_BIN = "/Work/pipeline/software/meta/bracken/v2.6.1/bin/"
KrakenTools = "/Work/pipeline/software/meta/KrakenTools/v1.2/"
BRONA_BIN = "/Work/pipeline/software/meta/krona/v2.8.1/bin/"
KRAKEN_DB = "/Work/user/liyilin/database/kraken"
#KRAKEN_DB = "/Work/database/Standard/20210517/standard"
TAXONOMY = "/Work/database/taxonomy/202111/kraken.taxonomy.gz"

##assemble
SPADES_BIN = "/Work/pipeline/software/meta/spades/v3.15.3/bin/"
MEGAHIT_BIN = "/Work/pipeline/software/meta/megahit/v1.2.9/bin/"
VIRSORTER_BIN = "/Work/pipeline/software/meta/virsorter/lastest/bin/"
CHECKV_BIN = "/Work/pipeline/software/meta/checkv/v0.8.1/bin/"
DEEPVIRFINDER_BIN = "/Work/pipeline/software/meta/deepvirfinder/v2020.11.21/bin/"


##annotation
GMHMMP_BIN = "/Work/pipeline/software/meta/MetaGeneMark/v.3.38/"
ABRICATE_BIN = "/Work/pipeline/software/meta/abricate/v1.0.1/bin/"
ABRICATE_ENV = ABRICATE_BIN

##virusreseq
SNIPPY_BIN = "/Work/pipeline/software/meta/snippy/v4.5.9/bin/"
VCFTOOLS_BIN = ""


METAGENEMARK_BIN = "/export/personal/software/software/MetaGeneMark/v.3.38/"
CDHIT_BIN = "/export/personal/software/software/cdhit/v4.8.1/"
BEDTOOLS_BIN = "/export/personal/software/software/bedtools/v2.92.2/"
DIRSEQ_BIN = "/export/personal/software/software/dirseq/v0.4.1/bin/"



GENETIC_CODES = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, ]

#TAXONOMY="/export/personal/software/database/kraken/20200420/microbe/taxonomy/names.dmp"

PORECHOP_BIN = "/export/personal/software/software/porechop/v0.2.4/bin/"
#NANOFILT_BIN
NANOFILT_BIN = "/export/personal/software/software/nanofilt/v2.6.0/bin/"
#FILTLONG_BNI
FILTLONG_BNI = "/export/personal/software/software/Filtlong/bin"
SAMBLASTER_BIN = "/export/personal/software/software/samblaster/v0.1.25"

#METATAX
#CURRENT_BIN="/export2/master2/wangn/software/03_TaxAnnotaion/current"
CENTRIFUGR_BIN= "/export/personal/software/software/centrifuge/v1.0.4/bin"
CENTRIFUGR_DB= "/export/personal/software/database/centrifuge/20200318/prokaryotic/hpv"


#METAASM
FLYE_BIN = "/export/personal/software/software/flye/v2.8/bin/"
QUICKMERGE_BIN = "/export/personal/software/software/quickmerge/v0.3/bin/"
PBMM2_BIN = "/export/personal/software/software/miniconda/v4.8.2/envs/pbmm2/bin/"
GCPP_BIN = "/export/personal/software/software/miniconda/v4.8.2/envs/pbgcpp/bin/"
QUAST_BIN = "/export/personal/software/software/quast/v5.2.0/bin/"
MEDAKA_BIN = "/export/personal/software/software/miniconda/v4.8.2/envs/medaka/bin/"

#MINIMAP
BWA_BIN = "/Work/pipeline/software/Base/bwa/lastest/"

#METAANT
PRODIGAL_BIN = "/export/personal/software/software/prodigal"
PROKKA_BIN = "/export/personal/software/software/prokka/v1.14.6/bin/"
PROKKA_SOFT = "/export/personal/software/software/prokka/v1.14.6/binaries/linux"
DIAMOND_BIN = "/export/personal/software/software/diamond/v0.9.34/"

#METABIN
METABAT_BIN = "/export/personal/software/software/metabat/v2.12.1/"
PRODIGAL_BIN = "/export/personal/software/software/prodigal/"
HMMER_BIN = "/export/personal/software/software/hmmer/3.2.1/bin/"
PPLACER_BIN = "/export/personal/software/software/pplacer/v1.1/"
CHECKM_BIN = "/export/personal/software/software/checkm/v1.1.3/bin/"


# SOFTWARE VERSION

SOFTWARE_VERSION = {
    "fastp":{
        "GETVER": "%s/fastp --version 2>&1 |grep 'fastp'" % FASTP_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.20.0",
    },
    "fastqc":{
        "GETVER": "%s/fastqc --version 2>&1 |grep 'FastQC'" % FASTQC_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.11.8",
    },
    "blastn":{
        "GETVER": "%s/blastn -version 2>&1| grep 'blastn'" % BLAST_BIN,
        "REGEXP": "\d+\.\d+\.\d+\+",
        "MINVER": "2.7.1+"
    },
    "hisat2":{
        "GETVER": "%s/hisat2 --version 2>&1 |grep 'hisat'" % HISAT2_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2.1.0",
    },
    "samtools": {
        "GETVER": "%s/samtools --version" % SAMTOOLS_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "1.4",
    },
    "minimap2": {
        "GETVER": "%s/minimap2 --version" % MINIMAP_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "2.17",
    },
    "kraken2": {
        "GETVER": "%s/kraken2 --version|grep 'version'" % KRAKEN2_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2.0.8",
    },
    "spades": {
        "GETVER": "%s/spades.py --version 2>&1|grep SPAdes" % SPADES_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "v3.15.3",
    },
    "megahit" : {
        "GETVER": "%s/megahit --version 2>&1|grep MEGAHIT" % MEGAHIT_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "1.2.9",
    },
    "virsorter":{
        "GETVER": "%s/virsorter -h 2>&1| grep 'VirSorter'" % VIRSORTER_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2.2.2"
    },
    "deepvirfinder":{
        "GETVER": "ls %s/dvf.py" % DEEPVIRFINDER_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2020.11.21"
    },
    "checkv":{
        "GETVER": "ls %s/checkv" % CHECKV_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.8.1"
    },
    "MetaGeneMark" : {
        "GETVER": "%s/gmhmmp 2>&1|grep 'version'" % METAGENEMARK_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "3.38",
    },
    "cd-hit" : {
        "GETVER": "%s/cd-hit --version|grep 'CD-HIT version'" % CDHIT_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "4.8.1",
    },
    "kraken2": {
        "GETVER": "%s/kraken2 --version|grep 'version'" % KRAKEN2_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2.0.8",
    },
    "centrifuge": {
        "GETVER": "%s/centrifuge --version|grep 'centrifuge-class version'" % CENTRIFUGR_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "1.0.4",
    },
    "bracken": {
        "GETVER": "ls %s/bracken" % BRACKEN_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "2.2",
    },
    "flye": {
        "GETVER": "%s/flye --version 2>&1" % FLYE_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "2.6",
    },
    "quickmerge" : {
        "GETVER": "%s/quickmerge --version 2>&1|grep quickmerge" % QUICKMERGE_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "0.3",
    },
    "bwa": {
        "GETVER": "%s/bwa 2>&1|grep -i '^Version:'" % BWA_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.7.17"
    },
    "pbmm2": {
        "GETVER": "%s/pbmm2 --version 2>&1|cut -d '(' -f1" % PBMM2_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "1.0.0",
    },
    "gcpp":{
       "GETVER": "%s/gcpp --version 2>&1|cut -d '-' -f1" % GCPP_BIN,
       "REGEXP": "\d+\.\d+\.\d+",
       "MINVER": "1.9.0"
    },
    "medaka": {
        "GETVER": "%s/medaka --version" % MEDAKA_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.10.0"
    },
    "metabat": {
        "GETVER": "%s/metabat -h 2>&1|grep 'version'" % METABAT_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2.6.3",
    },
    "prodigal": {
        "GETVER": "%s/prodigal -v 2>&1|grep 'Prodigal'" % PRODIGAL_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2.6.3",
    },
    "prokka": { 
        "GETVER": "%s/prokka -v 2>&1|grep 'prokka'" % PROKKA_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "1.14.5",
    },  
    "hmmer": {
        "GETVER": "%s/hmmsearch -h|grep 'HMMER'" % HMMER_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "3.2.1",
    },
    "pplacer": {
        "GETVER": "%s/pplacer --version" % PPLACER_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "1.1",
    },
    "checkm": {
        "GETVER": "%s/checkm -h |grep ': CheckM'" % CHECKM_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "1.0.18",
    },
    "prodigal": {
        "GETVER": "%s/prodigal -v 2>&1| grep 'Prodigal'" % PRODIGAL_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2.6.3"
    },
    "blastp": {
        "GETVER": "%s/blastp -version 2>&1| grep 'blastp'" % BLAST_BIN,
        "REGEXP": "\d+\.\d+\.\d+\+",
        "MINVER": "2.7.1+"
    },
    "abricate": {
        "GETVER": "%s/abricate --version 2>&1| grep 'abricate'" % ABRICATE_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.9.8"
    },
}

# DATABASE CONFIGURE
NT_DATABASE = "/Work/database/nt/20210922/"
NT = os.path.join(NT_DATABASE, "nt")
NT_TAXON = {
    "meta": NT,
}


CAZY = "/export/personal/software/database/CAZy/v9/"
CAZY_SUBFAM = os.path.join(CAZY, "CAZy.subfam.ec.txt")
CAZY_TAXON = {
    "Bacteria": os.path.join(CAZY, "CAZy"),
    "Viruses": os.path.join(CAZY, "CAZy"),
    "plasmid": os.path.join(CAZY, "CAZy"),
    "Archaea": os.path.join(CAZY, "CAZy"),
    "Mitochondria": "",
    "meta": os.path.join(CAZY, "CAZy")
}

KEGG = "/export/personal/software/database/KEGG/20180701/"
KEGG_KEG = os.path.join(KEGG, "ko00001.keg")
KEGG_TAXON = {
    "Bacteria": os.path.join(KEGG, "prokaryotes.kegg.dmnd"),
    "Viruses": "",
    "plasmid": os.path.join(KEGG, "prokaryotes.kegg.dmnd"),
    "Archaea": os.path.join(KEGG, "prokaryotes.kegg.dmnd"),
    "Mitochondria": "",
    "meta": os.path.join(KEGG, "prokaryotes.kegg.dmnd")
}

NOG = "/export/personal/software/database/NOG/20190302/"
NOG_AN = os.path.join(NOG, "nog.xls")
NOG_FUNC = os.path.join(NOG, "fun.tab")
NOG_NAME = os.path.join(NOG, "nognames.tab")
NOG_TAXON = {
    "Bacteria": os.path.join(NOG, "nog.dmnd"),
    "Viruses": "",
    "plasmid": os.path.join(NOG, "nog.dmnd"),
    "Archaea": os.path.join(NOG, "nog.dmnd"),
    "Mitochondria": "",
    "meta": os.path.join(NOG, "nog.dmnd")
}

ABRICATE_DB = "/export/personal/software/software/abricate/v1.0.1/db/"
ABRICATE_TAXON = {
    "card": os.path.join(ABRICATE_DB, "card/sequences"),
    "ncbi": os.path.join(ABRICATE_DB, "ncbi/sequences"),
    "vfdb": os.path.join(ABRICATE_DB, "vfdb/sequences"),
    "ecoh": os.path.join(ABRICATE_DB, "ecoh/sequences"),
    "plasmidfinder": os.path.join(ABRICATE_DB, "plasmidfinder/sequences"),
    "ecoli_vf": os.path.join(ABRICATE_DB, "ecoli_vf/sequences"),
    "resfinder": os.path.join(ABRICATE_DB, "resfinder/sequences"),
    "argannot": os.path.join(ABRICATE_DB, "argannot/sequences"),
    "abricate": os.path.join(ABRICATE_DB, "abricate/sequences")
}
