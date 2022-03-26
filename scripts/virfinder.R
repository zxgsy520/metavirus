#!/usr/bin/env Rscript
#/Work/pipeline/software/meta/virfinder/v1.1/bin/Rscript

library(VirFinder)


run_virfinder <- function(fasta, minlen=10000, score=0.9, pvalue=0.01){

    result <- VF.pred(fasta)

    for(i in 1:length(t(result[1]))){
       seqid <- t(result[1])[i]
       seqlen <- t(result[2])[i]
       seqscore <- t(result[3])[i]
       seqpvalue <- t(result[4])[i]

       if(seqlen < minlen){
           next;
       }
       if(seqscore < score){
           next;
       }
       if(seqpvalue > pvalue){
          next;
       }
       fo <- paste(seqid, seqlen, seqscore, seqpvalue, "\n", sep="\t")
       cat(fo)
    }
}


add_help_args <- function(args){

    if(length(args) != 1) {
        cat("Version: v1.0.0\n")
        cat("Author:Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:Predicted viral sequences.\n")
        cat("Example:virfinder.R genome.fasta |grep -v score >stat.tsv\n")
        quit()
    }

}

args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
run_virfinder(args[1], minlen=5000, score=0.8, pvalue=0.05)
