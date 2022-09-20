#!/usr/bin/env Rscript

library("ggplot2")
library("ggpubr")
#library("argparse")

options(bitmapType='cairo') #关闭服务器与界面的互动响应


get_compare_group <- function(groups){

    groups <- unique(groups)

    compare <- list()
    n <- 0
    for (i in 1:(length(groups)-1)) {
        for (j in (i+1):length(groups)){
            n <- n+1
            compare[[n]] <- c(groups[i], groups[j])
        }
    }

    return(compare)

}


run_plot_test_box <- function(data, group, prefix){
    
     data <- read.table(data, header=T, check.names=FALSE)
     group <- read.delim(group, sep="\t", stringsAsFactors=FALSE, header=TRUE)
     colnames(group) <- c("sample", "group")
     rownames(group) <- group$sample #对组进行排序对应

     sample_id <- row.names(data)
     group <- group[sample_id, ]
     group$Chao1 <- data$Chao1

     grouplen <- length(unique(group$group))
     if(grouplen <= 2){
         width <- 10
     }else if(grouplen <= 4){
         width <- 12
     }else {
         width <- 12 + (grouplen-4)*0.6
     }

     compare <- get_compare_group(group$group)
     p <- ggboxplot(group, x="group", y="Chao1", color="group", palette="jco", add="jitter", legend="none") +
         stat_compare_means(comparisons=compare, method="t.test", label="p.signif")
     ggsave(paste(prefix, "alpha_test_box.png", sep="."), units="cm", width=width, height=10)
     ggsave(paste(prefix, "alpha_test_box.pdf", sep="."), units="cm", width=width, height=10)
}


add_help_args <- function(args){

    #parser <- ArgumentParser()
    #parser$add_argument("-i", "--input", 
    #    help = "Input alpha diversity statistics(stat_alpha_diversity.tsv).")
    #parser$add_argument("-g", "--group",
    #    help = "Input group file(group.list).")
    #parser$add_argument("-p", "--prefix", default="output",
    #    help = "Output file prefix(prefix).")

    if(length(args) != 3) {
        cat("Version: v1.0.0\n")
        cat("Author:Boya Xu, Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:anosim analysis.\n")
        cat("Example:plot_test_box.R stat_alpha_diversity.tsv group.list prefix\n")
        cat("Input file format:\n")
        cat("\tRichness\tShannon\tSimpson\tChao1\tACE\tCoverage\n")
        cat("A1\t1471\t2.2472\t0.6335\t1471\t1471\t100.0000\n")
        quit()
    }
}

args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
run_plot_test_box(args[1], args[2], args[3])
