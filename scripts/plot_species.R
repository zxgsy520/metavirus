#!/usr/bin/env Rscript

library(ggplot2)
options(bitmapType='cairo')


add_help_args <- function(args){

    if(length(args) != 2) {
        cat("Version: v1.0.0\n")
        cat("Author:Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:Plot a species distribution pie chart.\n")
        cat("Example:plot_species.R species.tsv species_distribution\n")
        cat("Input file format:\n")
        cat("\tSpecies\tGene number\n")
        cat("Galerina marginata\t2726")
        quit()
    }
}

plot_species <- function(file, prefix){

    data <- read.table(file, header=T, sep="\t")
    colnames(data) <- c("species_name", "gene_number")
    gnumber <- as.vector(data$gene_number)
    gspecies <- as.vector(data$species_name)
    gtotal <- sum(gnumber)
    percent <- paste(round(gnumber*1.0/gtotal*100, 2), "%", sep="")
    mylabel <- paste(gspecies, "(", percent, ")", sep="")

    p <- ggplot(data, aes(x="", y=gnumber, fill=factor(gspecies, levels=rev(gspecies)))) +
        geom_bar(stat='identity', width=1) + coord_polar(theta = "y") +
        labs(x = "", y = "", title = " ") +
        theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
        theme(plot.title=element_text(hjust=0.5)) +
        theme(panel.background=element_blank()) +
        scale_fill_discrete(name="Species", breaks=gspecies, labels=mylabel)

        ggsave(file=paste0(prefix, ".pdf"), plot=p)
        ggsave(file=paste0(prefix, ".png"), plot=p, height=5)
}


args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
plot_species(args[1], args[2])
