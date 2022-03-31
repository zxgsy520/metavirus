#!/usr/bin/env Rscript

library(plyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
options(bitmapType='cairo')

add_help_args <- function(args){

    if(length(args) != 3) {
        cat("Version: v1.0.1\n")
        cat("Author:Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:Plot the abundance of each sample.\n")
        cat("Example:plot_abundance_bar.R abundance_species_top20.tsv level prefix\n")
        cat("Input file format:\n")
        cat("\tsample1\tsample2\tsample3\n")
        cat("OTU1\t9.276046\t9.339209\t9.063164")
        quit()
    }

}


COLORS <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
           "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "gray","#FFB6C1",
           "#3A5FCD","#7CCD7C","#8A2BE2","#CDAD00","#EEEE00",
           "#8B6914","#FF00FF","#8F8F8F","#698B22","#87CEFF","#C71585","#EE9A49",
           "#00FFFF","#FFB6C1","#DC143C","#DB7093","#FF1493","#DA70D6","#8B008B",
           "#BA55D3","#9400D3","#4B0082","#9370DB","#6A5ACD","#483D8B","#0000CD",
           "#00008B","#4169E1","#B0C4DE","#FF0000")


plot_abundance_bar <- function(file, level, prefix){

    data <- read.table(file, header=T, sep="\t", check.names=FALSE)
    data <- dplyr::filter(data, !grepl('unclassified|other', data[,1])) #过滤掉没分类的和分类到其他的
    rownames(data) <- data[,1]
    colnames(data)[1] <- c("taxid")

    if(ncol(data) <= 6){
        width <- 16
    }else if(ncol(data) <= 11){
        width <- 28
    }else{
        width <- 35
    }
    if(level=="species"){
        width <- width + 2
    }
    
    data_melt <- melt(data)
    colnames(data_melt) <- c("taxid", "sample", "value")
    data_melt$taxid <- factor(data_melt$taxid, levels = rev(data$taxid))

    ggplot(data_melt, aes(x=sample, y=value, fill=taxid)) + geom_bar(stat="identity", width=0.8) +
        xlab(NULL) + ylab("Relative abundance ") + theme_bw()+ scale_y_continuous(expand=c(0, 0)) +
        theme(axis.title.y=element_text(size=14), axis.text.y=element_text(size=12),
            axis.text.x=element_text(angle=45, size=12, hjust=0.5, vjust=0.5), 
            legend.title=element_text(size=12), legend.text=element_text(size=11),
            legend.key.width=unit(0.3, "cm"), legend.key.height=unit(0.2, "cm"),
            panel.border=element_blank(), panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(), axis.line=element_line(colour="black")
        ) + guides(fill=guide_legend(ncol=1, title=NULL)) + #geom_text(aes(y=lable_y, label=value), colour="white", size=2, hjust=0.6) +
        scale_fill_manual(values =COLORS)

    ggsave(paste(prefix, paste0(level, ".pdf"), sep='_'), units="cm", width=width, height=15)

    ggsave(paste(prefix, paste0(level, ".png"), sep='_'), units="cm", width=width, height=15)

}

args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
plot_abundance_bar(args[1], args[2], args[3])
