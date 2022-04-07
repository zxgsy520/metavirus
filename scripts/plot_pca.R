#!/usr/bin/env Rscript

library(ggforce)
library(ggplot2)
library(vegan)

options(bitmapType='cairo') #关闭服务器与界面的互动响应

add_help_args <- function(args){

    if(length(args) != 3) {
        cat("Author:XuChang\n")
        cat("Email:xc@bioyigene.com\n")
        cat("Function:Plot a pca.\n")
        cat("Example:Rscript pca.R abundance_species.xls group.list prefix\n")
        cat("Input file format:\n")
        cat("Tax Id\tsample1\tsample2\tsample3\n")
        cat("OTU1\t10\t9\t9")
        quit()
    }
}


plot_pca <- function(data, group, prefix){

    data <- read.table(data, sep="\t", row.names=1, head=TRUE, check.names=FALSE, quote="")
    groups <- read.delim(group, header=T, sep="\t", check.names=F, stringsAsFactors=F)
    colnames(groups) <- c("sample", "group")
    rownames(groups)=groups[,1]
    data <- data[,rownames(groups)]
    data_t <- t(data)
    data_t2 = merge(groups, data_t, by="row.names")

    data_nogroup.pca <- prcomp(data_t2[,-1:-3])
    data_nogroup.pca <- data.frame(data_nogroup.pca$x)

    row.names(data_nogroup.pca) <- groups[,1]
    data_nogroup.pca <- merge(data_nogroup.pca, groups, all.x=TRUE, by='row.names')

    mx <- max(data_nogroup.pca$PC1)*1.5#定义坐标轴长
    nx <- min(data_nogroup.pca$PC1)*1.5
    my <- max(data_nogroup.pca$PC2)*1.5
    ny <- min(data_nogroup.pca$PC2)*1.5

    #画图
    print(data_nogroup.pca)
    p <- ggplot(data_nogroup.pca, aes(PC1, PC2)) + geom_point(aes(colour=group, shape=group), alpha=0.5, size=3) +
                labs(x="PCA1", y="PCA2") + geom_mark_ellipse(aes(colour=group)) +
                theme_classic() + geom_vline(xintercept=0, color='gray', size=0.4) +
                geom_hline(yintercept=0, color='gray', size=0.4) +
                scale_x_continuous(limits=c(nx, mx)) + scale_y_continuous(limits=c(ny,my)) +
                theme(legend.text=element_text(size=rel(0.6)), legend.key.size=unit(0.8,'lines'), legend.title=element_blank())

    ggsave(paste(prefix, ".PCA.pdf",sep=""), p, width=120, height=100, units="mm")
    ggsave(paste(prefix, ".PCA.png",sep=""), p, width=120, height=100, units="mm")

}


args <- commandArgs(T)

add_help_args(args)

plot_pca(args[1], args[2], args[3])
