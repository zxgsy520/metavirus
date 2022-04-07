#!/usr/bin/env Rscript

library(ggforce)
library(ggplot2)
library(vegan)

options(bitmapType='cairo') #关闭服务器与界面的互动响应

add_help_args <- function(args){

    if(length(args) != 3) {
        cat("Author:XuChang\n")
        cat("Email:xc@bioyigene.com\n")
        cat("Function:Plot a pcoa.\n")
        cat("Example:Rscript plot_pcoa.R Bacteria.abundance_species.xls group.list prefix\n")
        cat("Input file format:\n")
        cat("Tax Id\tsample1\tsample2\tsample3\n")
        cat("OTU1\t10\t9\t9")
        quit()
    }

}


plot_pcoa <- function(data, group, prefix){

    #数据准备
    raw_data <- read.table(data, sep= "\t", row.names=1, head=TRUE, check.names=FALSE, quote="")
    groups <- read.delim(group, header=T, sep="\t", check.names=F, stringsAsFactors=F)
    colnames(groups) <- c("sample", "group")
    rownames(groups)=groups[,1]


    data <- apply(raw_data, 2, as.numeric)
    row.names(data) <- row.names(raw_data)
    distance <- vegdist(t(data), method="bray")
    pcoa <- cmdscale(distance, eig=TRUE)

    poing <- data.frame(pcoa$points)#查看主要排序轴的特征值和各样本在各排序轴中的坐标值
    write.csv(poing, paste(prefix, ".pcoa.csv", sep ="")) #可将样本坐标转化为数据框后导出，例如导出为csv格式

    pcoa_eig <- (pcoa$eig)[1:2]/sum(pcoa$eig)#坐标轴解释量（前两轴)
    sample_site<- data.frame(pcoa$point)[1:2]#提取样本点坐标（前两轴)
    names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
    sample_site <- merge(sample_site, groups, by = 'row.names', all.x = TRUE)
    rownames(sample_site)=sample_site[,1]

    write.csv(sample_site, paste(prefix, ".compare_site.csv", sep=""), quote=F)

    #绘图
    mx <- max(sample_site$PCoA1)*1.5#定义坐标轴长度
    nx <- min(sample_site$PCoA1)*1.5
    my <- max(sample_site$PCoA2)*1.5
    ny <- min(sample_site$PCoA2)*1.5
    p <- ggplot(sample_site, aes(PCoA1,PCoA2)) + geom_point(aes(colour=group, shape=group), alpha=0.5, size=3) +
                xlab(paste("PCoA1=", round(100 * pcoa_eig[1], 2), "%", sep="")) + ylab(paste("PCoA2=", round(100 * pcoa_eig[2], 2), "%", sep="")) +
                geom_mark_ellipse(aes(colour=group)) + theme_classic() + geom_vline(xintercept=0, color="gray", size=0.4) +
                geom_hline(yintercept=0, color="gray", size=0.4) + geom_mark_ellipse(aes(colour=group)) +
                scale_x_continuous(limits=c(nx,mx)) + scale_y_continuous(limits=c(ny,my)) +
                theme(legend.text=element_text(size = rel(0.6)), legend.key.size=unit(0.8, "lines"), legend.title=element_blank())

    ggsave(paste(prefix, ".PCoA.pdf", sep=""), p, width=120, height=100, units="mm")
    ggsave(paste(prefix, ".PCoA.png", sep=""), p, width=120, height=100, units="mm")
}


args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)

plot_pcoa(args[1], args[2], args[3])
