#!/usr/bin/env Rscript

library(vegan)
library(ggplot2)
library(ggforce)

options(bitmapType='cairo') #关闭服务器与界面的互动响应
#导入数据

nmds_analyze <- function(abundance, group){

    raw_data <- read.table(abundance, sep ="\t", row.names=1, head=TRUE, check.names=FALSE, quote="")
    groups <- read.delim(group, header=T, sep="\t", row.names=1, check.names=F, stringsAsFactors=F)
    colnames(groups) <- c("group")

    data <- apply(raw_data, 2, as.numeric)
    row.names(data) <- row.names(raw_data)
    groups$value <- groups$group
    index <- rownames(groups) %in% colnames(data)
    groups <- groups[index,]
    data <- data[,rownames(groups)]

    #NMDS分析
    nmds1 <- metaMDS(t(data))
    nmds1.stress <- nmds1$stress
    nmds1.point <- data.frame(nmds1$point)
    nmds1.species <- data.frame(nmds1$species)
    sample_site <- nmds1.point[1:2]
    sample_site$names <- rownames(sample_site)
    sample_site <- merge(sample_site, groups, by="row.names", all.x=TRUE)
    rownames(sample_site) <- sample_site[,1]
    sample_site <- sample_site[,-1]
    colnames(sample_site)[1:2] <- c('NMDS1', 'NMDS2')
    result <- list(sample=sample_site, nmds=nmds1)
    return(result)

}


plot_nmds <- function(abundance, group, prefix){

    #NMDS图绘制
    a <- nmds_analyze(abundance=abundance, group=group)
    sample_site <- a$sample
    nmds1 <- a$nmds
    #设置坐标轴范围
    x_n <- min(sample_site$NMDS1)*2
    x_m <- max(sample_site$NMDS1)*2
    y_n <- min(sample_site$NMDS2)*2
    y_m <- max(sample_site$NMDS2)*2

    p <- ggplot(sample_site, aes(NMDS1, NMDS2)) + geom_point(aes(color=group, shape=group), size=3, alpha=0.8) + #修改点的透明度、大小
                theme(panel.grid=element_blank(), panel.background=element_rect(color="black", fill="transparent")) + #去掉背景
                theme(legend.key=element_rect(fill="transparent"), legend.title=element_blank(), legend.key.size=unit(0.9, "lines"),
                      legend.text=element_text(size= rel(0.7)), legend.position=c(1.1,1), plot.margin=unit(c(1.5, 7, 0.5, 0.5), "lines")) + #修改图例
                labs(x="NMDS1", y="NMDS1", title=paste("NMDS Stress =", round(nmds1$stress, 4))) +
                scale_x_continuous(limits=c(x_n,x_m)) + scale_y_continuous(limits=c(y_n,y_m)) + #设置坐标轴范围
                theme(plot.title=element_text(hjust=0.5)) +  #标题居中
                theme(panel.background=element_blank(), axis.line=element_line(color="black"), plot.margin=unit(c(1.5, 7, 0.5, 0.5), "lines")) +
                geom_mark_ellipse(aes(color=group))

    ggsave(paste(prefix, ".nmds.png", sep=""),  p, width=120, height=100, units="mm")
    ggsave(paste(prefix, ".nmds.pdf", sep=""), p, width=120, height=100, units="mm")

}


add_help_args <- function(args){

  if(length(args) != 3) {
    cat("Version: v1.0.0\n")
    cat("Author:Boya Xu\n")
    cat("Email:xby@bioyegene.com\n")
    cat("Function:Draw NMDS Analysis Picture.\n")
    cat("Example:plot_nmds.r abundance_species.xls group.list prefix\n")
    cat("Input file format:\n")
    cat("Tax Id\tsample1\tsample2\tsample3\n")
    cat("OTU1\t10\t9\t9\n")
    quit()
  }
}


args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)

plot_nmds(args[1], args[2], args[3])
