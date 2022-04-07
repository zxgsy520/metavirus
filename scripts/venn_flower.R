#!/usr/bin/env Rscript

library(rio)
library(VennDiagram)
library(plotrix)
library(RColorBrewer)

options(bitmapType='cairo') #关闭服务器与界面的互动响应


get_color <- function(samples){

    #根据样本数目产生对应得颜色列表
    sam_num = length(samples)
    if (sam_num <= 2){
        colors <- c("#8DD3C7", "#FFFFB3")
    }else if(sam_num <= 12){
        colors <- brewer.pal(n=sam_num, name="Paired")
    }else{
        colors <- brewer.pal(n=12, name="Paired")
        pal <- colorRampPalette(colors)
        colors <- pal(sam_num)
    }

    return(colors)
}


stat_core_otu <- function(data){

    data <- data[which(rowSums(data) > 0),]   #过滤所有丰度都为0得行
    otus <- rownames(data)
    otus <- paste("OTU", 1:length(otus), sep="")
    rownames(data) <- otus
    samples <- colnames(data)

    core_otu <- otus
    assign(samples[1], otus)
    temp <- list()
    cores <- c(temp, list(get(samples[1])))
    otu_num <- length(otus)

    for (i in 2:ncol(data)) {
      otus <- unique(rownames(data[data[,i]>0,]))

      assign(samples[i], otus)
      cores <- c(cores, list(get(samples[i])))
      core_otu <- intersect(core_otu, otus)
      otu_num <- c(otu_num, length(otus))
    }

    core_num <- length(core_otu)
    return(list(cores, otu_num, core_num))

}


flower_plot <- function(sample, otu_num, core_otu, ellipse_col, prefix="A1",
                        start=90, a=0.5, b=2, r=1, circle_col="#cafafa"){

    #绘制花瓣图
    pdf(paste(prefix, ".venn_flower.pdf", sep=""), width = 15, height = 15)
    ptemp <- dev.cur()
    png(paste(prefix, ".venn_flower.png", sep=""), width=600, height=600)
    dev.control("enable")

    par(bty='n', ann=F, xaxt='n', yaxt='n', mar=c(1,1,1,1))
    plot(c(0,10), c(0,10), type='n')

    n <- length(sample)
    deg <- 360/n

    res <- lapply(1:n, function(t){

        draw.ellipse(x = 5 + cos((start + deg*(t-1))*pi/180),
                     y = 5 + sin((start + deg*(t-1))*pi/180),
                     col = ellipse_col[t],
                     border = ellipse_col[t],
                     a = a, b = b, angle = deg*(t-1)
        )

        text(x = 5 + 2.5*cos((start + deg*(t-1))*pi/180),
             y = 5 + 2.5*sin((start + deg*(t-1))*pi/180),
             otu_num[t], cex = 1.5
        )

        if (deg*(t-1) < 180 && deg*(t-1) > 0 ){
            text(x = 5 + 3.3*cos((start + deg*(t-1))*pi/180),
                 y = 5 + 3.3*sin((start + deg*(t-1))*pi/180),
                 sample[t],
                 srt = deg*(t-1)-start,
                 adj = 1,
                 cex = 1.5
            )
        }else{
            text(x = 5 + 3.3*cos((start + deg*(t-1))*pi/180),
                y = 5 + 3.3*sin((start + deg*(t-1))*pi/180),
                sample[t],
                srt = deg*(t-1)+start,
                adj = 0,
                cex = 1.5
            )
        }
    })

    draw.circle(x=5, y=5, r=r, col=circle_col, border=NA)
    text(x=5, y=5, paste("Core:", core_otu), cex=1.5)
    dev.copy(which=ptemp)
    dev.off()
    dev.off()
}


venn_flower <- function(data, prefix){

    data <- read.table(data, sep="\t", row.names=1, head=TRUE, check.names=FALSE, quote="")
    samples <- colnames(data)
    colors <- get_color(samples)
    result <- stat_core_otu(data)

    if (length(samples) <= 6){
        venn <- venn.diagram(x=result[[1]], category.names=samples, fill=colors,
                             filename=NULL, col="black", lwd=3, cex=2,
                             cat.cex=2.5, output=TRUE)

        pdf(paste(prefix, ".venn_flower.pdf", sep=""), width = 15, height = 15)
        ptemp <- dev.cur()
        png(paste(prefix, ".venn_flower.png", sep=""), width=600, height=600)
        dev.control("enable")

        grid.draw(venn)

        dev.copy(which=ptemp)
        dev.off()
        dev.off()
    }else{
        flower_plot(sample=samples, otu_num=result[[2]], core_otu=result[[3]],
                    ellipse_col=colors, prefix=prefix)
    }

}


add_help_args <- function(args){

    if(length(args) != 2) {
        cat("Version: v1.0.0\n")
        cat("Author:Boya Xu, Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:Draw venn or flower diagram.\n")
        cat("Example:venn_flowerr.R abundance_species.tsv prefix\n")
        cat("Input file format:\n")
        cat("Tax Id\tsample1\tsample2\tsample3\n")
        cat("OTU1\t10\t9\t9")
        quit()
    }
}


args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)

venn_flower(args[1], args[2])
