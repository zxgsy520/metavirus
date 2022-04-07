#!/usr/bin/env Rscript

library(vegan)
library(RColorBrewer)

options(bitmapType='cairo') #关闭服务器与界面的互动响应


get_color <- function(samples){

    #根据样本数目产生对应得颜色列表
    sam_num = length(samples)
    if (sam_num <= 12){
        colors <- brewer.pal(n=sam_num, name="Paired")
    }else{
        colors <- brewer.pal(n=12, name="Paired")
        pal <- colorRampPalette(colors)
        colors <- pal(sam_num)
    }
    return(colors)
}


plot_anosim <- function(abundance, group, prefix){

    #anosim画图
    otu <- read.delim(abundance, sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
    otu <- data.frame(t(otu), stringsAsFactors = FALSE) #将丰度表倒置
    otu <- otu[rowSums(otu[])>0,]
    group <- read.delim(group, sep="\t", stringsAsFactors=FALSE, header=TRUE)
    colnames(group) <- c("sample", "group")
    sample_id <- row.names(otu)
    part_group <- data.frame(group, stringsAsFactors = FALSE)

    group <- part_group[which (part_group[,1] %in% sample_id),]
    colors <- get_color(unique(group$group))
    anosim_result_otu <- anosim(otu, group$group, permutations=999, distance="bray") #跑anosim结果

    pdf(paste(prefix, ".anosim.pdf", sep=""), width=14, height=7)
    ptemp <- dev.cur()
    png(paste(prefix, ".anosim.png", sep=""), width=800, height=400) #画所有分组图
    dev.control("enable")

    plot(anosim_result_otu, xlab="Sample grouping", ylab="Bray-curtis distance", col=colors, bty="l")
    dev.copy(which=ptemp)
    dev.off() 
    dev.off()

    group_name <- unique(group$group)
    anosim_result_two <- NULL
    print(group_name)
    for (i in 1:(length(group_name)-1)) {  #两两对比画图
        for (j in (i+1):length(group_name)) {
            group_ij <- subset(group, group %in% c(group_name[i], group_name[j]))
            otu_ij <- otu[group_ij$sample, ]
            anosim_result_otu_ij <- anosim(otu_ij, group_ij$group, permutations = 999, distance = 'bray')	#随机置换检验 999 次
            
            print(paste(group_name[i], group_name[j], sep="|"))
            if (anosim_result_otu_ij$signif <= 0.001) Sig <- "***"
            else if (anosim_result_otu_ij$signif <= 0.01) Sig <- "**"
            else if (anosim_result_otu_ij$signif <= 0.05) Sig <- "*"
            else Sig <- NA
            anosim_result_two <- rbind(anosim_result_two, c(paste(group_name[i], group_name[j], sep="|"), "Bray-Curtis", anosim_result_otu_ij$statistic, anosim_result_otu_ij$signif, Sig))

            #每次循环输出图片
            pdf(paste(prefix, ".", group_name[i], "_", group_name[j], ".anosim.pdf", sep=""), width=12, height=5)
            ptemp <- dev.cur()
            png(paste(prefix, ".", group_name[i], "_", group_name[j], ".anosim.png", sep=""), width=600, height=300)
            plot(anosim_result_otu_ij, xlab="Sample grouping", ylab="Bray-curtis distance", col=c("gray","red","blue"), bty="7")
            dev.copy(which=ptemp)
            dev.off()
            dev.off()
        }
    }
    return(anosim_result_two)
}


run_anosims <- function(abundance, group, prefix){

    anosim_result_two <- plot_anosim(abundance, group, prefix)
    anosim_result_two <- data.frame(anosim_result_two, stringsAsFactors=FALSE)
    names(anosim_result_two) <- c("group", "distance", "R", "P_value", "Sig")

    write.table(anosim_result_two,  paste(prefix, ".ANOSIM.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE, na ="")
}


add_help_args <- function(args){

    if(length(args) != 3) {
        cat("Version: v1.0.0\n")
        cat("Author:Boya Xu, Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:anosim analysis.\n")
        cat("Example:anosim_group.r abundance_species.tsv group.listprefix\n")
        cat("Input file format:\n")
        cat("Tax Id\tsample1\tsample2\tsample3\n")
        cat("OTU1\t10\t9\t9\n")
        quit()
    }
}


args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
run_anosims(args[1], args[2], args[3])
