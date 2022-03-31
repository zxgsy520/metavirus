#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
options(bitmapType='cairo') #关闭服务器与界面的互动响应

GROUP_COLORS <- c("red", "#CC6699", "#00FFFF", "blue", "#BA55D3", "#9400D3",
                  "#4B0082", "#9370DB", "#6A5ACD", "#483D8B", "#0000CD", "#00008B",
                  "#4169E1", "#B0C4DE")

COLORS <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
            "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "gray","#FFB6C1",
            "#3A5FCD","#7CCD7C","#8A2BE2","#CDAD00","#EEEE00",
            "#8B6914","#FF00FF","#8F8F8F","#698B22","#87CEFF","#C71585","#EE9A49",
            "#00FFFF","#FFB6C1","#DC143C","#DB7093","#FF1493","#DA70D6","#8B008B",
            "#BA55D3","#9400D3","#4B0082","#9370DB","#6A5ACD","#483D8B","#0000CD",
            "#00008B","#4169E1","#B0C4DE","#FF0000")


add_help_args <- function(args){

    if(length(args) != 4) {
        cat("Version: v1.1.0\n")
        cat("Author:Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:Plot a cluster histogram.\n")
        cat("Example:plot_tree_bar.R abundance_species_top20.tsv group.lsit level prefix\n")
        cat("Input file format:\n")
        cat("\tsample1\tsample2\tsample3\n")
        cat("OTU1\t9.276046\t9.339209\t9.063164")
        cat("group.lsit format:\n")
        cat("Sample\tGroup\n")
        cat("sample1\tgroup1")
        quit()
    }

}


treeline <- function(pos1, pos2, height, col1, col2){ #聚类树绘制，按分组给分支上色

    meanpos <- (pos1[1] + pos2[1]) / 2
    segments(y0 = pos1[1] - 0.4, x0 = -pos1[2], y1=pos1[1] - 0.4, x1= -height,  col = col1,lwd = 2)
    segments(y0 = pos1[1] - 0.4, x0 = -height,  y1 = meanpos - 0.4, x1 = -height,  col = col1,lwd = 2)
    segments(y0 = meanpos - 0.4, x0 = -height,  y1 = pos2[1] - 0.4, x1 = -height,  col = col2,lwd = 2)
    segments(y0 = pos2[1] - 0.4, x0 = -height,  y1 = pos2[1] - 0.4, x1 = -pos2[2], col = col2,lwd = 2)
}


maxlength <- function(names){

    maxlen <- 0
    for (i in names){
        if (nchar(i) >= maxlen){
            maxlen <- nchar(i)
        }
    }
    return(maxlen)

}


zoom_ratio <- function(names){
 
    zor <- 2.4   
    maxlen <- maxlength(names)

    if(maxlen <= 3){
        zor <- 2.4
    }else if(maxlen <= 5){
        zor <- 2.3
    }else if(maxlen <= 8){
        zor <- 2.2
    }else{
        zor <- 2.1
    }
    return(zor)
}


plot_tree_bar <- function(file, group, level, prefix){

    data <- read.delim(file, row.names=1, sep="\t", head=TRUE, check.names=FALSE)
    data <- dplyr::filter(data, !grepl('unclassified|other', data[,1])) #过滤掉没分类的和分类到其他的

    dis_bray <- vegan::vegdist(t(data), method="bray") #计算各样本直接的距离
    tree <- hclust(dis_bray, method="average")  #进行乘次聚类
    plot(tree)

    if(ncol(data) <= 8){
        dheight <- 6
        pheight <- 400
    }else if(ncol(data) <= 12){
        dheight <- 8
        pheight <- 600
    }else if(ncol(data) <= 18){
        dheight <- 10
        pheight <- 800
    }else{
        dheight <- 14
        pheight <- 1000
    }

    
    group <- read.delim(group, row.names=1, sep="\t", head=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
    #group <- read.table(group, row.names=1, header=T, sep="\t")
    colnames(group) <- c("group")
    
    grp <- t(group[1])
    group_col <- GROUP_COLORS[1:length(grp)]
    names(group_col) <- c(1:length(grp)) #为每种颜色名命为数字
    group_name <- unique(t(grp))

    #根据分组组名生成一列数字，同组的数字相同，为了方便同组的上同种颜色
    for(step in 1:length(group_name)){
        group$group_num[group$group==group_name[step]] <- step
    }
    grp <- group[2]

    pdf(paste(prefix, paste0(level, ".pdf"), sep="_"), width=18, height=dheight, onefile= FALSE)
    a <- dev.cur()
    png(paste(prefix, paste0(level, ".png"), sep="_"), width=1100, height=pheight)
    dev.control("enable")

    layout(t(c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3)))
    par(mar = c(5, 5, 5, 0))

    plot(0, type='n', xaxt='n', yaxt='n', frame.plot=FALSE, xlab='', ylab='',
        xlim = c(-max(tree$height), 0), ylim=c(0, length(tree$order))
    )
    legend("topleft", legend=group_name, pch=15, col=group_col, bty="n", cex=2.4)

    meanpos <- matrix(rep(0, 2*length(tree$order)), ncol=2)
    meancol <- rep(0, length(tree$order))

    for(step in 1:nrow(tree$merge)) {
        if(tree$merge[step, 1] < 0){
            pos1 <- c(which(tree$order == -tree$merge[step, 1]), 0)
            col1 <- group_col[as.character(grp[tree$labels[-tree$merge[step, 1]],1])]
        }else{
            pos1 <- meanpos[tree$merge[step, 1], ]
            col1 <- meancol[tree$merge[step, 1]]
        }
        if (tree$merge[step, 2] < 0) {
            pos2 <- c(which(tree$order == -tree$merge[step, 2]), 0)
            col2 <- group_col[as.character(grp[tree$labels[-tree$merge[step, 2]],1])]
        } else {
            pos2 <- meanpos[tree$merge[step, 2], ]
            col2 <- meancol[tree$merge[step, 2]]
        }
        height <- tree$height[step]
        treeline(pos1, pos2, height, col1, col2)
        meanpos[step, ] <- c((pos1[1] + pos2[1]) / 2, height)

        if (col1 == col2) {
            meancol[step] <- col1
        } else {
            meancol[step] <- "grey"
        }
    }

    ##堆叠柱形图
    #样本顺序调整为和聚类树中的顺序一致
    data <- data[ ,tree$order]
    names(COLORS) <- rownames(data)
    interval <- maxlength(colnames(data))+2

    par(mar = c(5, interval, 5, 0)) #堆叠柱形图
    py <- barplot(as.matrix(data), col=COLORS, space=1, width=0.5,
                 cex.axis = 2.3, horiz=TRUE, cex.lab=3, xlab="Relative Abundance",
                 yaxt="n", las=0.5, ylim=c(0, ncol(data)), #family="mono",
        )

    zor <- zoom_ratio(colnames(data))
    text(x=-interval*1.2, y=py, labels=colnames(data),
         col=group_col[group[tree$order, 2]], cex=zor, xpd=TRUE
    )

    par(mar=c(5, 0, 5, 0)) #柱形图图例
    plot(0, type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
    legend("left", pch=15, col=COLORS, legend=rownames(data),
           bty="n", x.intersp=0.5, y.intersp=0.9, cex=2.3
    ) #设置图框与文字的距离x.intersp, 文字之前的距离, 字体大小cex。

    dev.copy(which=a)
    dev.off()
    dev.off()
}


args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)
plot_tree_bar(args[1], args[2], args[3], args[4])
