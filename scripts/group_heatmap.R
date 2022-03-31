#!/usr/bin/env Rscript

library(grid)
library(circlize)
library(ComplexHeatmap)
options(bitmapType='cairo')

add_help_args <- function(args){

    if(length(args) != 3) {
        cat("Version: v1.0.0\n")
        cat("Author:Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:Plot grouped heatmaps.\n")
        cat("Example:group_heatmap.R abundance_zscore_top50.xls group.list group.heatmap\n")
        cat("Input file format:\n")
        cat("Taxid\tsample1\tsample2\tsample3\n")
        cat("OTU1\t0.276046\t-1.339209\t1.063164")
        cat("group.lsit format:\n")
        cat("Sample\tGroup\n")
        cat("sample1\tgroup1")
        quit()
    }

}


read_group <- function(group){
   
}

args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)

file <- args[1]
group <- args[2]
prefix <- args[3]

data1 <- read.table(file, header=T, row.names=1, stringsAsFactors=F) #读取数据
data <- as.matrix(data1)	#将数据转化为矩阵
rows <- nrow(data)	#计算行数
columns <- length(data[1,])	#计算列数
group <- read.delim(group, row.names=1, sep="\t", head=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
colnames(group) <- c("group")

if(rows<=20) {
    heights <- 6
}else if(rows<=50) {
    heights <- 14
}else {
    heights <- 18
}

if(columns<=4) {
    widths <- 6
}else if(columns<=8) {
    widths <- 10
}else {
    widths <- 14
}

pdf(paste(prefix, ".pdf", sep=""), width=widths, height=heights)
a <- dev.cur()
png(paste(prefix, ".png", sep=""), width=widths*60, height=650)

ha <- HeatmapAnnotation(df=group)

format <- colorRamp2(seq(min(data), max(data), length=3), c("blue", "#EEEEEE", "red"), space="RGB")
dev.control("enable")
Heatmap(data, name=" ", col=format,  top_annotation=ha) #show_column_dend=FALSE)#column_names_side = "top", column_names_rot=0) #column_names_rot=0旋转字体方向
dev.copy(which=a)
dev.off()
dev.off()
