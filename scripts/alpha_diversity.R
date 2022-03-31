#!/usr/bin/env Rscript

library(vegan)
library(picante)
options(bitmapType='cairo') #关闭服务器与界面的互动响应


stat_alpha_diversity <- function(data, tree=NULL){

    stat_data <- estimateR(data)

    richness <- stat_data[1, ]
    chao1 <- stat_data[2, ]
    ACE <- stat_data[4, ]

    shannon <- diversity(data, index="shannon", base=2)
    simpson <- diversity(data, index="simpson")#注意，这里是Gini-Simpson 指数
    coverage <- (1-rowSums(data<1)/rowSums(data))*100

    #保留四位小数
    shannon <- sprintf("%0.4f", shannon)
    simpson <- sprintf("%0.4f", simpson)
    coverage <- sprintf("%0.4f", coverage)

    result <- data.frame(Richness=richness, Shannon=shannon, Simpson=simpson,
                         Chao1=chao1, ACE=ACE, Coverage=coverage)

    if (!is.null(tree)){
        pd_tree <- pd(data, tree, include.root=FALSE)[1]
        names(pd_tree) <- 'pd_tree'
        result <- cbind(result, pd_tree)
        result <- data.frame(Richness=richness, Shannon=shannon, Simpson=simpson, 
                             Chao1=chao1, ACE=ACE, Coverage=coverage, Pielou=pd_tree,
                             Coverage=coverage)
    }

    return(result)
}


add_help_args <- function(args){

    if(length(args) == 0) {
        cat("Version: v1.0.0\n")
        cat("Author:Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:Alpha Diversity Analysis.\n")
        cat("Example:alpha_diversity.R abundance.xls\n")
        cat("Input file format:\n")
        cat("#Tax Id\tsample1\tsample2\tsample3\n")
        cat("OTU1\t934\t932\t910")
        cat("Input must be an integer")
        quit()
    }

}


args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)

data <- read.delim(args[1], row.names=1, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
data <- t(data)

tree <- NULL
if(length(args) == 2) {
    tree <- read.tree(args[2])
}

result <- stat_alpha_diversity(data, tree)

print(result)
