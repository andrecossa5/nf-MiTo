
#!/usr/bin/env Rscript

library(argparse)
library(dplyr)
library(reshape2)
library(Matrix)
library(ape)
library(PATH)


##


# cov <- 'GBC'
# sample <- 'MDA_clones'
# path_tree <- '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/results/MI_TO_bench/phylo_inference/final_trees/MDA_clones/job1/annotated_tree.newick'
# path_meta <- '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/cells_meta.csv'

parser <- ArgumentParser(description = 'PATH wrapper')

parser$add_argument("path_tree", type = "character", help = "Path to .newick tree")
parser$add_argument("path_meta", type = "character", help = "Path to cells meta data")
parser$add_argument("lineage_column", type = "character", help = "Lineage column to evaluate phylocorr")
args <- parser$parse_args()


##


# Processing
tree <- ape::read.tree(args$path_tree)
cov_df <- read.csv(args$path_meta, row.names = 1)[tree$tip.label,]
x_factor <- as.factor(cov_df[[args$lineage_column]])
X <- state2mat.sparse(x_factor)

Winv <- inv.tree.dist(tree, node=TRUE, norm=FALSE)
modxcor <- xcor(X, Winv)

Idf <- reshape2::melt(modxcor$Morans.I, value.name = "I")
Zdf <- reshape2::melt(modxcor$Z.score, value.name = "Z")
pdf <- reshape2::melt(modxcor$one.sided.pvalue, value.name = "p")
dfs <- list(Idf, Zdf, pdf)

# Save
results <- Reduce(function(x, y) full_join(x, y, by = c("Var1", "Var2")), dfs)
results$Var1 <- factor(results$Var1, levels = 1:length(levels(x_factor)), labels = levels(x_factor)) 
results$Var2 <- factor(results$Var2, levels = 1:length(levels(x_factor)), labels = levels(x_factor)) 

write.csv(results, 'phylocorr.csv')


##

