library(AUCell)
library(SCENIC)
library("reticulate")
library(tidyr)
library(dplyr)
library(GSEABase)

final_regulons <- read.csv("results/results_041023/final_regulons_forR.csv")
# rownames(final_regulons) = final_regulons[,1]

# final_regulons <- final_regulons %>% subset(select=-c(TF))

final_regulons_list = apply(final_regulons,1, function(x) list(as.character(x['gene'])))

genes_list <- do.call(c, final_regulons_list)
names(genes_list) = final_regulons$TF
#### load exp. Data
library(SeuratDisk)
library(Seurat)
sdata <- LoadH5Seurat("./data/SeuratProject.h5Seurat")
# exprMatrix = sdata@assays$RNA@counts
exprMatrix = as.matrix(read.csv("/Users/danabarilan/Documents/forMasters/Intenship/code/devo-grn/data/expression_mat.csv"))
rownames(exprMatrix) <- exprMatrix[, 1]  ## set rownames
exprMatrix <- exprMatrix[, -1]           ## remove the first variable

#geneSets <- AUCell::subsetGeneSets(genes_list, rownames(exprMatrix)) 
#cbind(nGenes(geneSets))


cells_AUC <- AUCell::AUCell_run(exprMatrix, genes_list). # does not work
#########################

