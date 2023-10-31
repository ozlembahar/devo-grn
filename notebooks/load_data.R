library(SeuratDisk)
library(dplyr)
library(Seurat)
library(Matrix)

sdata <- LoadH5Seurat("./data/SeuratProject.h5Seurat")
# # counts_df <- sdata@assays$RNA@counts %>% as.matrix %>% t %>% as.data.frame

# # write.csv(counts_df, "data/expression_mat.csv", row.names = TRUE)

# cell_type <- as.data.frame(sdata@active.ident)
# colnames(cell_type) <- c("type")
# write.csv(cell_type, "./data/cell_type.csv", row.names = TRUE)

# anndata <- Convert("./data/SeuratProject.h5Seurat", dest= "h5ad")  # convert Seurat to AnnData to be analysed via Scipy

#### subset ####
# subset <- subset(x=sdata, idents = "MG", invert=TRUE) 
# subset_counts_df <- subset@assays$RNA@counts %>% as.matrix %>% t %>% as.data.frame
# write.csv(subset_counts_df, "data/expression_mat_noMG.csv", row.names = TRUE)


######
# convert Seurat to python 

# DimPlot(sdata, reduction='umap', label=TRUE)

# counts_matrix = GetAssayData(sdata, assay='RNA', slot='counts')
# writeMM(counts_matrix, file=paste0(file='./data/to_py/matrix.mtx'))

# write.csv(sdata@reductions$pca@cell.embeddings, file='./data/to_py/pca.csv', quote=F, row.names=F)
# write.table(data.frame('gene' = rownames(counts_matrix)), file='./data/to_py/gene_names.csv', quote = F, row.names=F, col.names = F)

# sdata$barcode <- colnames(sdata)
# sdata$UMAP_1 <- sdata@reductions$umap@cell.embeddings[,1]
# sdata$UMAP_2 <- sdata@reductions$umap@cell.embeddings[,2]
# write.csv(sdata@meta.data, file="./data/to_py/metadata.csv", quote=F, row.names=F)


# write.csv(sdata@active.ident, file='./data/to_py/types.csv', quote=F, row.names=F)


