library(SeuratDisk)
library(dplyr)

sdata <- LoadH5Seurat("./data/SeuratProject.h5Seurat")
# counts_df <- sdata@assays$RNA@data %>% as.matrix %>% t %>% as.data.frame

# write.csv(counts_df, "data/expression_mat.csv", row.names = TRUE)

cell_type <- as.data.frame(sdata@active.ident)
colnames(cell_type) <- c("type")
write.csv(cell_type, "./data/cell_type.csv", row.names = TRUE)

