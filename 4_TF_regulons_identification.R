######## figure5 TF regulons identification of day 16 epithelial subtypes and heatmap
library(SCENIC)
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(ggsci)
library(circlize)
library(stringr)
library(writexl)

# In vitro 
load(file = 'epi.RData')
d16 <- counts[, counts$day %in% 'D16'] 
loom <- open_loom('out_SCENIC.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")   
regulons <- regulonsToGeneLists(regulons_incidMat) 
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)     
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])   
embeddings <- get_embeddings(loom)  
close_loom(loom)
sub_regulonAUC <- regulonAUC[,match(colnames(d16),colnames(regulonAUC))]
identical(colnames(sub_regulonAUC), colnames(d16))
cellTypes <- data.frame(row.names = colnames(d16), 
                        celltype = d16$cell.type)
selectedResolution <- "celltype"
cellsPerGroup <- split(rownames(cellTypes),cellTypes[,selectedResolution])
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)
regulonActivity_byGroup <- sapply(cellsPerGroup,function(cells)rowMeans(getAUC(sub_regulonAUC)[,cells]))
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
melt <- melt(regulonActivity_byGroup_Scaled)
melt %>%
  group_by(Var2) %>%
  top_n(n = 10, wt = value) -> top10
plot.data <- regulonActivity_byGroup_Scaled[top10$Var1,]
p1 <- Heatmap(
  plot.data,
  name                         = "z-score",
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = T,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

# In vivo
load(file = 'epi.subset.rdata')
output <- epi.subset[, epi.subset$subtype %in% c('Fetal_stomach_fundus','Fetal_stomach_antrum')] 
week14 <- output[,output$Age_week %in% '14'] 
output <- week14

invivo_loom <- open_loom('out_SCENIC.loom')
invivo_regulons_incidMat <- get_regulons(invivo_loom, column.attr.name="Regulons")   
invivo_regulons <- regulonsToGeneLists(invivo_regulons_incidMat) 

invivo_regulonAUC <- get_regulons_AUC(invivo_loom,column.attr.name='RegulonsAUC')
invivo_regulonAucThresholds <- get_regulon_thresholds(invivo_loom)     
tail(invivo_regulonAucThresholds[order(as.numeric(names(invivo_regulonAucThresholds)))])   
close_loom(invivo_loom)
invivo_sub_regulonAUC <- invivo_regulonAUC[,match(colnames(output),colnames(invivo_regulonAUC))]
identical(colnames(invivo_sub_regulonAUC), colnames(output))

invivo_cellTypes <- data.frame(row.names = colnames(output), 
                               celltype = output$subtype)
selectedResolution <- "celltype"
invivo_cellsPerGroup <- split(rownames(invivo_cellTypes),invivo_cellTypes[,selectedResolution])
invivo_sub_regulonAUC <- invivo_sub_regulonAUC[onlyNonDuplicatedExtended(rownames(invivo_sub_regulonAUC)),] 
dim(invivo_sub_regulonAUC)
invivo_regulonActivity_byGroup <- sapply(invivo_cellsPerGroup,function(cells)rowMeans(getAUC(invivo_sub_regulonAUC)[,cells]))
invivo_regulonActivity_byGroup_Scaled <- t(scale(t(invivo_regulonActivity_byGroup),
                                                 center = T, scale=T)) 

length(intersect(names(regulons),names(invivo_regulons)))
common_regulon <- intersect(names(regulons),names(invivo_regulons))

plot.data <- invivo_regulonActivity_byGroup_Scaled[common_regulon,]
p1 <- pheatmap(plot.data,cluster_rows = T,name = 'z-score',show_rownames = T,cluster_cols = F)
row_cluster <- cutree(p1$tree_row, k=2)
invivoOrder <- as.data.frame(plot.data[p1$tree_row$order,])
invivoOrder$cluster <- as.character(row_cluster[match(rownames(invivoOrder), names(row_cluster))])

plot.data <- regulonActivity_byGroup_Scaled[common_regulon,]
p2 <- pheatmap(plot.data,cluster_rows = T,name = 'z-score',show_rownames = T,cluster_cols = F)
row_cluster <- cutree(p2$tree_row, k=2)
Order <- as.data.frame(plot.data[p2$tree_row$order,])
Order$cluster <- as.character(row_cluster[match(rownames(Order), names(row_cluster))])

fundus_tf <- intersect(rownames(subset(Order,cluster %in% '2')),rownames(subset(invivoOrder,cluster %in% '1')))
antrum_tf <- intersect(rownames(subset(Order,cluster %in% '1')),rownames(subset(invivoOrder,cluster %in% '2')))
tf <- c(fundus_tf,antrum_tf)

colnames(invivo_regulonActivity_byGroup) <- c('Fetal_stomach_antrum','Fetal_stomach_fundus')
fundus_top10 <- invivo_regulonActivity_byGroup[fundus_tf,]|>
  as.data.frame()|>
  arrange(desc(Fetal_stomach_fundus))
antrum_top10 <- invivo_regulonActivity_byGroup[antrum_tf,]|>
  as.data.frame()|>
  arrange(desc(Fetal_stomach_antrum))

plot <- c(rownames(fundus_top10)[1:10],rownames(antrum_top10)[1:10])
invivo_regulonActivity_byGroup_Scaled <- subset(invivo_regulonActivity_byGroup_Scaled,select=c(2,1)) 
# Figure R25
pheatmap(invivo_regulonActivity_byGroup_Scaled[plot,],cluster_rows = F,name = 'z-score',show_rownames = T,cluster_cols = F,show_colnames = F)
pheatmap(regulonActivity_byGroup_Scaled[plot,],cluster_rows = F,name = 'z-score',show_rownames = T,cluster_cols = F,show_colnames = F)







