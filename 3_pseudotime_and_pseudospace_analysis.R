######## differentiation states among cells of epithelial, neural, and mesenchymal lineages
######## CytoTRACE  -------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(ggsci)
library(CytoTRACE)
library(readxl)

load(file = './integrated.RData')
col=pal_igv('default',alpha = 1)(51)

# mesenchymal lineage
sce <- sample.integrated[, sample.integrated$cell_type %in% c('hPSC','Mesenchymal1','Mesenchymal2','Mesenchymal3')]
DefaultAssay(sce) <- 'integrated'
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 
sce <- RunUMAP(sce, dims = 1:20)
DimPlot(sce, reduction = 'umap')

data <- as.matrix(GetAssayData(sce, assay = "RNA", slot = "counts"))
results <- CytoTRACE(data, ncores = 12, subsamplesize = 1000)
pheno <- sce$lab
pheno <- as.character(pheno)
names(pheno) <- rownames(sce@meta.data)
VEC=sce@reductions$umap@cell.embeddings
# Extended Data Figure8
plotCytoTRACE(results, phenotype = pheno,emb = VEC)

# Supplementary Figure5 Ecm gene heatmap
highlight <- read_excel('ECM_gene.xlsx') 
DefaultAssay(sce) <- 'RNA'
sce <- ScaleData(sce)
mt <- as.data.frame(t(as.matrix(GetAssayData(sce, assay = "RNA", slot = "scale.data"))))
group.by <- 'cell_type'
mt <- aggregate(mt, by=list(sce@meta.data[[group.by]]), FUN="mean")
rownames(mt) <- mt$Group.1
mt <- t(mt[,-1])
highlight <- highlight[order(-highlight[,'mesen.cytotrce.cor']),]
cts_mesen <- as.matrix(mt[highlight$Gene,])
pheatmap(cts_mesen,show_colnames =T,show_rownames = T,color = viridis(8),
         cluster_rows = F,
         cluster_cols = F,
         name= 'Scaled Expression')

# Extended Data Figure8 Absolute time and relative time
data <- sce[,sce$day %in% c('D4','D7','D10','D13','D16')]@meta.data[,c('cell.type','day','cyto_pseudotime')]

ggplot(data,aes(x=cell.type,y=cyto_pseudotime))+
  geom_violin(width =0.8,fill='grey90',color='grey90')+
  geom_quasirandom(aes(color=day),width = 0.2,size=0.5)+ 
  stat_summary(aes(fill=day), geom="point",
               fun = mean, shape=21, size=2.5,stroke=0.5)+ 
  scale_color_manual(name = 'day',
                     values = c('#57C3F3',"#E5D2DD","#53A85F","#F1BB72",'#E63863'),
                     labels = c('D4','D7','D10','D13','D16'))+ 
  scale_fill_manual(values = c('#57C3F3',"#E5D2DD","#53A85F","#F1BB72",'#E63863'),
                    guide = 'none')+ 
  scale_y_continuous(limits = c(0,1))+
  theme_classic(base_size = 15)+
  coord_cartesian(clip = 'off')+ 
  geom_point(aes(x=2.5,y=6),pch=21,size=4,stroke=1.3)+ 
  geom_point(aes(x=5,y=0.7),pch=21,size=2.5,stroke=0.5)+ 
  annotate(geom = 'text',label = 'Mean',x=5.2,y=0.7,size=4)+ 
  theme(plot.margin = margin(10,50,10,10), 
        legend.position = c(1,0.3),
        axis.text = element_text(color='black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1)) +ylab('Predicted ordering')+xlab('Cell type')

saveRDS(sce,file = 'mesen.rds')

# neural lineage
sce <- sample.integrated[, sample.integrated$cell_type %in% c('hPSC','Early NPC','Premigratory ENC','NPC','Migratory ENCC','Neuron')]
DefaultAssay(sce) <- 'integrated'
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 
sce <- RunUMAP(sce, dims = 1:20)
DimPlot(sce, reduction = 'umap')

data <- as.matrix(GetAssayData(sce, assay = "RNA", slot = "counts"))
results <- CytoTRACE(data, ncores = 12, subsamplesize = 1000)
pheno <- sce$lab
pheno <- as.character(pheno)
names(pheno) <- rownames(sce@meta.data)
VEC=sce@reductions$umap@cell.embeddings
# Extended Data Figure8
plotCytoTRACE(results, phenotype = pheno,emb = VEC)

# Supplementary Figure5 Ecm gene heatmap
highlight <- read_excel('ECM_gene.xlsx') 
DefaultAssay(sce) <- 'RNA'
sce <- ScaleData(sce)
mt <- as.data.frame(t(as.matrix(GetAssayData(sce, assay = "RNA", slot = "scale.data"))))
group.by <- 'cell_type'
mt <- aggregate(mt, by=list(sce@meta.data[[group.by]]), FUN="mean")
rownames(mt) <- mt$Group.1
mt <- t(mt[,-1])
highlight <- highlight[order(-highlight[,'neuron.cytotrce.cor']),]
cts_neuron <- as.matrix(mt[highlight$Gene,])
pheatmap(cts_neuron,show_colnames =T,show_rownames = T,color = viridis(8),
         cluster_rows = F,
         cluster_cols = F,
         name= 'Scaled Expression')

# Extended Data Figure8 Absolute time and relative time
sample <- subset(sce,downsample=2000)
data <- sample[,sample$day %in% c('D4','D7','D10','D13','D16')]@meta.data[,c('cell.type','day','cyto_pseudotime')]

data <- sce[,sce$day %in% c('D4','D7','D10','D13','D16')]@meta.data[,c('cell.type','day','cyto_pseudotime')]

ggplot(data,aes(x=cell.type,y=cyto_pseudotime))+
  geom_violin(width =0.8,fill='grey90',color='grey90')+
  geom_quasirandom(aes(color=day),width = 0.2,size=0.2)+ 
  stat_summary(aes(fill=day), geom="point",
               fun = mean, shape=21, size=2.5,stroke=0.5)+ 
  scale_color_manual(name = 'day',
                     values =  c('#57C3F3',"#E5D2DD","#53A85F","#F1BB72",'#E63863'),
                     labels = c('D4','D7','D10','D13','D16'))+ 
  scale_fill_manual(values =  c('#57C3F3',"#E5D2DD","#53A85F","#F1BB72",'#E63863'),
                    guide = 'none')+ 
  scale_y_continuous(limits = c(0,1))+
  theme_classic(base_size = 15)+
  coord_cartesian(clip = 'off')+ 
  geom_point(aes(x=2.5,y=6),pch=21,size=4,stroke=1.3)+
  theme(plot.margin = margin(10,50,10,10), 
        legend.position = c(1,0.3),
        axis.text = element_text(color='black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))+ 
  geom_point(aes(x=6.8,y=0.7),pch=21,size=2.5,stroke=0.5)+ 
  annotate(geom = 'text',label = 'Mean',x=7.1,y=0.7,size=4) +ylab('Predicted ordering')+xlab('Cell type')

saveRDS(sce,file = 'neuron.rds')

# epithelial lineage
sce <- sample.integrated[, sample.integrated$cell.type %in% c('hPSC','DE','Partial Epi','Epithelium','Enteroendocrine')]
DefaultAssay(sce) <- 'integrated'
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 
sce <- RunUMAP(sce, dims = 1:20)
DimPlot(sce, reduction = 'umap')

data <- as.matrix(GetAssayData(sce, assay = "RNA", slot = "counts"))
results <- CytoTRACE(data, ncores = 12, subsamplesize = 500)
pheno <- sce$lab
pheno <- as.character(pheno)
names(pheno) <- rownames(sce@meta.data)
VEC=sce@reductions$umap@cell.embeddings
# Extended Data Figure8
plotCytoTRACE(results, phenotype = pheno,emb = VEC)


# Supplementary Figure5 Ecm gene heatmap 
highlight <- read_excel('ECM_gene.xlsx') 
DefaultAssay(sce) <- 'RNA'
sce <- ScaleData(sce)
mt <- as.data.frame(t(as.matrix(GetAssayData(sce, assay = "RNA", slot = "scale.data"))))
group.by <- 'cell_type'
mt <- aggregate(mt, by=list(sce@meta.data[[group.by]]), FUN="mean")
rownames(mt) <- mt$Group.1
mt <- t(mt[,-1])
highlight <- highlight[order(-highlight[,'epi.cytotrce.cor']),]
cts_epi <- as.matrix(mt[highlight$Gene,])
pheatmap(cts_epi,show_colnames =T,show_rownames = T,color = viridis(8),
         cluster_rows = F,
         cluster_cols = F,
         name= 'Scaled Expression')

# Extended Data Figure8 Absolute time and relative time
sample <- subset(sce,downsample=2000)
data <- sample[,sample$day %in% c('D4','D7','D10','D13','D16')]@meta.data[,c('cell.type','day','cyto_pseudotime')]

data <- sce[,sce$day %in% c('D4','D7','D10','D13','D16')]@meta.data[,c('cell.type','day','cyto_pseudotime')]
ggplot(data,aes(x=cell.type,y=cyto_pseudotime))+
  geom_violin(width =0.8,fill='grey90',color='grey90')+
  geom_quasirandom(aes(color=day),width = 0.2,size=0.2)+ 
  stat_summary(aes(fill=day), geom="point",
               fun = mean, shape=21, size=2.5,stroke=0.5)+ 
  scale_color_manual(name = 'day',
                     values = c('#57C3F3',"#E5D2DD","#53A85F","#F1BB72",'#E63863'),
                     labels = c('D4','D7','D10','D13','D16'))+ 
  scale_fill_manual(values = c('#57C3F3',"#E5D2DD","#53A85F","#F1BB72",'#E63863'),
                    guide = 'none')+ 
  scale_y_continuous(limits = c(0,1))+
  theme_classic(base_size = 15)+
  coord_cartesian(clip = 'off')+ 
  geom_point(aes(x=2.5,y=6),pch=21,size=4,stroke=1.3)+
  theme(plot.margin = margin(10,50,10,10), 
        legend.position = c(1,0.3),
        axis.text = element_text(color='black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))+ 
  geom_point(aes(x=5.5,y=0.7),pch=21,size=2.5,stroke=0.5)+ 
  annotate(geom = 'text',label = 'Mean',x=5.8,y=0.7,size=4) +ylab('Predicted ordering')+xlab('Cell type') 

saveRDS(sce,file = 'epi.rds')


############# Pseudotime and pseudospace of epithelial cells 
library(monocle)
library(URD)
library(clusterProfiler)
library(org.Hs.eg.db)

###### monocle2 ---------------------------------------------------------------
load(file = 'epi.RData')
DefaultAssay(counts) <- 'RNA'
data <- as(as.matrix(counts@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = counts@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
disp_table <- dispersionTable(mycds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
mycds <- setOrderingFilter(mycds, unsup_clustering_genes$gene_id)
plot_ordering_genes(mycds)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)


#  supp_figure9 Epithelium DPT
library(reticulate)
library(viridisLite)

load(file = 'seu.rdata')
col <- c('#e7c23e',"#A6CEE3","#1F78B4","#FB9A99","#E95C59","#6A3D9A","#33A02C" )
DimPlot(seu, reduction = "diffmap", group.by = "celltype", cols = col) 

pseu <- read.csv(file = 'epiremovegland_pseudotime.csv')
seu$dpt_pseudotime <- pseu$dpt_pseudotime
tmp <- as.data.frame(seu@reductions$diffmap@cell.embeddings)
tmp <- data.frame(tmp,seu@meta.data)
# Figure4
ggplot(tmp, aes(x = DC_1, y = DC_2, colour = dpt_pseudotime)) + geom_point() + facet_wrap(tmp$day)+
        scale_color_gradientn(colors = viridis(8))
ggplot(tmp, aes(x = DC_1, y = DC_2, colour = dpt_pseudotime)) + geom_point(size=0.3)+
        scale_color_gradientn(colors = viridis(8))+theme_classic() 


# heatmap of dynamic genes along development trajectory  ---------
load(file = 'seu.rdata')

cds_sub <- mycds[,pData(mycds)$celltype %in% c('Precursor',"Antrum1",'Antrum2','Antrum3','Fundus1','Fundus2')] 
identical(colnames(seu),colnames(cds_sub))
pData(cds_sub)$Pseudotime <- seu$dpt_pseudotime
pseudotime_de <- differentialGeneTest(cds_sub[rownames(cds_sub),], fullModelFormulaStr = "~sm.ns(Pseudotime)")

Idents(seu) <- seu$celltype
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
de <- pseudotime_de[unique(markers$gene), ]
de <- subset(de, qval < 0.05)

fundus <- subset(seu, subset = celltype %in%  c("Precursor","Fundus1",'Fundus2') )
antrum <- subset(seu, subset = celltype %in%  c("Precursor","Antrum1",'Antrum2','Antrum3'))

fundus@meta.data <- fundus@meta.data %>% arrange(-dpt_pseudotime)
antrum@meta.data <- antrum@meta.data %>% arrange(dpt_pseudotime)

mat1 <- fundus@assays$RNA@data[de$gene_short_name,rownames(fundus@meta.data)]
mat2 <- antrum@assays$RNA@data[de$gene_short_name,rownames(antrum@meta.data)]
mat <- cbind(mat1,mat2)
mat <- t(apply(mat,1,function(x){smooth.spline(x,df=3)$y}))
mat <- t(apply(mat,1,function(x){(x-mean(x))/sd(x)}))

meta1 <- fundus@meta.data[,c("celltype","day")]
meta2 <- antrum@meta.data[,c("celltype","day")]
meta <- rbind(meta1,meta2)
meta$celltype <- factor(meta$celltype,levels = c("Precursor","Fundus1","Fundus2",'Antrum1','Antrum2','Antrum3'))

type_col <- setNames(col <- c("#FFD147FF","#A6CEE3","#1F78B4","#FB9A99","#E95C59","#6A3D9A" ),
                     c("Precursor","Fundus1","Fundus2",'Antrum1','Antrum2','Antrum3'))
col_anno <- HeatmapAnnotation(Cell_type=meta$celltype,
                              col=list(Cell_type = type_col))
htmp <- Heatmap(
  mat,
  name                         = "z-score",
  km = 5,
  col                          = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
  show_row_names               = F,
  show_column_names            = FALSE,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation = col_anno,
  column_split = c(rep("1_Fundus",17443),rep("2_Antrum",31913))) 
dht <- draw(htmp)

number <- c(as.numeric(dht@ht_list[["z-score"]]@row_order_list[["3"]]),#1164
            as.numeric(dht@ht_list[["z-score"]]@row_order_list[["4"]]),#552
            as.numeric(dht@ht_list[["z-score"]]@row_order_list[["5"]]),#307
            as.numeric(dht@ht_list[["z-score"]]@row_order_list[["1"]]),#334
            as.numeric(dht@ht_list[["z-score"]]@row_order_list[["2"]]))#273
order <- data.frame(number=number)
dt <- data.frame(gene=rownames(mat),number=1:2630)
dt <- full_join(order,dt,by = 'number')
matrix <- mat[dt$gene,]
# Extended Data Figure9 heatmap
Heatmap(
  matrix,
  name                         = "Smoothed expression",
  km = 0,
  col                          = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
  show_row_names               = F,
  show_column_names            = FALSE,
  cluster_rows                 = F,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation = col_anno,
  column_split = c(rep("1_Fundus",17443),rep("2_Antrum",31913)),
  row_split = c(rep("C1",1164),rep("C2",552),rep("C3",307),rep("C4",334),rep("C5",273))) 

# GO analysis
dt$cluster <- c(rep("C1",1164),rep("C2",552),rep("C3",307),rep("C4",334),rep("C5",273))
gcSample=split(dt$gene,dt$cluster)
bp <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     keyType = 'SYMBOL',
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1
)
dotplot(bp) + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))

# Absolute time and relative time
seu <- seu[,seu$day %in% c('D7','D10','D13','D16')]
seu$day <- ordered(seu$day,levels =c('D16', 'D13', 'D10', 'D7'))
Idents(seu) <- seu$celltype
sample <- subset(seu,downsample=500)
data <- sample@meta.data[,c('celltype1','day','dpt_pseudotime')]
data <- data %>%  
  arrange(factor(day, levels = c('D16', 'D13', 'D10', 'D7')))   

col <- c('#E63863',"#F1BB72","#53A85F","#E5D2DD")
# Extended Data Figure9
ggplot(data,aes(x=celltype1,y=dpt_pseudotime))+
  geom_violin(width =0.8,fill='grey90',color='grey90')+
  geom_quasirandom(aes(color=day),width = 0.2,size=0.01)+ 
  stat_summary(aes(fill=day), geom="point",
               fun = mean, shape=21, size=2.5,stroke=0.5)+ 
  scale_color_manual(name = 'day',
                     values = col,
                     labels = c('D16', 'D13', 'D10', 'D7'))+ 
  scale_fill_manual(values = col,
                    guide = 'none',
                    labels = c('D16', 'D13', 'D10', 'D7'))+ 
  scale_y_continuous(limits = c(0,1))+
  theme_classic(base_size = 15)+
  coord_cartesian(clip = 'off')+ 
  geom_point(aes(x=2.5,y=6),pch=21,size=4,stroke=1.3)+
  theme(plot.margin = margin(10,50,10,10), 
        legend.position = c(1,0.3),
        axis.text = element_text(color='black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))+ 
  geom_point(aes(x=6.8,y=0.7),pch=21,size=2.5,stroke=0.5)+ 
  annotate(geom = 'text',label = 'Mean',x=7.1,y=0.7,size=4) +
  xlab('Cell type')+RotatedAxis()


######## figure5 URD on day 16 epithelial cells ------------------------------
d16 <- epiremovegland[, epiremovegland$day %in% "D16" ] 
d16$cell.type <- ordered(d16$cell.type,levels = c("Fundus1","Fundus2","Antrum1","Antrum2","Antrum3"))

de.gland <- d16
count <- de.gland@assays$RNA@counts
meta <- de.gland@meta.data
object <- createURD(count.data=count, meta=meta, min.cells = 20, min.counts=20)
rm(list=c("count", "meta"))
shhhh <- gc()
stages <- unique(object@meta$cell.type)
cells.each.stage <- lapply(stages, function(stage) rownames(object@meta)[which(object@meta$cell.type == stage)])
var.genes.by.stage <- lapply(1:length(stages), function(n) findVariableGenes(object, cells.fit=cells.each.stage[[n]], set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, main.use=stages[[n]], do.plot=T))
names(var.genes.by.stage) <- stages
var.genes <- sort(unique(unlist(var.genes.by.stage))) 
object@var.genes <- var.genes

object <- calcPCA(object)
set.seed(18)
object <- calcTsne(object, perplexity = 30, theta=0.5)
set.seed(17)
object <- graphClustering(object, dim.use="pca", num.nn=c(15,20,30), do.jaccard=T, method="Louvain")
object <- calcKNN(object, nn=100)

object <- calcDM(object, knn=200, sigma.use=8)
root.cells <- rownames(object@meta)[object@meta$cell.type=="Fundus1"]
flood.result <- floodPseudotime(object, root.cells=root.cells, n=10, minimum.cells.flooded=2, verbose=T)
object <- floodPseudotimeProcess(object, flood.result, floods.name="pseudotime", max.frac.NA=0.4, pseudotime.fun=mean, stability.div=20)
gg.data <- cbind(object@pseudotime, object@meta[rownames(object@pseudotime),])
# Figure4
col <- c("#A6CEE3","#1F78B4","#FB9A99","#E95C59","#6A3D9A","#33A02C" )
ggplot(gg.data, aes(x=pseudotime, color=cell.type, fill=cell.type))+ geom_density(alpha=0.5)+ 
  theme_bw()+scale_fill_manual(values=col)+ scale_color_manual(values=col)

cds_sub <- mycds[,pData(mycds)$cell.type %in% c("Antral Epi","Fundic Epi")]
cds_sub <- cds_sub[,pData(cds_sub)$day %in% 'D16']
pse <- as.data.frame(object@pseudotime)
pse$cell <- rownames(pse)
pse[is.na(pse)] <- 0
pdata <- pData(cds_sub)
pdata$cell <- rownames(pdata)
mg <- full_join(pdata,pse,by = 'cell')
pData(cds_sub)$Pseudotime <- mg$pseudotime

pdata <- pData(cds_sub)
identical(colnames(d16),rownames(pdata))
d16$Pseudotime <- pdata$Pseudotime
d16@meta.data <- d16@meta.data %>% arrange(Pseudotime)
mat <- d16@assays$RNA@data[,rownames(d16@meta.data)]
mat <- t(apply(mat,1,function(x){smooth.spline(x,df=3)$y}))
mat <- t(apply(mat,1,function(x){(x-mean(x))/sd(x)}))
meta <- d16@meta.data[,c("cell.type",'day')]
meta$cell.type <- factor(meta$cell.type,levels = c("Fundus1","Fundus2",'Antrum1','Antrum2','Antrum3'))
type_col <- setNames(c("#A6CEE3","#1F78B4","#FB9A99","#E95C59","#6A3D9A"),
                     c("Fundus1","Fundus2",'Antrum1','Antrum2','Antrum3'))
col_anno <- HeatmapAnnotation(Cell_type=meta$cell.type,col=list(Cell_type = type_col))

# Extended Data Figure9 HOX gene heatmap
hox <- c('HOXA1',"HOXC5",'HOXB6',"HOXA2","HOXB3","HOXB4","HOXB2","HOXB7") 
matrix <- mat[hox,]
Heatmap(
  matrix,
  name                         = "Smoothed expression",
  km = 0,
  col                          = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
  show_row_names               = T,
  show_column_names            = FALSE,
  cluster_rows                 = F,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation = col_anno) 

# Extended Data Figure9 ECM gene heatmap
highlight <- read_excel('ECM_gene.xlsx') 
highlight <- as.data.frame(highlight)
highlight <- highlight[order(-highlight[,'epi.cytotrce.cor']),]

pseudotime_de <- differentialGeneTest(cds_sub[highlight$Gene,], fullModelFormulaStr = "~sm.ns(Pseudotime)") 
pseudotime_de <- subset(pseudotime_de, qval < 1e-2) 
ecm <- pseudotime_de$gene_short_name
matrix <- mat[ecm,]
Heatmap(
  matrix,
  name                         = "Smoothed expression",
  km = 0,
  col                          = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
  show_row_names               = T,
  show_column_names            = FALSE,
  cluster_rows                 = F,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation = col_anno) 

# Extended Data Figure9 TF heatmap
tf.plot <- read.table(file = 'plotgene.txt')
matrix <- mat[tf.plot,]
Heatmap(
  matrix,
  name                         = "Smoothed expression",
  km = 0,
  col                          = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
  show_row_names               = T,
  show_column_names            = FALSE,
  cluster_rows                 = F,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation = col_anno,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete")









