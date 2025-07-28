## Integrate all samples 
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(ggsci)
library(ComplexHeatmap)

load(file = "D4.O.S1.1.filter.rdata") 
load(file = "D4.O.S1.2.filter.rdata") 
load(file = "D4.O.S2.1.filter.rdata") 
load(file = "D4.O.S2.2.filter.rdata") 

load(file = "D7.O.S1.1.filter.rdata") 
load(file = "D7.O.S1.2.filter.rdata") 
load(file = "D7.O.S2.1.filter.rdata") 
load(file = "D7.O.S2.2.filter.rdata") 


load(file = "D10.O.S1.1.filter.rdata") 
load(file = "D10.O.S2.1.filter.rdata") 
load(file = "D10.O.S2.2.filter.rdata") 

load(file = "D13.O.S1.1.filter.rdata") 
load(file = "D13.O.S1.2.filter.rdata") 
load(file = "D13.O.S2.1.filter.rdata") 
load(file = "D13.O.S2.2.filter.rdata") 

load(file = 'PDMS_1.filter.rdata') 
load(file = 'PDMS_2.filter.rdata') 
load(file = 'dish.filter.rdata') 


sample.list <- list(D4.O.S1.1.f,D4.O.S1.2.f,D4.O.S2.1.f,D4.O.S2.2.f,
                    D7.O.S1.1.f,D7.O.S1.2.f,D7.O.S2.1.f,D7.O.S2.2.f,D10.O.S1.1.f,D10.O.S2.1.f,
                    D10.O.S2.2.f,D13.O.S1.1.f,D13.O.S1.2.f,D13.O.S2.1.f,D13.O.S2.2.f,PDMS_1.f,PDMS_2.f,dish.f)
for (i in 1:length(sample.list)){
  sample.list[[i]] <- NormalizeData(sample.list[[i]], verbose = FALSE)
  sample.list[[i]] <- FindVariableFeatures(sample.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

sample.anchors <- FindIntegrationAnchors(object.list = sample.list,k.filter=200, anchor.features = 2000, verbose = F)
sample.integrated <- IntegrateData(anchorset = sample.anchors, verbose = F) 

DefaultAssay(sample.integrated) <- "integrated"
sample.integrated <- ScaleData(sample.integrated, verbose = F)
sample.integrated <- RunPCA(sample.integrated,verbose = F)
ElbowPlot(sample.integrated, ndims = 50)

sample.integrated <- RunUMAP(sample.integrated,reduction = "pca", dims = 1:20, verbose = F)
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:20, verbose = F)
sample.integrated <- FindClusters(sample.integrated, resolution = seq(0.1,1,0.1), verbose = F)

save(sample.integrated,file = "./integrated.RData")

######## Extended Figure 6 
col=pal_igv('default',alpha = 1)(51)
p1 <- DimPlot(sample.integrated, reduction = "umap", group.by = "day",label = F,raster=FALSE,cols = col[35:51])
p2 <- DimPlot(sample.integrated, reduction = "umap", group.by = "cell.type",label = F,raster=FALSE,cols = col)

DimPlot(sample.integrated, reduction = "umap", group.by = "cell.type",split.by = 'day',label = F,raster=FALSE,cols = col)

Cellratio <- prop.table(table(sample.integrated$cell.type, sample.integrated$day), margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio$cell_type <- Cellratio$Var1
Cellratio$day <- Cellratio$Var2
ggplot(Cellratio) + 
  geom_bar(aes(x =day, y= Freq, fill = cell_type),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Day',y = 'Cell type')+
  scale_fill_manual(values = col)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

# marker genes dotplot
DefaultAssay(sample.integrated) <- "RNA"
gene <- c('POU5F1','NANOG', #hPSC
          'SOX17','GSC',#DE
          'CLDN4','CLDN18','GATA4','CDH1',#Epithelium
          'COL1A2','COL3A1',#Mesenchymal
          'PCGF6',#NE
          'SOX2','PAX6','NES','OLIG2',#NPC
          'SOX10','PAX3','EDNRB',#Premigratory ENNC
          'CDH6','SOX11','B3GAT1',#migratory ENNC
          'TUBB3','MAP2','NGFR',#Neuron
          'CDH5','FLT1', #Endothelial
          'SST','GHRL')#Enteroendocrine
DotPlot(sample.integrated, features =gene, group.by = "cell.type",cols = c("lightgrey",'#FF0000'))+ RotatedAxis() 

FeaturePlot(sample.integrated,features = c('EPCAM','PAX6','COL3A1'),cols = c('lightgrey','red'),ncol = 3,label = F)

######## epithelial subtypes ------------------------------- 
dt <- subset(sample.integrated, cell.type == "Epithelium")
counts <- dt@assays$RNA@counts
counts <- CreateSeuratObject(counts = counts)
counts$orig.ident <- dt$orig.ident
counts$day <- dt$day
list <- SplitObject(counts, split.by = "orig.ident")
for (i in 1:length(list)){
  list[[i]] <- NormalizeData(list[[i]], verbose = FALSE)
  list[[i]] <- FindVariableFeatures(list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = 2000)
counts <- IntegrateData(anchorset = anchors) #26186
DefaultAssay(counts) <- "integrated"
counts <- counts %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:20)
counts <- FindNeighbors(counts, reduction = "pca", dims = 1:20, verbose = F)
counts <- FindClusters(counts,resolution = 0.1, verbose = F)
counts <- FindClusters(counts,resolution = 0.4, verbose = F)
counts <- FindClusters(counts,resolution = 0.6, verbose = F)
counts <- FindClusters(counts,resolution = 1.2, verbose = F)

counts$cell.type <- plyr::mapvalues(x=counts$cell.type,from=c('0','1','2','3','4',
                                                              '5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21'),
                                    to=c('Antrum1','Antrum1','Precursor','Antrum1','Antrum1',
                                         'Fundus1','Antrum1','Fundus2','Fundus1','Antrum1','Antrum2',
                                         'Fundus1','Gland','Antrum1','Precursor','Gland','Antrum1','Fundus1','Fundus1','Antrum2','Precursor','Antrum2'))
counts$cell.type <- ordered(counts$cell.type,levels=c('Precursor','Fundus1','Fundus2','Antrum1','Antrum2','Antrum3','Gland'))


# Supplementary Figure 5
col <- c('#e7c23e',"#A6CEE3","#1F78B4","#FB9A99","#E95C59","#6A3D9A","#33A02C" )
DimPlot(counts,group.by = 'cell.type',cols =col)

# Extended Figure9 featureplot 
FeaturePlot(counts,features = c('PDX1','SOX2'),cols = c('lightgrey','red'),split.by = 'day',label = F)

# Supplementary Figure 5 DEG heatmap
epiremovedgland <- subset(counts, cell.type %in% c("Fundus1", "Fundus2", "Antrum1", "Antrum2", "Antrum3")) 

DefaultAssay(epiremovedgland) <- 'RNA'
markers <- FindAllMarkers(epiremovedgland, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- filter(markers,p_val_adj < 0.05)
markers %>%
  group_by(cluster) %>%
  top_n(n = 10 ,wt = avg_log2FC) -> top10

mt <- as.data.frame(t(as.matrix(GetAssayData(epiremovedgland, assay = "RNA", slot = "scale.data"))))
group.by <- 'cell.type'
mt <- aggregate(mt, by=list(epiremovedgland@meta.data[[group.by]]), FUN="mean")
rownames(mt) <- mt$Group.1
mt <- t(mt[,-1])

cts <- as.matrix(mt[unique(markers$gene),])
bk <- c(seq(-1,0,by=0.005),seq(0.001,0.4,by=0.001))
ComplexHeatmap::pheatmap(cts,show_colnames =T,show_rownames = F,breaks = bk,legend_breaks=seq(-1,0.4,0.2),
                         color =colorRampPalette(rev(brewer.pal(n = 35, name ="RdYlBu")))(length(bk)), 
                         cluster_rows = F,
                         cluster_cols = F,
                         name= 'Scaled Expression')


######## Supplementary Figure 3 Partial Epi and Epithelium DEG heatmap -------------------
subset <- subset(sample.integrated, cell.type %in% c('Partial Epi', 'Epithelium'))
subset$cell.type <- ordered(subset$cell.type,levels=c('Partial Epi','Epithelium'))

DefaultAssay(subset) <- 'RNA'
subset <- ScaleData(subset,features = rownames(subset))
Idents(subset) <- subset$cell.type
markers <- FindAllMarkers(subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- filter(markers,p_val_adj < 0.05)

mt <- as.data.frame(t(as.matrix(GetAssayData(subset, assay = "RNA", slot = "scale.data"))))
group.by <- 'cell.type'
mt <- aggregate(mt, by=list(subset@meta.data[[group.by]]), FUN="mean")
rownames(mt) <- mt$Group.1
mt <- t(mt[,-1])
cts <- as.matrix(mt[unique(markers$gene),])
cts <- na.omit(cts)
bk <- c(seq(-0.5,0,by=0.01),seq(0.005,0.2,by=0.005))
ComplexHeatmap::pheatmap(cts,show_colnames =T,show_rownames = F,breaks = bk,legend_breaks=seq(-0.5,0.2,0.1),
                         color =colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(length(bk)), 
                         cluster_rows = F,
                         cluster_cols = F,
                         name= 'Scaled Expression')

gene <- c('POU5F1','SOX17','CXCR4','SFRP1','CDH1','FOS','CDH6')
DotPlot(subset,group.by = 'cell.type',features = rev(gene),cols = c('lightgray','red'))+RotatedAxis()+coord_flip()


######## Mesenchymal subtypes ------------------------------- 
mesen <- subset(sample.integrated, cell.type == "Mesenchymal")
DefaultAssay(mesen) <- 'integrated'
mesen <- RunPCA(mesen) 
mesen <- RunUMAP(mesen, reduction = "pca", dims = 1:20)
mesen <- FindNeighbors(mesen, dims = 1:10)
mesen <- FindClusters(mesen, resolution = seq(0.1,1,0.1) )
DimPlot(mesen, reduction = "umap", group.by = "integrated_snn_res.0.1",raster=FALSE,label = TRUE)
mesen$cell.type <- mesen$integrated_snn_res.0.1
mesen$cell.type <- plyr::mapvalues(x=mesen$cell.type, from=c("0","2","1"), 
                                   to=c("Mesenchymal1","Mesenchymal2","Mesenchymal3"))
# Supplementary Figure 4 DEG heatmap
DefaultAssay(mesen) <- 'RNA'
marker <- FindAllMarkers(mesen, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
marker <- filter(marker,p_val_adj < 0.05)
marker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

mt <- as.data.frame(t(as.matrix(GetAssayData(mesen, assay = "RNA", slot = "scale.data"))))
group.by <- 'cell.type'
mt <- aggregate(mt, by=list(mesen@meta.data[[group.by]]), FUN="mean")
rownames(mt) <- mt$Group.1
mt <- t(mt[,-1])
cts <- as.matrix(mt[unique(marker$gene),])
cts <- na.omit(cts)
bk <- c(seq(-0.6,0,by=0.01),seq(0.01,0.5,by=0.01))
ComplexHeatmap::pheatmap(cts,show_colnames =T,show_rownames = F,breaks = bk,legend_breaks=seq(-1,0.5,0.5),
                         color =colorRampPalette(rev(brewer.pal(n = 35, name ="RdYlBu")))(length(bk)), 
                         cluster_rows = F,
                         cluster_cols = F,
                         name= 'Scaled Expression')


save(mesen,file = 'mesen.rdata')

######## Neuronal subtypes ------------------------------- 
neuron <- subset(sample.integrated, Major_cell_type == 'Neuronal') 

DefaultAssay(neuron) <- 'integrated'
neuron <- RunPCA(neuron) 
neuron <- RunUMAP(neuron, reduction = "pca", dims = 1:20)

# Neuron1 and Neuron2
Neu <- subset(neuron,cell.type=='Neuron')
Neu <- RunPCA(Neu, npcs = 30, verbose = T)
Neu <- RunUMAP(Neu, seed.use = -1, reduction = "pca", dims = 1:20,return.model = TRUE)
Neu <- FindNeighbors(Neu, reduction = "pca", dims = 1:20)
Neu  <- FindClusters(Neu,resolution = 0.025)
unique(Neu$integrated_snn_res.0.025)
DimPlot(Neu, reduction = "umap", group.by = "integrated_snn_res.0.025",label = T,repel = T,raster=FALSE)

Neu1 <- subset(Neu,integrated_snn_res.0.025==0)
num1 <- match(colnames(Neu1),colnames(neuron))
neuron$cell.type <- as.character(neuron$cell.type)
neuron$cell.type[num1] <- 'Neuron1'
Neu2 <- subset(Neu,integrated_snn_res.0.025==1)
num2 <- match(colnames(Neu2),colnames(neuron))
neuron$cell.type[num2] <- 'Neuron2'
DimPlot(neuron,group.by = 'cell.type',label = F,cols = c("#A6CEE3","#1F78B4","#e7c23e","#802268FF", "#6BD76BFF","#FF4500"))

# Extended Figure 7 heatmap 
mt.lab <- as.data.frame(t(as.matrix(GetAssayData(neuron, assay = "RNA", slot = "scale.data"))))
group.by <- 'cell.type'
mt.lab <- aggregate(mt.lab, by=list(neuron@meta.data[[group.by]]), FUN="mean")
rownames(mt.lab) <- mt.lab$Group.1
mt.lab <- t(mt.lab[,-1])

cts <- as.matrix(mt.lab[rev(c('GBX2','SOX3','NKX6-1','HOXA4','HOXA5','ERBB3','SOX10','TBX3')),])
bk <- c(seq(-1,1,by=0.1),seq(1,2,by=0.1))
pheatmap(cts, cluster_cols = F,cluster_rows = F,breaks = bk,legend_breaks=seq(-1,2,1),
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(length(bk))
)



