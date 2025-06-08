######## compare with in vivo data  
library(Matrix)
library(Seurat)
library(ggsci)
library(pheatmap)
library(RColorBrewer)

# Article:Charting human development using a multiendodermal organ atlas and organoid models
invivo <- readRDS(file = 'Dat_HIO_Cl12_and_d132_stem_cell_fetal_pt_gene_average_expr_as_ref_for_qp.rds')
index <- read.csv(file = 'Table_fetal_atlas_cell_index.csv',header = F)
symbol <- read.csv(file = 'Table_fetal_atlas_gene_symbol.csv',header = F)
info <- read.csv(file = 'Table_fetal_atlas_meta_info.csv')
count <- Matrix::readMM('Table_fetal_atlas_count.mtx')
load(file = 'integrated.RData')

rownames(count) <- symbol$V1 
colnames(count) <- index$V1
all <- CreateSeuratObject(count,meta.data = info) 
all <- all %>%
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(dims = 1:20)%>%
  RunUMAP(dims = 1:20)
cell.embeddings <- info[,c('UMAP_X','UMAP_Y')]
cell.embeddings <- as.matrix(cell.embeddings)
colnames(cell.embeddings)[1] <- 'UMAP_1'
colnames(cell.embeddings)[2] <- 'UMAP_2'
all@reductions$umap@cell.embeddings <- cell.embeddings
col <- levels(as.factor(all$cluster_col))
DimPlot(all, reduction = "umap", group.by = "Corrected_tissue",raster=FALSE,split.by = 'Corrected_tissue',ncol = 4,label = F,cols = col) 

stomach <- all[, all$Corrected_organ_group %in% 'Stomach']
p1 <- DimPlot(stomach, reduction = "umap", group.by = "Age_week",raster=FALSE,label = F,cols = col[15:27])
p2 <- DimPlot(stomach, reduction = "umap", group.by = "Major_cell_type",raster=FALSE,label = TRUE,cols = col) 
p3 <- DimPlot(stomach, reduction = "umap", group.by = "Cell_type",raster=FALSE,label = F,cols = col)
p1|p2|p3

# Figure3 Dotplot
load(file = 'integrated.RData')

d16 <- sample.integrated[, sample.integrated$day %in% 'D16'] 
d16$Corrected_organ_group <- d16$day
all.sample <- all[, all$Corrected_tissue %in% c('Esophagus','Duodenum','Stomach','Stomach-antrum','Stomach-corpus')] 
all.sample_d16 <- merge(all.sample,d16)
all.sample_d16 <- NormalizeData(all.sample_d16)
all.sample_d16 <- ScaleData(all.sample_d16, features = rownames(all.sample_d16))

all.sample_d16_epi <- all.sample_d16[, all.sample_d16$Major_cell_type %in% c('Epi','Epithelial')] 
all.sample_d16_epi$Corrected_organ_group <- ordered(all.sample_d16_epi$Corrected_organ_group,levels=c('Intestine','Esophagus','Stomach','d16'))
DefaultAssay(all.sample_d16_epi) <- 'RNA'
gene <- c('CLDN18','SOX2','CDH1', # Stomach epithelium
          'TP63','KRT5', # Esophagus epithelium
          'CDX2') # Intestine epithelium
DotPlot(all.sample_d16_epi,features = gene,cols = c('lightgrey','red'),group.by = 'Corrected_organ_group')+RotatedAxis()

# in vivo/in vitro UMAP and correlation heatmap
set.seed(20230518)

sample.integrated$Major_cell_type <- result$celltype
Idents(sample.integrated) <- 'Major_cell_type'
sample<-subset(sample.integrated,downsample=2000)#23068
vivo$celltype <- paste0('lite_ ',vivo$Major_cell_type)

vivo$sample <- vivo$orig.ident
vivo$day <- 'Vivo'
da1 <- data.frame(colnames(sample),sample$day)
da2 <- data.frame(unique(sample$day),c('D4','D7','Other','Other','Other'))
colnames(da2) <- c('sample.day','day')
result <- join(da1,da2)
sample$sample <- result$day

all_sample.combined <- merge(vivo,sample)
alldata.list <- SplitObject(all_sample.combined, split.by = "sample")
alldata.list <- lapply(X = alldata.list, FUN = function(x) {
  DefaultAssay(x)="RNA"  
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 2000)
})

all_sample.anchors <- FindIntegrationAnchors(object.list = alldata.list , dims = 1:20)
all_sample.combined <- IntegrateData(anchorset = all_sample.anchors, dims = 1:20)

DefaultAssay(all_sample.combined) <- "integrated"
all_sample.combined <- ScaleData(all_sample.combined,verbose = T)
all_sample.combined <- RunPCA(all_sample.combined, npcs = 30, verbose = T)
ElbowPlot(all_sample.combined)
all_sample.combined <- RunUMAP(all_sample.combined, seed.use = -1, reduction = "pca", dims = 1:20,return.model = TRUE)
all_sample.combined <- FindNeighbors(all_sample.combined, reduction = "pca", dims = 1:20)

d <- data.frame(table(all_sample.combined$orig.ident))
all_sample.combined <- FindClusters(all_sample.combined,  resolution = c(2))
P <- DimPlot(all_sample.combined, reduction = "umap", group.by = "integrated_snn_res.2",label = T,repel = T,raster=FALSE)

all_sample.combined$Major_cell_type[num] <- paste0('Fetal_stomach_',vivo$Major_cell_type)
cols <- c('#ffff00','#1ce6ff','#ff34ff','#ff4a46','#008941','#006fa6','#a30059',
          '#ffdbe5','#7a4900','#0000a6','#63ffac','#b79762',"#ff4a46","#008941","#ffdbe5","#0000a6","#eec3ff","#456d75")
all_sample.combined$Major_cell_type<- factor(all_sample.combined$Major_cell_type,levels = c("hPSC","DE","Partial Epi","Gastirc Epi","Mesenchymal","NE",'NPC',
                                                                                            'Neuron','ENCC','Endothelial','Enteroendocrine','Unidentified',
                                                                                            'Fetal_stomach_Epithelial','Fetal_stomach_Mesenchymal','Fetal_stomach_Neuronal',
                                                                                            "Fetal_stomach_Endothelial",'Fetal_stomach_Erythroid',"Fetal_stomach_Immune"))
all_sample.combined$day<- factor(all_sample.combined$day,levels = c("D4","D7","D10","D13","D16","Vivo"))

# Extended Data Figure 6 
DimPlot(all_sample.combined, group.by = "Major_cell_type",split.by = 'day',label = F,repel = T,raster=FALSE,cols = cols)+
  facet_wrap(~ day, nrow = 2)

data <- subset(all_sample.combined,day=='Vivo'|day=='D16')
num <- na.omit(match(colnames(vivo),colnames(data)))
data$Major_cell_type[num] <- paste0('Fetal_stomach_',vivo$Major_cell_type)
data <- subset(data,Major_cell_type!='Enteroendocrine'&Major_cell_type!="Partial Epi"&Major_cell_type!="hPSC"&Major_cell_type!="DE"&Major_cell_type!="Erythroid"&Major_cell_type!="Immune"&Major_cell_type!="Unidentified"&Major_cell_type!="Fetal_stomach_Erythroid"&Major_cell_type!="Fetal_stomach_Immune")

DEG <- FindAllMarkers(vivo,group.by = 'Major_cell_type')
DEG <- dplyr::filter(DEG,DEG$p_val_adj<0.05&DEG$avg_log2FC>0)

alldata.list <- SplitObject(data, split.by = "orig.ident")
alldata.list <- lapply(X = alldata.list, FUN = function(x) {
  DefaultAssay(x)="RNA"  
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 2000)
})

all_sample.anchors <- FindIntegrationAnchors(object.list = alldata.list , dims = 1:20,anchor.features = 5000)
all_sample.combined <- IntegrateData(anchorset = all_sample.anchors,features = unique(DEG$gene))

DefaultAssay(all_sample.combined) <- 'integrated'
all_sample.combined <- FindVariableFeatures(all_sample.combined,nfeatures = 10000)
all_sample.combined <- ScaleData(all_sample.combined)
aver <- AverageExpression(all_sample.combined,group.by = 'Major_cell_type',slot = 'scale.data')
exp <- data.frame(aver$integrated)

num=names(tail(sort(apply(exp, 1, sd)),10000))
num <- match(num,rownames(exp))
exp <- exp[num,]
exp <- data.frame(exp)
exp <- cor(exp, method= "spearman")

exprTable_t <- as.data.frame(t(exp))
col_dist = dist(exprTable_t)
hclust_1 <- hclust(col_dist)
manual_order = c("Fetal.stomach.Mesenchymal",'Mesenchymal','Fetal.stomach.Epithelial','Gastirc.Epi','Fetal.stomach.Endothelial','Endothelial',
                 'Fetal.stomach.Neuronal','Neuron','ENCC','NPC','NE')
dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable_t))), agglo.FUN = max)
col_cluster <- as.hclust(dend)
bk <- c(seq(0,0.2,by=0.05),seq(0.21,0.3,by=0.03),seq(0.31,0.4,by=0.03),seq(0.41,0.64,by=0.23),seq(0.64,0.82,by=0.18),seq(0.82,0.98,by=0.16),seq(0.98,1,by=0.02))###0.02(有颜色)或者0.03(白)
# Figure3
pheatmap(exp,cluster_cols = col_cluster,cluster_rows = col_cluster,display_numbers = F,fontsize = 20,angle_col = '315') 

# Figure3 D16+vivo integrated 
data <- subset(all_sample.combined,day=='Vivo'|day=='D16')
num <- na.omit(match(colnames(vivo),colnames(data)))

data$Sample[num] <- 'vivo'
data$Sample[-num] <- 'vitro'

S32 <- subset(data,orig.ident=='S32')
num2 <- match(colnames(S32),colnames(data))
data$Sample[num2] <- 'S32'

alldata.list <- SplitObject(data, split.by = "Sample")
alldata.list <- lapply(X = alldata.list, FUN = function(x) {
  DefaultAssay(x)="RNA"  
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 2000)
})

all_sample.anchors <- FindIntegrationAnchors(object.list = alldata.list , dims = 1:20)
all_sample.combined <- IntegrateData(anchorset = all_sample.anchors, dims = 1:20)

DefaultAssay(all_sample.combined) <- 'integrated'
all_sample.combined <- ScaleData(all_sample.combined)
all_sample.combined <- RunPCA(all_sample.combined, npcs = 20, verbose = T)
ElbowPlot(all_sample.combined)
all_sample.combined <- RunUMAP(all_sample.combined, seed.use = -1, reduction = "pca", dims = 1:20,return.model = TRUE)
all_sample.combined <- FindNeighbors(all_sample.combined, reduction = "pca", dims = 1:20)

num <- match(colnames(vivo),colnames(all_sample.combined))
all_sample.combined$sample[num] <- 'vivo'
all_sample.combined$sample[-num] <- 'vitro'

cols <- c('#ffff00','#1ce6ff','#ff34ff','#ff4a46','#008941','#006fa6','#a30059',
          '#ffdbe5','#7a4900','#0000a6','#63ffac','#b79762',"#ff4a46","#008941","#ffdbe5","#0000a6","#eec3ff","#456d75")
all_sample.combined$Major_cell_type<- factor(all_sample.combined$Major_cell_type,levels = c("hPSC","DE","Partial Epi","Gastirc Epi","Mesenchymal","NE",'NPC',
                                                                                            'Neuron','ENCC','Endothelial','Enteroendocrine','Unidentified',
                                                                                            'Fetal_stomach_Epithelial','Fetal_stomach_Mesenchymal','Fetal_stomach_Neuronal',
                                                                                            "Fetal_stomach_Endothelial",'Fetal_stomach_Erythroid',"Fetal_stomach_Immune"))

DimPlot(all_sample.combined, reduction = "umap", group.by = 'Major_cell_type',split.by = 'sample',label = T,repel = T,raster=FALSE,pt.size = 1.5,cols = cols)

# Extended Data Figure 6 
DefaultAssay(all_sample.combined) <- 'RNA'
vivo <- subset(all_sample.combined,Sample=='vivo')
vitro <-  subset(all_sample.combined,Sample=='vitro')

FeaturePlot(vivo,features = c('CDH1'),cols = c('lightgrey','red'),pt.size = 1)
FeaturePlot(vitro,features = c('CDH1'),cols = c('lightgrey','red'),pt.size = 1)
FeaturePlot(vivo,features = c('COL1A2'),cols = c('lightgrey','red'),pt.size = 1)
FeaturePlot(vitro,features = c('COL1A2'),cols = c('lightgrey','red'),pt.size = 1)
FeaturePlot(vivo,features = c('FLT1'),cols = c('lightgrey','red'),pt.size = 1)
FeaturePlot(vitro,features = c('FLT1'),cols = c('lightgrey','red'),pt.size = 1)
FeaturePlot(vivo,features = c('MAP2'),cols = c('lightgrey','red'),pt.size = 0.1)
FeaturePlot(vitro,features = c('MAP2'),cols = c('lightgrey','red'),pt.size = 0.5)



#### In vivo and In vitro epithelium 
# Article:Charting human development using a multiendodermal organ atlas and organoid models
index <- read.csv(file = './Table_fetal_stomach_epithelium_cell_index.csv',header = F)
symbol <- read.csv(file = './Table_fetal_stomach_epithelium_gene_symbol.csv',header = F)
info <- read.csv(file = './Table_fetal_stomach_epithelium_meta_info.csv')
count <- Matrix::readMM('./Table_fetal_stomach_epithelium_count.mtx')

rownames(count) <- symbol$V1 
colnames(count) <- index$V1
epi <- CreateSeuratObject(count,meta.data = info) 
epi <- epi %>%
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(dims = 1:20)%>%
  RunUMAP(dims = 1:20)
cell.embeddings <- info[,c('UMAP_X','UMAP_Y')]
cell.embeddings <- as.matrix(cell.embeddings)
colnames(cell.embeddings) <- c('umap_1','umap_2')
epi@reductions$umap@cell.embeddings <- cell.embeddings
DimPlot(epi, reduction = "umap", group.by = "Cell_type",
              raster=FALSE,label = T,cols = levels(as.factor(epi$Cell_type_color))) 

epi$subtype <- epi$Cell_type 
epi$subtype <- plyr::mapvalues(x=epi$subtype,from=c('Proliferative gastric epithelial precursor',
                                                    'Antrum stem cell','Antrum proliferative mucous cell','Corpus stem cell',
                                                    'Corpus mucous neck cell','Antral gland cell','Antral surface mucous cell','Corpus surface mucous cell'),
                               to=c('Fetal_stomach_precursor','Fetal_stomach_antrum','Fetal_stomach_antrum','Fetal_stomach_fundus','Fetal_stomach_fundus',
                                    'Fetal_stomach_gland','Fetal_stomache_gland','Fetal_stomach_gland'))
epi.subset <- epi[, epi$subtype %in% c('Fetal_stomach_precursor','Fetal_stomach_antrum','Fetal_stomach_fundus','Fetal_stomach_gland')] #784
epi.subset$split <- 'In vivo' 

# Figure4 Dotplot
load(file = 'allepi.reintegrate.rdata')

gene <- c('EZH2','FOXO3','TCF4', # Precursor
          'SOX2','IRX2','IRX3','IRX5','NR2F2', # Fundus
          'PDX1','MEIS2',# Antrum
          'TFF2','UPK1B', 'PLAC8'# Secreting
) 
DefaultAssay(counts) <- 'RNA'

p1 <- DotPlot(counts,features = gene,group.by = 'lab')+RotatedAxis()+coord_flip()+ 
  scale_color_gradientn(colors = c("grey85", brewer.pal(9, "OrRd"))) 
p2 <- DotPlot(epi.subset,features = gene,group.by = 'subtype')+RotatedAxis()+coord_flip()+ 
  scale_color_gradientn(colors = c("grey85", brewer.pal(9, "OrRd")))
p1|p2

# Figure4 FeaturePlot
d16 <- subset(counts,day=='D16')
# In vivo
FeaturePlot(epi.subset,features = 'PDX1',slot = 'scale.data',pt.size = 0.2,order = T)+
  scale_colour_gradientn(colours = colorRampPalette(c("#DCDCDC",'#FDF5E6','#FFA07A',"#fc2c14","#fc2c14",'#CD0000','#CD0000'))(30),labels = c("low", "high"),
                         breaks=c(-0.56,4.6))
FeaturePlot(epi.subset,features = 'LAMC2',slot = 'scale.data',pt.size = 0.2,order = T)+
  scale_colour_gradientn(colours = colorRampPalette(c("#DCDCDC",'#FDF5E6',"#fc2c14",'#CD0000','#CD0000','#CD0000','#CD0000'))(30),labels = c("low", "high"),
                         breaks=c(-0.36,6.6))
FeaturePlot(epi.subset,features = 'HOXB2',slot = 'scale.data',pt.size = 0.2,order = T)+
  scale_colour_gradientn(colours = colorRampPalette(c("#DCDCDC","#fc2c14","#fc2c14",'#CD0000','#CD0000','#CD0000','#CD0000'))(30),labels = c("low", "high"),
                         breaks=c(-0.18,10))
FeaturePlot(epi.subset,features = 'HOXC5',slot = 'scale.data',pt.size = 0.2)+
  scale_colour_gradientn(colours = colorRampPalette(c("#DCDCDC","#fc2c14",'#CD0000','#CD0000','#CD0000','#CD0000','#CD0000','#CD0000','#CD0000','#CD0000','#CD0000','#CD0000'))(30),
                         labels = c("low", "high"),breaks=c(-0.1,10))
FeaturePlot(epi.subset,features = 'NR2F2',slot = 'scale.data',pt.size = 0.2)+
  scale_colour_gradientn(colours = colorRampPalette(c("#DCDCDC",'#FDF5E6','#FFA07A',"#fc2c14",'#CD0000'))(30),labels = c("low", "high"),
                         breaks=c(-0.95,3))
# In vitro
FeaturePlot(d16,features = 'PDX1',slot = 'scale.data',pt.size = 0.2,order = T)+
  scale_colour_gradientn(colours = colorRampPalette(c("#DCDCDC",'#FDF5E6','#FFA07A',"#fc2c14",'#CD0000'))(30),labels = c("low", "high"),
                         breaks=c(-0.43,8.4))
FeaturePlot(d16,features = 'LAMC2',slot = 'scale.data',pt.size = 0.2,order = T)+
  scale_colour_gradientn(colours = colorRampPalette(c("#DCDCDC",'#FDF5E6','#FFA07A',"#fc2c14",'#CD0000','#CD0000'))(30),labels = c("low", "high"),
                         breaks=c(-0.37,10))
FeaturePlot(d16,features = 'HOXB2',slot = 'scale.data',pt.size = 0.2,order = T)+
  scale_colour_gradientn(colours = colorRampPalette(c("#DCDCDC",'#FDF5E6','#FFA07A',"#fc2c14",'#CD0000'))(30),labels = c("low", "high"),
                         breaks=c(-0.55,5.8))
FeaturePlot(d16,features = 'HOXC5',slot = 'scale.data',pt.size = 0.2)+
  scale_colour_gradientn(colours = colorRampPalette(c("#DCDCDC",'#FDF5E6','#FFA07A',"#fc2c14",'#CD0000','#CD0000','#CD0000','#CD0000','#CD0000'))(30),
                         labels = c("low", "high"),breaks=c(-0,10))
FeaturePlot(d16,features = 'NR2F2',slot = 'scale.data',pt.size = 0.2)+
  scale_colour_gradientn(colours = colorRampPalette(c("#DCDCDC","#DCDCDC",'#FDF5E6','#FFA07A',"#fc2c14",'#CD0000','#CD0000'))(30),labels = c("low", "high"),
                         breaks=c(-0.73,4.3))

# Figure4 integrate
epi.subset$orig.ident <- 'lite'
counts$subtype <- counts$lab
all.data <- merge(counts,lite_data)

alldata.list <- SplitObject(all.data, split.by = "orig.ident")
alldata.list <- lapply(X = alldata.list, FUN = function(x) {
  DefaultAssay(x)="RNA"  
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000)
})

DEG <- SetIdent(all.data,value = 'orig.ident')
DefaultAssay(DEG) <- 'RNA'
DEG <- FindAllMarkers(DEG,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5)
DEG <- filter(DEG,DEG$p_val_adj<0.05)
gene <- unique(DEG$gene)

all_sample.anchors <- FindIntegrationAnchors(object.list = alldata.list,dims =1:20,anchor.features =  gene)
all_sample.combined <- IntegrateData(anchorset = all_sample.anchors, dims = 1:20)
DefaultAssay(all_sample.combined) <- "integrated"

all_sample.combined <- ScaleData(all_sample.combined,verbose = T)
all_sample.combined <- RunPCA(all_sample.combined, npcs = 30, verbose = T)
ElbowPlot(all_sample.combined)
all_sample.combined <- RunUMAP(all_sample.combined, seed.use = -1, reduction = "pca", dims = 1:20,return.model=T)
all_sample.combined <- FindNeighbors(all_sample.combined, reduction = "pca", dims = 1:20)
all_sample.combined <- FindClusters(all_sample.combined,  resolution = c(1))
all_sample.combined <- SetIdent(all_sample.combined,value = 'subtype')
DimPlot(all_sample.combined,split.by  = 'cell.type',label = T)

all_sample.combined$subtype <- ordered(all_sample.combined$subtype,
                                       levels=c('Precursor','Fundic Epi','Antral Epi','Gland','Fetal_stomach_precursor','Fetal_stomach_fundus','Fetal_stomach_antrum','Fetal_stomach_gland'))

vitro <- subset(all_sample.combined,split=='in vitro')
p1 <- DimPlot(vitro,group.by = 'celltype',cols = c('#E95C59', '#E4C755', '#58A4C3', '#23452F'),pt.size = 1)+ scale_x_reverse()

vivo <- subset(all_sample.combined,split=='In vivo')
p2 <- DimPlot(vivo,group.by = 'celltype',cols = c('#E95C59', '#E4C755', '#58A4C3', '#23452F'),pt.size = 2.5)+ scale_x_reverse()
p1+p2

# Supplementary Figure 5 spearman correlation
Idents(all_sample.combined) <- 'celltype'
exp<- AverageExpression(all_sample.combined)
data <- data.frame(exp$integrated)
numeric_cols <- sapply(data, is.numeric)
data[numeric_cols] <- lapply(data[numeric_cols], function(x) ifelse(is.na(x), 0, x))

num=names(tail(sort(apply(data, 1, sd)),88))  
num <- na.omit(match(num,rownames(data)))
data <- cor(data, method= "spearman")
spearman_dist_matrix <- as.dist(1 - data)

hc <- hclust(spearman_dist_matrix)
bk<-c(seq(0.6,0.9,by=0.03),seq(0.91,0.934,by=0.001),seq(0.935,0.96,by=0.01),seq(0.961,0.98,by=0.001),seq(0.981,1,by=0.01))
pheatmap(data, 
         display_numbers = F,
         cluster_rows = hc,
         cluster_cols = hc,angle_col = '315',
         legend_breaks=c(0.6,1),
         breaks=bk,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(60))



###### In vitro and in vivo expression of ECM and secreting protein genes in Epi lineage 
# Article:Single-cell transcriptomic characterization of a gastrulating human embryo 
anno <- readRDS(file = 'annot_umap.rds')
mat <- readRDS(file = 'expression_values.rds')
mat <- t(mat)
colnames(mat) <- anno$cell_name

lite <- CreateSeuratObject(counts = mat,meta.data = anno)
lite <- NormalizeData(lite)
lite@assays$RNA$data <- lite@assays$RNA$counts
lite <- FindVariableFeatures(lite, selection.method = "vst", nfeatures = 2000)
lite <- ScaleData(lite, features = rownames(lite))
lite <- RunPCA(lite, features = VariableFeatures(object = lite))
ElbowPlot(lite)
lite <- RunUMAP(lite, dims = 1:20)
umap <- anno[,c(2,3)]
rownames(umap) <- anno$cell_name
colnames(umap) <- c('umap_1','umap_2')
umap <- as.matrix(umap)
lite@reductions$umap@cell.embeddings <- umap

DimPlot(lite,group.by = 'cluster_id',label = T)
DimPlot(lite,group.by = 'sub_cluster',label = T)

# Article:Charting human development using a multiendodermal organ atlas and organoid models
load(file = 'stomach_epi.rdata')
DimPlot(epi, reduction = "umap",group.by = 'Cell_type',label = T,raster = F,cols = my36colors)
epi$lab <- epi$Cell_type
epi$lab <- plyr::mapvalues(x=epi$lab,from=c('Antral gland cell','Antral surface mucous cell','Antrum proliferative mucous cell','Antrum stem cell',
                                            'Corpus mucous neck cell','Corpus stem cell','Corpus surface mucous cell','Gastric epithelial precursor','Proliferative gastric epithelial precursor'),
                           to=c('Fetal_stomach','Fetal_stomach','Fetal_stomach','Fetal_stomach',
                                'Fetal_stomach','Fetal_stomach','Fetal_stomach','Fetal_stomach','Fetal_stomach'))

de <- lite[,lite$sub_cluster %in% c('DE(NP)','DE(P)')] 
mg <- merge(de,epi)
mg$lab[1:54] <- 'DE'
mg$lab <- ordered(mg$lab,levels=c('DE','Fetal_stomach','Enteroendocrine'))
DefaultAssay(mg) <- 'RNA'
mg <- ScaleData(mg,features = rownames(mg))

epi.lineage <- readRDS(file = 'epi.rds') # in vitro
DefaultAssay(epi.lineage) <- 'RNA'

# Extended Data Figure8
secreted.gene <- c('BMP2','FGF2','FGF14','IGF2')
p1 <- DotPlot(epi.lineage,features = rev(secreted.gene),
              group.by = 'cell.type',cols = c('lightgray','red'))+RotatedAxis()+coord_flip()
p2 <- DotPlot(mg,features = rev(secreted.gene),
              group.by = 'lab',cols = c('lightgray','red'))+RotatedAxis()+coord_flip()
p1|p2

ecm.gene <- c('COL1A2','FBLN1','COL4A2','COL12A1','HSPG2')
p1 <- DotPlot(epi.lineage,features = rev(ecm.gene),
              group.by = 'cell.type',cols = c('lightgray','red'))+RotatedAxis()+coord_flip()
p2 <- DotPlot(mg,features = rev(ecm.gene),
              group.by = 'lab',cols = c('lightgray','red'))+RotatedAxis()+coord_flip()
p1|p2



###### Extended Data Figure7 ICC in multi-lineage PFG 
# Article:Primate gastrulation and early organogenesis at single-cell resolution
meta <- read.csv(file = './MFE56636-meta.csv')
rawdata <- Read10X('./filtered_feature_bc_matrix')
mfe <- Read10X('./MFE56636MTX', gene.column = 1)
sce <- CreateSeuratObject(counts = mfe, min.cells = 3, min.features = 200) 

identical(meta$X,colnames(sce))
rownames(meta) <- meta$X
meta <- meta[,-1]
sce <- AddMetaData(sce,metadata = meta)
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
ElbowPlot(sce)
sce <- RunUMAP(sce, dims = 1:20)
cell.embeddings <- meta[,c('UMAP_1','UMAP_2')]
cell.embeddings <- as.matrix(cell.embeddings)
colnames(cell.embeddings) <- c('umap_1','umap_2')
sce@reductions$umap@cell.embeddings <- cell.embeddings
col <- levels(as.factor(sce$color))
DimPlot(sce, reduction = "umap", group.by = "cell_type",raster=FALSE,label = T,cols = col) 

load(file = 'integrated.RData')
d7 <- sample.integrated[, sample.integrated$day %in% 'D7'] 
d7_negative <- subset(x = d7, subset = FOXA2 == 0) 
d7_negative <- subset(x = d7_negative, subset = CDH1 == 0) 
sox2_gata4 <- subset(x = d7_negative, subset = GATA4 > 0|SOX2>0) 
data <- subset(sox2_gata4,subset=GATA4>0&SOX2>0)
num <- match(colnames(data),colnames(sox2_gata4))
data2 <- sox2_gata4[,-num]

# MapQuery
query_obj <- data2
DimPlot(query_obj, reduction = "umap", group.by = "cell.type",
        label = F,repel = T,label.size = 6,raster=FALSE,cols = c("#5050FFFF","#CE3D32FF","#466983FF","#A6CEE3","#1A0099FF"))

ref_obj <- sce
ref_obj[["umap.new"]] <- CreateDimReducObject(embeddings = ref_obj[["umap"]]@cell.embeddings, key = "UMAPnew_")
# set UMAP models
umap.new.model <- list()
umap.new.model$n_epochs <- 500
umap.new.model$alpha <-1
umap.new.model$method <- "umap"
umap.new.model$negative_sample_rate <- 5
umap.new.model$gamma <- 1
umap.new.model$approx_pow <- 0
umap.new.model$metric$cosine <- list()
umap.new.model$embedding <- ref_obj[["umap.new"]]@cell.embeddings
ab_param <- uwot:::find_ab_params(spread = 1, min_dist = 0.3)
umap.new.model$a <- ab_param["a"]
umap.new.model$b <- ab_param["b"]
ref_obj[["umap.new"]]@misc$model <- umap.new.model
Misc(object =ref_obj[['umap.new']], slot = "model")$n_neighbors <- 30
Misc(ref_obj[["umap.new"]], slot="model")$num_precomputed_nns <- 1

ref_obj <- FindClusters(ref_obj,resolution = 0.5)
DimPlot(ref_obj,group.by = 'integrated_snn_res.0.5',label=T,repel=T)
DimPlot(ref_obj,group.by = 'C38A2',label=T,repel=T)

name <- data.frame(colnames(ref_obj),ref_obj$integrated_snn_res.0.5)
cell <- data.frame(c(0:19),c('Ectoderm','Mesoderm','Mesoderm','Endoderm','Mesoderm','EPI','Mesoderm','Endoderm','Mesoderm','Ectoderm','Mesoderm','Mesoderm','Ectoderm','Mesoderm','Mesoderm','Ectoderm','Endoderm','Endoderm','Mesoderm','Mesoderm'))
colnames(cell) <- c("ref_obj.integrated_snn_res.0.5",'celltype')
name <- join(name,cell)
ref_obj$celltype <- name$celltype
DimPlot(ref_obj,group.by = 'celltype',label=T,repel=T)

anchors <- FindTransferAnchors(
  reference = ref_obj,
  query = query_obj,
  k.filter = NA ,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:20)

test_obj <- MapQuery(
  anchorset = anchors,
  query = query_obj,
  reference = ref_obj,
  refdata = list(
    Major = "celltype"),
  reference.reduction = "pca",
  reduction.model = "umap.new"
)

ref_obj$id <- 'reference'
test_obj$id <- 'query'

refquery <- merge(ref_obj,test_obj)
refquery[["umap"]] <- merge(ref_obj[["umap"]],
                            test_obj[["ref.umap"]])
refquery$cell <- 'Other cells'
D <- subset(refquery,cell.type=="DE")
E <- subset(refquery,cell.type=="NE")
h <- subset(refquery,cell.type=="hPSC")
M<- subset(refquery,cell.type=="Mesenchymal")
U<- subset(refquery,cell.type=="Unidentified")
D <- match(colnames(D),colnames(refquery))
refquery$cell[D] <- "DE"
E <- match(colnames(E),colnames(refquery))
refquery$cell[E] <- "NE"
h <- match(colnames(h),colnames(refquery))
refquery$cell[h] <- "hPSC"
M <- match(colnames(M),colnames(refquery))
refquery$cell[M] <- "Mesenchymal"
U <- match(colnames(U),colnames(refquery))
refquery$cell[U] <- "Unidentified"
p1 <- DimPlot(refquery, group.by = 'cell', shuffle = TRUE,label=F,repel = T,cols = c('#D3D3D3',"#CE3D32FF","#5050FFFF","#1A0099FF","#466983FF","#A6CEE3"),pt.size = 1 ,order = c("NE" ,"Mesenchymal","Unidentified","hPSC" ,"DE","Other cells"))
      
refquery$cell <- 'Other cells'
EC <- subset(refquery,predicted.Major=="Ectoderm")
EP <- subset(refquery,predicted.Major=="EPI")
M <- subset(refquery,predicted.Major=="Mesoderm")
EN <- subset(refquery,predicted.Major=="Endoderm")
EC <- match(colnames(EC),colnames(refquery))
refquery$cell[EC] <- "Ectoderm"
EP <- match(colnames(EP),colnames(refquery))
refquery$cell[EP] <- "EPI"
M <- match(colnames(M),colnames(refquery))
refquery$cell[M] <- "Mesoderm"
EN <- match(colnames(EN),colnames(refquery))
refquery$cell[EN] <- "Endoderm"

p2 <- DimPlot(ref_obj, group.by = 'celltype', shuffle = TRUE,label=F,cols = c("#A6CEE3",'#D8E88D',"#5050FFFF","#466983FF"),pt.size = 1,)




###### Extended Data Figure 7 NE in gastroids and the neuroepithelium in week 3 human embryos
# Article:The single-cell and spatial transcriptional landscape of human gastrulation and brain
data <- Read10X(data.dir = "GSM4695826_PCW3-01/")
data <- CreateSeuratObject(counts = data, project = "PC3", min.cells = 3, min.features = 200)
meta <- read.csv('meta.csv',header = T)

num <- match(meta$barcode,colnames(data))
data <- data[,num]
data$celltype <- meta$cell_type

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
data<- subset(data, subset = nFeature_RNA > 200  & percent.mt < 5)

data <- NormalizeData(data,verbose = TRUE)
data <- FindVariableFeatures(data,selection.method = 'vst',nfeatures = 2000)
data <- ScaleData(data)
data <- RunPCA(data, npcs = 30, verbose = T)
ElbowPlot(data)
data <- RunUMAP(data, seed.use = -1, reduction = "pca", dims = 1:20,return.model = TRUE)
data <- FindNeighbors(data, reduction = "pca", dims = 1:20)
data <- FindClusters(data,  resolution = 2)
DimPlot(data, reduction = "umap", group.by = "celltype",label = T,repel = T,raster=FALSE)
Neural <- subset(data,celltype=='Neural tbue 1'|celltype=='Neural tbue 2'|celltype=='Neural tbue 3')

DefaultAssay(Neural) <- 'RNA'
Neural <- FindVariableFeatures(Neural,selection.method = 'vst',nfeatures = 2000)
Neural <- ScaleData(Neural)
Neural <- RunPCA(Neural, npcs = 30, verbose = T)
ElbowPlot(Neural)
Neural <- RunUMAP(Neural, seed.use = -1, reduction = "pca", dims = 1:20,return.model = TRUE)
Neural <- FindNeighbors(Neural, reduction = "pca", dims = 1:20)
Neural <- FindClusters(Neural,  resolution = 0.5)
DimPlot(Neural, reduction = "umap", group.by = "RNA_snn_res.0.5",label = T,repel = T,raster=FALSE)

Neural <- subset(Neural,RNA_snn_res.0.5==0|RNA_snn_res.0.5==1|RNA_snn_res.0.5==2|RNA_snn_res.0.5==3|RNA_snn_res.0.5==6)
ta1 <- data.frame(colnames(Neural),Neural$RNA_snn_res.0.5)
ta2 <- data.frame(c(0,1,2,3,6),c('Hindbrain','spinalcord','forebrain','midbrain','spinalcord'))
colnames(ta2) <- c("Neural.RNA_snn_res.0.5",'cell.type')
ta <- join(ta1,ta2)
Neural$cell.type <- ta$cell.type
DimPlot(Neural, reduction = "umap", group.by = "cell.type",label = T,repel = T,raster=FALSE)

# gastroids data
load('integrated.rdata')
data <- subset(sample.integrated,cell.type=='NPC'|cell.type=='NE'|cell.type=='Migratory ENCC'|cell.type=='Premigratory ENCC'|cell.type=='Neuron')
data <- subset(data,cell.type=='NE')
data <- subset(data,downsample=222)
data@assays$integrated <- NULL 
data$orig.ident <- 'NE'
DefaultAssay(data) <- 'RNA'

# integrate
all.data <- merge(Neural,data)
alldata.list <- SplitObject(all.data, split.by = "orig.ident")
alldata.list <- lapply(X = alldata.list, FUN = function(x) {
  DefaultAssay(x)="RNA"  
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000)
})
all_sample.anchors <- FindIntegrationAnchors(object.list = alldata.list,dims =1:20,anchor.features = 1000)
Neural <- SetIdent(Neural,value = 'cell.type')
DEG <- FindAllMarkers(Neural,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.1)
DEG <- filter(DEG, DEG$p_val_adj < 0.05)
all_sample.combined <- IntegrateData(anchorset = all_sample.anchors, dims = 1:20,features.to.integrate = DEG$gene)
DefaultAssay(all_sample.combined) <- "integrated"

all_sample.combined <- ScaleData(all_sample.combined,verbose = T)
all_sample.combined <- RunPCA(all_sample.combined, npcs = 30, verbose = T)
ElbowPlot(all_sample.combined)
all_sample.combined <- RunUMAP(all_sample.combined, seed.use = -1, reduction = "pca", dims = 1:20,return.model=T)
all_sample.combined <- RunTSNE(all_sample.combined, seed.use = 2, reduction = "pca", dims = 1:20, check_duplicates = FALSE)
all_sample.combined <- FindNeighbors(all_sample.combined, reduction = "pca", dims = 1:20)
all_sample.combined <- FindClusters(all_sample.combined,  resolution = c(1))
all_sample.combined <- SetIdent(all_sample.combined,value = 'cell.type')
DimPlot(all_sample.combined,split.by  = 'cell.type',label = T)

# spearman correlation
Idents(all_sample.combined) <- all_sample.combined$cell.type
exp <- AverageExpression(all_sample.combined)

data <- exp$integrated
num=names(tail(sort(apply(data, 1, sd)),1000))
num <- match(num,rownames(data))
data <- data[num,]
data <- data.frame(data)
data <- cor(data, method= "spearman")

for (i in 1:5) {
  for (j in 1:5) {
    data[i,j] <- (data[i,j]-0.8961564)/(1-0.8961564)
  }  
}

bk <- c(seq(0,0.2,by=0.03),seq(0.21,0.3,by=0.03),seq(0.31,0.4,by=0.03),seq(0.41,0.5,by=0.03),seq(0.51,0.72,by=0.03),seq(0.73,0.74,by=0.03),seq(0.75,0.76,by=0.005),seq(0.77,0.78,by=0.005),seq(0.79,0.8,by=0.03),seq(0.81,0.90,by=0.03),seq(0.91,1,by=0.03))###0.02(鏈夐鑹?)鎴栬�?0.03(鐧?)
pheatmap(data, 
         display_numbers = F,
         cluster_rows = T,
         cluster_cols = T,angle_col = '315',
         breaks=bk,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(42))



