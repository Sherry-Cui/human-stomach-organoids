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
all <- CreateSeuratObject(count,meta.data = info) #28691 155232
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
load(file = 'all.sample_d16_epi.rdata')
DefaultAssay(all.sample_d16_epi) <- 'RNA'
all.sample_d16_epi$Corrected_organ_group <- ordered(all.sample_d16_epi$Corrected_organ_group,levels=c('Intestine','Esophagus','Stomach','d16'))
gene <- c('CLDN18','SOX2','CDH1', # Stomach epithelium
          'TP63','KRT5', # Esophagus epithelium
          'CDX2') # Intestine epithelium
DotPlot(all.sample_d16_epi,features = gene,cols = c('lightgrey','red'),group.by = 'Corrected_organ_group')+RotatedAxis()

# in vivo/in vitro UMAP and correlation heatmap
set.seed(20230518)
sample <- subset(sample.integrated,downsample=2000) #23068
list <- list(stomach,sample)
for (i in 1:length(list)){
  list[[i]] <- NormalizeData(list[[i]], verbose = FALSE)
  list[[i]] <- FindVariableFeatures(list[[i]], selection.method = "vst", nfeatures = 2000, verbose = T)
}
anchors <- FindIntegrationAnchors(object.list = list,k.filter=200, anchor.features = 2000)
integrated <- IntegrateData(anchorset = anchors)
integrated <- ScaleData(integrated, verbose = T)
integrated <- RunPCA(integrated,verbose = T)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:20)
integrated$lab1[1:13372] <- integrated$Major_cell_type[1:13372]
integrated$lab1[1:13372] <- paste('lite',integrated$lab1[1:13372],sep = '_')
integrated$lite <- 'in vitro'
integrated$lite[1:13372] <- 'in vivo'
integrated$day[1:13372] <- 'in vivo'

d16 <- integrated[, integrated$day %in% c("D16",'In vivo')] 
mt.lab <- as.data.frame(t(as.matrix(GetAssayData(d16, assay = "integrated", slot = "data"))))
group.by <- 'lab1'
mt.lab <- aggregate(mt.lab, by=list(d16@meta.data[[group.by]]), FUN="mean")
rownames(mt.lab) <- mt.lab$Group.1
mt.lab <- t(mt.lab[,-1])
mt.lab <- as.data.frame(mt.lab[,-c(1,2,3,11,12,17,18)]) 
cor <- cor(mt.lab,method = "spearman")
# Figure3
pheatmap(cor,display_numbers=F,cluster_cols=T,cluster_rows=T,clustering_method = 'average') 

vitro <- integrated[, integrated$lite %in% "in vitro"] 
vivo <- integrated[, Idents(integrated) %in% "in vivo"] 
cols <- c("#ff4a46","#008941","#ffdbe5","#0000a6","#eec3ff","#456d75")
d16 <- vitro[, vitro$day %in% "D16"] 
# Figure3
p1 <- DimPlot(d16, group.by = "lab1", pt.size=0.1,label = T,reduction = 'umap')+
  xlim(-13, 10) + ggtitle("In vitro(D16)") + scale_color_manual(values = col)
p2 <- DimPlot(vivo, group.by = "lab1", pt.size=0.1,label = T) + ggtitle("In vivo") + scale_color_manual(values = cols)
p1|p2

# Extended Data Figure5 UMAP
integrated$day <- ordered(integrated$day,levels=c('D4','D7','D10','D13','D16','In vivo'))
DimPlot(integrated, group.by = "lab1",split.by = 'day', pt.size=0.1,label = F,ncol = 3) + ggtitle("Day") + scale_color_manual(values = col)

# Figure3 FeaturePlot
DefaultAssay(d16) <- 'RNA'
DefaultAssay(vivo) <- 'RNA'
p1 <- FeaturePlot(d16,features = c('CDH1','COL1A2','MAP2','FLT1'),cols = c('lightgrey','red'),ncol = 1,label = F)
p1 <- FeaturePlot(d16,features = 'CDH1',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p2 <- FeaturePlot(d16,features = 'COL1A2',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p3 <- FeaturePlot(d16,features = 'MAP2',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p4 <- FeaturePlot(d16,features = 'FLT1',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p <- p1/p2/p3/p4
p5 <- FeaturePlot(vivo,features = c('CDH1','COL1A2','MAP2','FLT1'),cols = c('lightgrey','red'),ncol = 1,label = F)
p5 <- FeaturePlot(vivo,features = 'CDH1',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p6 <- FeaturePlot(vivo,features = 'COL1A2',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p7 <- FeaturePlot(vivo,features = 'MAP2',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p8 <- FeaturePlot(vivo,features = 'FLT1',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p|(p5/p6/p7/p8)


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

# Dotplot
load(file = 'allepi.reintegrate.rdata')

gene <- c('EZH2','FOXO3','TCF4', # Precursor
          'SOX2','IRX2','IRX3','IRX5','NR2F2', # Fundus
          'PDX1','MEIS2',# Antrum
          'TFF2','UPK1B', 'PLAC8'# Secreting
) 
DefaultAssay(counts) <- 'RNA'
# Figure4
p1 <- DotPlot(counts,features = gene,group.by = 'lab')+RotatedAxis()+coord_flip()+ 
  scale_color_gradientn(colors = c("grey85", brewer.pal(9, "OrRd"))) 
p2 <- DotPlot(epi.subset,features = gene,group.by = 'subtype')+RotatedAxis()+coord_flip()+ 
  scale_color_gradientn(colors = c("grey85", brewer.pal(9, "OrRd")))
p1|p2

# Figure4 UMAP
DimPlot(counts,group.by = 'lab',cols =c('#E95C59', '#E4C755', '#58A4C3', '#23452F'),label = F)
DimPlot(epi.subset,group.by = 'subtype',cols =c('#E95C59', '#E4C755', '#58A4C3', '#23452F'),label = F)

# FeaturePlot
d16 <- subset(counts,day=='D16')
# Figure4
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

# Figure6
FeaturePlot(epi.subset,features = 'NR2F2',slot = 'scale.data',pt.size = 0.2)+
  scale_colour_gradientn(colours = colorRampPalette(c("#DCDCDC",'#FDF5E6','#FFA07A',"#fc2c14",'#CD0000'))(30),labels = c("low", "high"),
                         breaks=c(-0.95,3))
FeaturePlot(d16,features = 'NR2F2',slot = 'scale.data',pt.size = 0.2)+
  scale_colour_gradientn(colours = colorRampPalette(c("#DCDCDC","#DCDCDC",'#FDF5E6','#FFA07A',"#fc2c14",'#CD0000','#CD0000'))(30),labels = c("low", "high"),
                         breaks=c(-0.73,4.3))

# integrate
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
Cols <- c('#E95C59', '#E4C755', '#58A4C3', '#23452F','#E95C59', '#E4C755', '#58A4C3' ,'#23452F')
names(Cols) <- levels(all_sample.combined$subtype)
# Figure4
DimPlot(all_sample.combined,group.by = 'subtype',cols = Cols,label = F,split.by = 'split')

# spearman correlation
Idents(all_sample.combined)<- all_sample.combined$subtype
exp<- AverageExpression(all_sample.combined)
data <- exp$integrated
num <- names(tail(sort(apply(data, 1, sd)),1000))
num <- match(num,rownames(data))
data <- data[num,]
data <- data.frame(data)
data <- cor(data, method= "spearman")

for (i in 1:8) {
  for (j in 1:8) {
    data[i,j] <- (data[i,j]-0.5443106)/(1-0.5443106)
  }  
}
colnames(data) <- c("Antral Epi","Precursor","Fundic Epi","Gland","Fetal_stomach_gland","Fetal_stomach_fundus","Fetal_stomach_antrum","Fetal_stomach_precursor")
rownames(data) <- c("Antral Epi","Precursor","Fundic Epi","Gland","Fetal_stomach_gland","Fetal_stomach_fundus","Fetal_stomach_antrum","Fetal_stomach_precursor")
bk <- c(seq(0,0.49,by=0.01),seq(0.5,1,by=0.01))
# Supplementary Figure6
pheatmap(data,display_numbers = F,color = colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))(100),
         legend_breaks=c(0,1),
         breaks=bk,
         legend_labels = c('low','high'),angle_col = '315',fontsize = 12)



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



###### ICC in multi-lineage PFG 
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
# Extended Data Figure7
DimPlot(sox2_gata4,group.by = 'cell.type',label = F,cols = c("#5050FFFF","#CE3D32FF","#466983FF","#A6CEE3","#1A0099FF"))+ggtitle(' ')

# MapQuery
query_obj <- sox2_gata4
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
refquery$cell <- ' '
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

refquery$cell <- ' '
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
# Extended Data Figure7
p1 <- DimPlot(refquery, group.by = 'cell', shuffle = TRUE,label=F,repel = T,cols = c('#D3D3D3',"#5050FFFF","#CE3D32FF","#466983FF","#1A0099FF","#A6CEE3"),order = c("Early NPC","Unidentified","Mesenchymal","DE","hPSC"),pt.size = 1)
p2 <- DimPlot(ref_obj, group.by = 'celltype', shuffle = TRUE,label=F,cols = c("#A6CEE3",'#D8E88D',"#5050FFFF","#466983FF"),pt.size = 1)



###### NE in gastroids and the neuroepithelium in week 3 human embryos
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
rev(brewer.pal(n = 11, name ="RdBu"))
bk <- c(seq(0,0.2,by=0.03),seq(0.21,0.3,by=0.03),seq(0.31,0.4,by=0.03),seq(0.41,0.5,by=0.03),seq(0.51,0.72,by=0.03),seq(0.73,0.74,by=0.03),seq(0.75,0.76,by=0.005),seq(0.77,0.78,by=0.005),seq(0.79,0.8,by=0.03),seq(0.81,0.90,by=0.03),seq(0.91,1,by=0.03))
# Extended Data Figure6
pheatmap(data,display_numbers = F,color = colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))(43),
              legend_breaks=c(0,0.999),
              breaks=bk,
              legend_labels = c('low','high'),fontsize = 20,angle_col = 315)



