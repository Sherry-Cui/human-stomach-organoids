######## CellChat of day 10 cells -----------
library(CellChat)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(ggsci)

load(file = 'integrated.RData')
sce <- subset(sample.integrated, day == 'D10')
sce <- subset(sce, cellchatuse %in% c('Precursor', 'Fundic Epi', 'Antral Epi', 'Mesenchymal', 'NE', 'NPC', 'Neuron', 'ENCC'))
data.input  <- sce@assays$RNA@data
identity = data.frame(group =sce$cell.type , row.names = names(sce$cell.type)) 
unique(identity$group) 
cellchat <- createCellChat(data <- data.input) 
summary(cellchat)

cellchat <- addMeta(cellchat, meta = identity, meta.name = "cell.type")
cellchat <- setIdent(cellchat, ident.use = "cell.type") 
levels(cellchat@idents) 
table(cellchat@idents)

CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
head(cellchat@LR$LRsig)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Supplementary Figure7
netVisual_aggregate(cellchat, signaling = 'WNT',pt.title=20,vertex.label.cex = 1.7,arrow.width=0.8,arrow.size = 0.5)
netVisual_aggregate(cellchat, signaling = 'ncWNT',pt.title=20,vertex.label.cex = 1.7,arrow.width=0.8,arrow.size = 0.5)
netVisual_aggregate(cellchat, signaling = 'CXCL',pt.title=20,vertex.label.cex = 1.7,arrow.width=0.8,arrow.size = 0.5)
netVisual_aggregate(cellchat, signaling = 'FGF',pt.title=20,vertex.label.cex = 1.7,arrow.width=0.8,arrow.size = 0.5)
netVisual_aggregate(cellchat, signaling = 'PDGF',pt.title=20,vertex.label.cex = 1.7,arrow.width=0.8,arrow.size = 0.5)
netVisual_aggregate(cellchat, signaling = 'BMP',pt.title=20,vertex.label.cex = 1.7,arrow.width=0.8,arrow.size = 0.5)
netVisual_aggregate(cellchat, signaling = 'IGF',pt.title=20,vertex.label.cex = 1.7,arrow.width=0.8,arrow.size = 0.5)
netVisual_aggregate(cellchat, signaling = 'GAS',pt.title=20,vertex.label.cex = 1.7,arrow.width=0.8,arrow.size = 0.5)


netVisual_bubble(cellchat, remove.isolate = FALSE, signaling = c("PDGF",'BMP'), 
                 sources.use = c('Precursor','Antral Epi',"Fundic Epi"),targets.use =c('Mesenchymal','NE','NPC','Neuron','ENCC'))

netAnalysis_signalingRole_network(cellchat, signaling = 'PDGF', width = 10, height = 5 ,font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = 'BMP', width = 10, height = 5 ,font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = 'ncWNT', width = 10, height = 5 ,font.size = 10)

plotGeneExpression(cellchat, signaling = c('PDGF','BMP'))
plotGeneExpression(cellchat, signaling = 'WNT')

# Supplementary Figure6
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",arrow.width=0.8,arrow.size = 0.5)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F,title.name = "Interaction weights/strength",arrow.width=0.8,arrow.size = 0.5)

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.heatmap = "GnBu", width = 7, height = 8)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.heatmap = "GnBu", width = 7, height = 8)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

netVisual_bubble(cellchat, remove.isolate = FALSE, signaling = 'WNT', 
                 sources.use = c('NE','NPC'),
                 targets.use =c('Precursor','Antral Epi',"Fundic Epi",'Mesenchymal','NE','NPC','Neuron','ENCC'))+coord_flip()

# Figure5
netAnalysis_signalingRole_network(cellchat, signaling = 'WNT', width = 10, height = 5 ,font.size = 10)


saveRDS(cellchat,file = 'cellchat.rds')





