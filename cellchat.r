library(reticulate)
ad <- import("anndata", convert = FALSE)
ad_object <- ad$read_h5ad("/data/AMLscRNA/mergeblood/20220501/cellchat/RESb.h5ad")
meta.data <- py_to_r(ad_object$obs)
data.input <- t(py_to_r(ad_object$X))
rownames(data.input) <- rownames(py_to_r(ad_object$var))
colnames(data.input) <- rownames(py_to_r(ad_object$obs))
library(Matrix)
data.input<- as(as.matrix(data.input), "dgCMatrix")
library(Seurat)
setwd("/data/AMLscRNA/mergeblood/20220501/cellchat")

library(CellChat)
library(dplyr)
library(tidyverse)
CellChatDB <- CellChatDB.human
Secreted.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 

RESb <- createCellChat(object = data.input, meta =meta.data, group.by = "CellTypeN")
RESb<- addMeta(RESb, meta =meta.data)
RESb<- setIdent(RESb, ident.use = "CellTypeN")
RESbSize <- as.numeric(table(RESb@idents))
RESb@DB <- Secreted.use
RESb <- subsetData(RESb)
RESb<- identifyOverExpressedGenes(RESb)
RESb <- identifyOverExpressedInteractions(RESb)
RESb <- projectData(RESb, PPI.human)  
RESb <- computeCommunProb(RESb, raw.use = TRUE)
RESb <- computeCommunProbPathway(RESb)
RESb <- aggregateNet(RESb)
RESb <- netAnalysis_computeCentrality(RESb, slot.name = "netP")
saveRDS(RESb,"RESb.rds")

ad <- import("anndata", convert = FALSE)
ad_object <- ad$read_h5ad("/data/AMLscRNA/mergeblood/20220501/cellchat/RESa.h5ad")
meta.data <- py_to_r(ad_object$obs)
data.input <- t(py_to_r(ad_object$X))
rownames(data.input) <- rownames(py_to_r(ad_object$var))
colnames(data.input) <- rownames(py_to_r(ad_object$obs))
data.input<- as(as.matrix(data.input), "dgCMatrix")

RESa <- createCellChat(object =data.input, meta =meta.data, group.by = "CellTypeN")
RESa<- addMeta(RESa, meta =meta.data)
RESa<- setIdent(RESa, ident.use = "CellTypeN")
RESaSize <- as.numeric(table(RESa@idents))
#Secreted.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
RESa@DB <- Secreted.use
RESa <- subsetData(RESa)
RESa<- identifyOverExpressedGenes(RESa)
RESa <- identifyOverExpressedInteractions(RESa)
RESa <- projectData(RESa, PPI.human)  
RESa <- computeCommunProb(RESa, raw.use = TRUE)
RESa <- computeCommunProbPathway(RESa)
RESa <- aggregateNet(RESa)
RESa <- netAnalysis_computeCentrality(RESa, slot.name = "netP")
saveRDS(RESa,"RESa.rds")

ad <- import("anndata", convert = FALSE)
ad_object <- ad$read_h5ad("/data/AMLscRNA/mergeblood/20220501/cellchat/UNb.h5ad")
meta.data <- py_to_r(ad_object$obs)
data.input <- t(py_to_r(ad_object$X))
rownames(data.input) <- rownames(py_to_r(ad_object$var))
colnames(data.input) <- rownames(py_to_r(ad_object$obs))
library(Matrix)
data.input<- as(as.matrix(data.input), "dgCMatrix")

UNb <- createCellChat(object =data.input, meta =meta.data, group.by = "CellTypeN")
UNb<- addMeta(UNb, meta =meta.data)
UNb<- setIdent(UNb, ident.use = "CellTypeN")
UNbSize <- as.numeric(table(UNb@idents))
UNb@DB <- Secreted.use
UNb <- subsetData(UNb)
UNb<- identifyOverExpressedGenes(UNb)
UNb <- identifyOverExpressedInteractions(UNb)
UNb <- projectData(UNb, PPI.human)  
UNb <- computeCommunProb(UNb, raw.use = TRUE)
UNb <- computeCommunProbPathway(UNb)
UNb <- aggregateNet(UNb)
UNb <- netAnalysis_computeCentrality(UNb, slot.name = "netP")
saveRDS(UNb,"UNb.rds")

ad <- import("anndata", convert = FALSE)
ad_object <- ad$read_h5ad("/data/AMLscRNA/mergeblood/20220501/cellchat/UNa.h5ad")
meta.data <- py_to_r(ad_object$obs)
data.input <- t(py_to_r(ad_object$X))
rownames(data.input) <- rownames(py_to_r(ad_object$var))
colnames(data.input) <- rownames(py_to_r(ad_object$obs))
library(Matrix)
data.input<- as(as.matrix(data.input), "dgCMatrix")

UNa <- createCellChat(object =data.input, meta =meta.data, group.by = "CellTypeN")
UNa<- addMeta(UNa, meta =meta.data)
UNa<- setIdent(UNa, ident.use = "CellTypeN")
UNaSize <- as.numeric(table(UNa@idents))
UNa@DB <- Secreted.use
UNa <- subsetData(UNa)
UNa<- identifyOverExpressedGenes(UNa)
UNa <- identifyOverExpressedInteractions(UNa)
UNa <- projectData(UNa, PPI.human)  
UNa <- computeCommunProb(UNa, raw.use = TRUE)
UNa <- computeCommunProbPathway(UNa)
UNa <- aggregateNet(UNa)
UNa <- netAnalysis_computeCentrality(UNa, slot.name = "netP")
saveRDS(UNa,"UNa.rds")

RESa<-readRDS(file = "RESa.rds")
RESb<-readRDS(file = "RESb.rds")
object.list1 <- list(RESb = RESb, RESa = RESa)
cellchat1 <- mergeCellChat(object.list1, add.names = names(object.list1))

UNa<-readRDS(file = "UNa.rds")
UNb<-readRDS(file = "UNb.rds")
group.new =union(levels(UNb@idents),levels(UNa@idents))
UNb <- liftCellChat(UNb, group.new)
UNa <- liftCellChat(UNa, group.new)
object.list <- list(UNb = UNb, UNa = UNa)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

group.new =union(levels(UNb@idents),levels(RESb@idents))
UNb <- liftCellChat(UNb, group.new)
RESb <- liftCellChat(RESb, group.new)
object.list <- list(RESb = RESb,UNb = UNb)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use=c("1"="black","2"="red"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",color.use=c("1"="black","2"="red"))
gg1 + gg2
pdf("pathway_change_circlemap_before.pdf",width=13,height=8)
par(mfrow = c(1,2),mar=c(.3,.3,.3,.3),xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
pdf("pathway_change_heatmap_RES.pdf",width=13,height=8)
pdf("interaction_change_heatmap_UN.pdf",width=13,height=8)
pdf("interaction_change_heatmap_before.pdf",width=13,height=8)
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
dev.off()

#pathwaycompare 
rankSimilarity(cellchat, type = "functional")
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use=c("RESb"="black","RESa"="red"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c("RESb"="black","RESa"="red"))
pdf("signalpathway_RES.pdf",width=14,height=8)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use=c("UNb"="black","UNa"="red"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c("UNb"="black","UNa"="red"))
pdf("signalpathway_UN.pdf",width=14,height=8)
pdf("signalpathway_before.pdf",width=14,height=8)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use=c("RESb"="black","UNb"="red"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c("RESb"="black","UNb"="red"))
gg1 + gg2
dev.off()

write.table(gg1$data,"/data/AMLscRNA/mergeblood/20220501/cellchat/UN_compare.bed",sep="\t",quote = F,row.names = F)

library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
source("/data/carT/lung5/netAnalysis_signalingRole_heatmap_mode.r")
ht1 = netAnalysis_signalingRole_heatmap_modi(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width =10, height = 14)
ht2 = netAnalysis_signalingRole_heatmap_modi(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height =14)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
source("/data/carT/lung5/netAnalysis_signalingRole_heatmap_3.r")
ht3 = ht2@matrix-ht1@matrix
g1<-netAnalysis_signalingRole_heatmap_ch(ht3,  width =10, height = 14,color.heatmap = c("#2166ac","red"), pattern = "outgoing")
write.table(ht3,"UN_cell_pathway_matrix_outgoing",sep="\t",quote = F,row.names = T)

ht1 = netAnalysis_signalingRole_heatmap_modi(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 14, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap_modi(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height =14, color.heatmap = "GnBu")
ht3 = ht2@matrix-ht1@matrix
#pdf("/data/RESa/PB5/incoming_patterndiff.pdf")
g2<-netAnalysis_signalingRole_heatmap_ch(ht3,  width =10, height = 14,color.heatmap = c("#2166ac","red"), pattern = "incoming")
#draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
pdf("pathway_change_heatmap_before.pdf",width=13,height=8)
draw(g1 + g2, ht_gap = unit(0.5, "cm"))
dev.off()
write.table(ht3,"UN_cell_pathway_matrix_incoming",sep="\t",quote = F,row.names = T)


gg1 <- netVisual_bubble(cellchat, sources.use = c(1,3,8,9,10,11,13,14,17:20), targets.use =  c(1,3,8,9,10,11,13,14,17:20),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in RESa", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in RESb", angle.x = 45, remove.isolate = T)
gg1 + gg2

pdf("CXCL_before.pdf",width=13,height=8)
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}

pdf("LSPC_out_UNup.pdf",width=14,height=5)
netVisual_bubble(cellchat, sources.use = c(7,8,9), targets.use = c(1:4,6,10,12:15,17:18,20:21),comparison = c(1, 2), max.dataset = 2, title.name = " Increased signaling in AML cells of UNa", angle.x = 45, remove.isolate = T)
dev.off()

aa<-netVisual_bubble(cellchat, targets.use = 7, sources.use = c(1:4,5,10,12:15,17:18,20:24),  comparison = c(1, 2), angle.x = 45)

write.table(aa$data,"LSPC_cycle_LR.bed",sep="\t",quote = F,row.names = F)
aa<-netVisual_bubble(cellchat, targets.use = 8, sources.use = c(1:4,5,10,12:15,17:18,20:24),  comparison = c(1, 2), angle.x = 45)

write.table(aa$data,"LSPC_Primed_LR.bed",sep="\t",quote = F,row.names = F)
aa<-netVisual_bubble(cellchat, targets.use = 9, sources.use = c(1:4,5,10,12:15,17:18,20:24),  comparison = c(1, 2), angle.x = 45)

write.table(aa$data,"LSPC_Quiescent_LR.bed",sep="\t",quote = F,row.names = F)
aa<-netVisual_bubble(cellchat, targets.use = 5, sources.use = c(1:4,5,10,12:15,17:18,20:24),  comparison = c(1, 2), angle.x = 45)

write.table(aa$data,"GMP_like_LR.bed",sep="\t",quote = F,row.names = F)
