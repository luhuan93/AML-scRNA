library(Seurat)
setwd("/data/AMLscRNA")
library(dplyr)
P1b <- Read10X(data.dir = "/data/AMLscRNA/1026/outs/filtered_feature_bc_matrix")
P1b <- CreateSeuratObject(counts = P1b, project = "P1b", min.cells = 3, min.features = 200)
P1b[["percent.mt"]] <- PercentageFeatureSet(P1b, pattern = "^MT-")
VlnPlot(P1b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P1b <- subset(P1b, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15&nCount_RNA > 1000 & nCount_RNA < 25000)
P1b@meta.data$state<-"before"
P1b@meta.data$sample<-"P1"
KA<-read.table("/data/AMLscRNA/scrublet/P1b_doublet.txt",header = T,sep = ",")
a<-P1b@meta.data
a$barcode<-rownames(P1b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P1b@meta.data<-c
P1b_sub=subset(P1b,subset=(predicted_doublets=="False"))
dim(P1b_sub)

P1a <- Read10X(data.dir = "/data/AMLscRNA/1163/outs/filtered_feature_bc_matrix")
P1a <- CreateSeuratObject(counts = P1a, project = "P1a", min.cells = 3, min.features = 200)
P1a[["percent.mt"]] <- PercentageFeatureSet(P1a, pattern = "^MT-")
VlnPlot(P1a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P1a <- subset(P1a, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15&nCount_RNA > 1000 & nCount_RNA < 40000)
P1a@meta.data$state<-"after"
P1a@meta.data$sample<-"P1"
KA<-read.table("/data/AMLscRNA/scrublet/P1a_doublet.txt",header = T,sep = ",")
a<-P1a@meta.data
a$barcode<-rownames(P1a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P1a@meta.data<-c
P1a_sub=subset(P1a,subset=(predicted_doublets=="False"))
dim(P1a_sub)

P2b <- Read10X(data.dir = "/data/AMLscRNA/1059/outs/filtered_feature_bc_matrix")
P2b <- CreateSeuratObject(counts = P2b, project = "P2b", min.cells = 3, min.features = 200)
P2b[["percent.mt"]] <- PercentageFeatureSet(P2b, pattern = "^MT-")
VlnPlot(P2b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P2b <- subset(P2b, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15&nCount_RNA > 1000 & nCount_RNA < 30000)
P2b@meta.data$state<-"before"
P2b@meta.data$sample<-"P2"
KA<-read.table("/data/AMLscRNA/scrublet/P2b_doublet.txt",header = T,sep = ",")
a<-P2b@meta.data
a$barcode<-rownames(P2b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P2b@meta.data<-c
P2b_sub=subset(P2b,subset=(predicted_doublets=="False"))
dim(P2b_sub)

P2a <- Read10X(data.dir = "/data/AMLscRNA/1176/outs/filtered_feature_bc_matrix")
P2a <- CreateSeuratObject(counts = P2a, project = "P2a", min.cells = 3, min.features = 200)
P2a[["percent.mt"]] <- PercentageFeatureSet(P2a, pattern = "^MT-")
VlnPlot(P2a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P2a <- subset(P2a, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 15&nCount_RNA > 1000 & nCount_RNA < 45000)
P2a@meta.data$state<-"after"
P2a@meta.data$sample<-"P2"
KA<-read.table("/data/AMLscRNA/scrublet/P2a_doublet.txt",header = T,sep = ",")
a<-P2a@meta.data
a$barcode<-rownames(P2a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P2a@meta.data<-c
P2a_sub=subset(P2a,subset=(predicted_doublets=="False"))
dim(P2a_sub)

P3b <- Read10X(data.dir = "/data/AMLscRNA/1077/outs/filtered_feature_bc_matrix")
P3b <- CreateSeuratObject(counts = P3b, project = "P3b", min.cells = 3, min.features = 200)
P3b[["percent.mt"]] <- PercentageFeatureSet(P3b, pattern = "^MT-")
VlnPlot(P3b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P3b <- subset(P3b, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15&nCount_RNA > 1000 & nCount_RNA < 20000)
P3b@meta.data$state<-"before"
P3b@meta.data$sample<-"P3"
KA<-read.table("/data/AMLscRNA/scrublet/P3b_doublet.txt",header = T,sep = ",")
a<-P3b@meta.data
a$barcode<-rownames(P3b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P3b@meta.data<-c
P3b_sub=subset(P3b,subset=(predicted_doublets=="False"))
dim(P3b_sub)

P3a <- Read10X(data.dir = "/data/AMLscRNA/1191/outs/filtered_feature_bc_matrix")
P3a <- CreateSeuratObject(counts = P3a, project = "P3a", min.cells = 3, min.features = 200)
P3a[["percent.mt"]] <- PercentageFeatureSet(P3a, pattern = "^MT-")
VlnPlot(P3a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P3a <- subset(P3a, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15&nCount_RNA > 1000 & nCount_RNA < 25000)
P3a@meta.data$state<-"after"
P3a@meta.data$sample<-"P3"
KA<-read.table("/data/AMLscRNA/scrublet/P3a_doublet.txt",header = T,sep = ",")
a<-P3a@meta.data
a$barcode<-rownames(P3a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P3a@meta.data<-c
P3a_sub=subset(P3a,subset=(predicted_doublets=="False"))
dim(P3a_sub)

P4b <- Read10X(data.dir = "/data/AMLscRNA/1266/outs/filtered_feature_bc_matrix")
P4b <- CreateSeuratObject(counts = P4b, project = "P4b", min.cells = 3, min.features = 200)
P4b[["percent.mt"]] <- PercentageFeatureSet(P4b, pattern = "^MT-")
VlnPlot(P4b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P4b <- subset(P4b, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 10&nCount_RNA > 1000 & nCount_RNA < 20000)
P4b@meta.data$state<-"before"
P4b@meta.data$sample<-"P4"
KA<-read.table("/data/AMLscRNA/scrublet/P4b_doublet.txt",header = T,sep = ",")
a<-P4b@meta.data
a$barcode<-rownames(P4b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P4b@meta.data<-c
P4b_sub=subset(P4b,subset=(predicted_doublets=="False"))
dim(P4b_sub)

P4a <- Read10X(data.dir = "/data/AMLscRNA/1421/outs/filtered_feature_bc_matrix")
P4a <- CreateSeuratObject(counts = P4a, project = "P4a", min.cells = 3, min.features = 200)
P4a[["percent.mt"]] <- PercentageFeatureSet(P4a, pattern = "^MT-")
VlnPlot(P4a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P4a <- subset(P4a, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10&nCount_RNA > 1000 & nCount_RNA < 25000)
P4a@meta.data$state<-"after"
P4a@meta.data$sample<-"P4"
KA<-read.table("/data/AMLscRNA/scrublet/P4a_doublet.txt",header = T,sep = ",")
a<-P4a@meta.data
a$barcode<-rownames(P4a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P4a@meta.data<-c
P4a_sub=subset(P4a,subset=(predicted_doublets=="False"))
dim(P4a_sub)

P5b <- Read10X(data.dir = "/data/AMLscRNA/1152/outs/filtered_feature_bc_matrix")
P5b <- CreateSeuratObject(counts = P5b, project = "P5b", min.cells = 3, min.features = 200)
P5b[["percent.mt"]] <- PercentageFeatureSet(P5b, pattern = "^MT-")
VlnPlot(P5b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P5b <- subset(P5b, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 15&nCount_RNA > 1000 & nCount_RNA < 35000)
P5b@meta.data$state<-"before"
P5b@meta.data$sample<-"P5"
KA<-read.table("/data/AMLscRNA/scrublet/P5b_doublet.txt",header = T,sep = ",")
a<-P5b@meta.data
a$barcode<-rownames(P5b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P5b@meta.data<-c
P5b_sub=subset(P5b,subset=(predicted_doublets=="False"))
dim(P5b_sub)

P5a <- Read10X(data.dir = "/data/AMLscRNA/1320/outs/filtered_feature_bc_matrix")
P5a <- CreateSeuratObject(counts = P5a, project = "P5a", min.cells = 3, min.features = 200)
P5a[["percent.mt"]] <- PercentageFeatureSet(P5a, pattern = "^MT-")
VlnPlot(P5a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P5a <- subset(P5a, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 20&nCount_RNA > 1000 & nCount_RNA < 15000)
P5a@meta.data$state<-"after"
P5a@meta.data$sample<-"P5"
KA<-read.table("/data/AMLscRNA/scrublet/P5a_doublet.txt",header = T,sep = ",")
a<-P5a@meta.data
a$barcode<-rownames(P5a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P5a@meta.data<-c
P5a_sub=subset(P5a,subset=(predicted_doublets=="False"))
dim(P5a_sub)

P6b <- Read10X(data.dir = "/data/AMLscRNA/1363/outs/filtered_feature_bc_matrix")
P6b <- CreateSeuratObject(counts = P6b, project = "P6b", min.cells = 3, min.features = 200)
P6b[["percent.mt"]] <- PercentageFeatureSet(P6b, pattern = "^MT-")
VlnPlot(P6b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P6b <- subset(P6b, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 20&nCount_RNA > 1000 & nCount_RNA < 15000)
P6b@meta.data$state<-"before"
P6b@meta.data$sample<-"P6"
KA<-read.table("/data/AMLscRNA/scrublet/P6b_doublet.txt",header = T,sep = ",")
a<-P6b@meta.data
a$barcode<-rownames(P6b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P6b@meta.data<-c
P6b_sub=subset(P6b,subset=(predicted_doublets=="False"))
dim(P6b_sub)

P6a <- Read10X(data.dir = "/data/AMLscRNA/1503/outs/filtered_feature_bc_matrix")
P6a <- CreateSeuratObject(counts = P6a, project = "P6a", min.cells = 3, min.features = 200)
P6a[["percent.mt"]] <- PercentageFeatureSet(P6a, pattern = "^MT-")
VlnPlot(P6a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P6a <- subset(P6a, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20&nCount_RNA > 1000 & nCount_RNA < 15000)
P6a@meta.data$state<-"after"
P6a@meta.data$sample<-"P6"
KA<-read.table("/data/AMLscRNA/scrublet/P6a_doublet.txt",header = T,sep = ",")
a<-P6a@meta.data
a$barcode<-rownames(P6a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P6a@meta.data<-c
P6a_sub=subset(P6a,subset=(predicted_doublets=="False"))
dim(P6a_sub)

P7b <- Read10X(data.dir = "/data/AMLscRNA/1309/outs/filtered_feature_bc_matrix")
P7b <- CreateSeuratObject(counts = P7b, project = "P7b", min.cells = 3, min.features = 200)
P7b[["percent.mt"]] <- PercentageFeatureSet(P7b, pattern = "^MT-")
VlnPlot(P7b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P7b <- subset(P7b, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt <15&nCount_RNA > 1000 & nCount_RNA < 20000)
P7b@meta.data$state<-"before"
P7b@meta.data$sample<-"P7"
KA<-read.table("/data/AMLscRNA/scrublet/P7b_doublet.txt",header = T,sep = ",")
a<-P7b@meta.data
a$barcode<-rownames(P7b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P7b@meta.data<-c
P7b_sub=subset(P7b,subset=(predicted_doublets=="False"))
dim(P7b_sub)

P7a <- Read10X(data.dir = "/data/AMLscRNA/1449/outs/filtered_feature_bc_matrix")
P7a <- CreateSeuratObject(counts = P7a, project = "P7a", min.cells = 3, min.features = 200)
P7a[["percent.mt"]] <- PercentageFeatureSet(P7a, pattern = "^MT-")
VlnPlot(P7a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P7a <- subset(P7a, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 15&nCount_RNA >1000 & nCount_RNA < 20000)
P7a@meta.data$state<-"after"
P7a@meta.data$sample<-"P7"
KA<-read.table("/data/AMLscRNA/scrublet/P7a_doublet.txt",header = T,sep = ",")
a<-P7a@meta.data
a$barcode<-rownames(P7a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P7a@meta.data<-c
P7a_sub=subset(P7a,subset=(predicted_doublets=="False"))
dim(P7a_sub)

P8b <- Read10X(data.dir = "/data/AMLscRNA/YMF-B-1677/filtered_feature_bc_matrix")
P8b <- CreateSeuratObject(counts = P8b, project = "P8b", min.cells = 3, min.features = 200)
P8b[["percent.mt"]] <- PercentageFeatureSet(P8b, pattern = "^MT-")
VlnPlot(P8b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P8b <- subset(P8b, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20&nCount_RNA > 1000 & nCount_RNA < 25000)
P8b@meta.data$state<-"before"
P8b@meta.data$sample<-"P8"
KA<-read.table("/data/AMLscRNA/scrublet/P8b_doublet.txt",header = T,sep = ",")
a<-P8b@meta.data
a$barcode<-rownames(P8b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P8b@meta.data<-c
P8b_sub=subset(P8b,subset=(predicted_doublets=="False"))
dim(P8b_sub)

P8a <- Read10X(data.dir = "/data/AMLscRNA/YMF-A-1817/filtered_feature_bc_matrix")
P8a <- CreateSeuratObject(counts = P8a, project = "P8a", min.cells = 3, min.features = 200)
P8a[["percent.mt"]] <- PercentageFeatureSet(P8a, pattern = "^MT-")
VlnPlot(P8a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P8a <- subset(P8a, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20&nCount_RNA > 1000 & nCount_RNA < 10000)
P8a@meta.data$state<-"after"
P8a@meta.data$sample<-"P8"
KA<-read.table("/data/AMLscRNA/scrublet/P8a_doublet.txt",header = T,sep = ",")
a<-P8a@meta.data
a$barcode<-rownames(P8a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P8a@meta.data<-c
P8a_sub=subset(P8a,subset=(predicted_doublets=="False"))
dim(P8a_sub)

P9b <- Read10X(data.dir = "/data/AMLscRNA/HZW-B-2238/filtered_feature_bc_matrix")
P9b <- CreateSeuratObject(counts = P9b, project = "P9b", min.cells = 3, min.features = 200)
P9b[["percent.mt"]] <- PercentageFeatureSet(P9b, pattern = "^MT-")
VlnPlot(P9b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P9b <- subset(P9b, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20&nCount_RNA > 1000 & nCount_RNA < 30000)
P9b@meta.data$state<-"before"
P9b@meta.data$sample<-"P9"
KA<-read.table("/data/AMLscRNA/scrublet/P9b_doublet.txt",header = T,sep = ",")
a<-P9b@meta.data
a$barcode<-rownames(P9b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P9b@meta.data<-c
P9b_sub=subset(P9b,subset=(predicted_doublets=="False"))
dim(P9b_sub)

P9a <- Read10X(data.dir = "/data/AMLscRNA/HZW-A-2381/filtered_feature_bc_matrix")
P9a <- CreateSeuratObject(counts = P9a, project = "P9a", min.cells = 3, min.features = 200)
P9a[["percent.mt"]] <- PercentageFeatureSet(P9a, pattern = "^MT-")
VlnPlot(P9a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P9a <- subset(P9a, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 15&nCount_RNA > 1000 & nCount_RNA < 30000)
P9a@meta.data$state<-"after"
P9a@meta.data$sample<-"P9"
KA<-read.table("/data/AMLscRNA/scrublet/P9a_doublet.txt",header = T,sep = ",")
a<-P9a@meta.data
a$barcode<-rownames(P9a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P9a@meta.data<-c
P9a_sub=subset(P9a,subset=(predicted_doublets=="False"))
dim(P9a_sub)

P10b <- Read10X(data.dir = "/data/AMLscRNA/ZYM-B-2215/filtered_feature_bc_matrix")
P10b <- CreateSeuratObject(counts = P10b, project = "P10b", min.cells = 3, min.features = 200)
P10b[["percent.mt"]] <- PercentageFeatureSet(P10b, pattern = "^MT-")
VlnPlot(P10b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P10b <- subset(P10b, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 10&nCount_RNA > 1000 & nCount_RNA < 35000)
P10b@meta.data$state<-"before"
P10b@meta.data$sample<-"P10"
KA<-read.table("/data/AMLscRNA/scrublet/P10b_doublet.txt",header = T,sep = ",")
a<-P10b@meta.data
a$barcode<-rownames(P10b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P10b@meta.data<-c
P10b_sub=subset(P10b,subset=(predicted_doublets=="False"))
dim(P10b_sub)

P10a <- Read10X(data.dir = "/data/AMLscRNA/ZYM-A-2328/filtered_feature_bc_matrix")
P10a <- CreateSeuratObject(counts = P10a, project = "P10a", min.cells = 3, min.features = 200)
P10a[["percent.mt"]] <- PercentageFeatureSet(P10a, pattern = "^MT-")
VlnPlot(P10a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
P10a <- subset(P10a, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 10&nCount_RNA > 1000 & nCount_RNA < 25000)
P10a@meta.data$state<-"after"
P10a@meta.data$sample<-"P10"
KA<-read.table("/data/AMLscRNA/scrublet/P10a_doublet.txt",header = T,sep = ",")
a<-P10a@meta.data
a$barcode<-rownames(P10a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
P10a@meta.data<-c
P10a_sub=subset(P10a,subset=(predicted_doublets=="False"))
dim(P10a_sub)

N26a <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/26_1")
dim(N26a)
N26a <- CreateSeuratObject(counts = N26a, project = "N26a", min.cells = 3, min.features = 200)
N26a[["percent.mt"]] <- PercentageFeatureSet(N26a, pattern = "^MT-")
VlnPlot(N26a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N26a <- subset(N26a, subset = nFeature_RNA >200 & nFeature_RNA < 5000 & percent.mt < 25&nCount_RNA > 1000 & nCount_RNA < 30000)
N26a@meta.data$state<-"a"
N26a@meta.data$sample<-"N26"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/26_1_doublet.txt",header = T,sep = ",")
a<-N26a@meta.data
a$barcode<-rownames(N26a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N26a@meta.data<-c
N26a_sub=subset(N26a,subset=(predicted_doublets=="False"))
dim(N26a_sub)

N26b <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/26_2")
dim(N26b)
N26b <- CreateSeuratObject(counts = N26b, project = "N26b", min.cells = 3, min.features = 200)
N26b[["percent.mt"]] <- PercentageFeatureSet(N26b, pattern = "^MT-")
VlnPlot(N26b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N26b <- subset(N26b, subset = nFeature_RNA >200 & nFeature_RNA < 1500 & percent.mt < 25&nCount_RNA > 1000 & nCount_RNA < 6000)
N26b@meta.data$state<-"b"
N26b@meta.data$sample<-"N26"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/26_2_doublet.txt",header = T,sep = ",")
a<-N26b@meta.data
a$barcode<-rownames(N26b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N26b@meta.data<-c
N26b_sub=subset(N26b,subset=(predicted_doublets=="False"))
dim(N26b_sub)

N26c <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/26_3")
dim(N26c)
N26c <- CreateSeuratObject(counts = N26c, project = "N26c", min.cells = 3, min.features = 200)
N26c[["percent.mt"]] <- PercentageFeatureSet(N26c, pattern = "^MT-")
VlnPlot(N26c, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N26c <- subset(N26c, subset = nFeature_RNA >200 & nFeature_RNA < 6000 & percent.mt < 25&nCount_RNA > 1000 & nCount_RNA < 25000)
N26c@meta.data$state<-"c"
N26c@meta.data$sample<-"N26"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/26_3_doublet.txt",header = T,sep = ",")
a<-N26c@meta.data
a$barcode<-rownames(N26c@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N26c@meta.data<-c
N26c_sub=subset(N26c,subset=(predicted_doublets=="False"))
dim(N26c_sub)

N26 <- merge(N26a_sub, y = c(N26b_sub,N26c_sub), add.cell.ids = c("N26a","N26b","N26c"), project = "N26")
SaveLoom(N26, filename = "/data/AMLscRNA/N26.loom", verbose = FALSE)

N27a <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/27_1")
dim(N27a)
N27a <- CreateSeuratObject(counts = N27a, project = "N27a", min.cells = 3, min.features = 200)
N27a[["percent.mt"]] <- PercentageFeatureSet(N27a, pattern = "^MT-")
VlnPlot(N27a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N27a <- subset(N27a, subset = nFeature_RNA >200 & nFeature_RNA < 6000 & percent.mt < 25&nCount_RNA > 1000 & nCount_RNA < 40000)
N27a@meta.data$state<-"a"
N27a@meta.data$sample<-"N27"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/27_1_doublet.txt",header = T,sep = ",")
a<-N27a@meta.data
a$barcode<-rownames(N27a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N27a@meta.data<-c
N27a_sub=subset(N27a,subset=(predicted_doublets=="False"))
dim(N27a_sub)

N27b <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/27_2")
dim(N27b)
N27b <- CreateSeuratObject(counts = N27b, project = "N27b", min.cells = 3, min.features = 200)
N27b[["percent.mt"]] <- PercentageFeatureSet(N27b, pattern = "^MT-")
VlnPlot(N27b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N27b <- subset(N27b, subset = nFeature_RNA >200 & nFeature_RNA < 5000 & percent.mt < 25&nCount_RNA > 1000 & nCount_RNA < 30000)
N27b@meta.data$state<-"b"
N27b@meta.data$sample<-"N27"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/27_2_doublet.txt",header = T,sep = ",")
a<-N27b@meta.data
a$barcode<-rownames(N27b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N27b@meta.data<-c
N27b_sub=subset(N27b,subset=(predicted_doublets=="False"))
dim(N27b_sub)

N27c <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/27_3")
dim(N27c)
N27c <- CreateSeuratObject(counts = N27c, project = "N27c", min.cells = 3, min.features = 200)
N27c[["percent.mt"]] <- PercentageFeatureSet(N27c, pattern = "^MT-")
VlnPlot(N27c, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N27c <- subset(N27c, subset = nFeature_RNA >200 & nFeature_RNA < 6000 & percent.mt < 25&nCount_RNA > 1000 & nCount_RNA < 25000)
N27c@meta.data$state<-"c"
N27c@meta.data$sample<-"N27"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/27_3_doublet.txt",header = T,sep = ",")
a<-N27c@meta.data
a$barcode<-rownames(N27c@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N27c@meta.data<-c
N27c_sub=subset(N27c,subset=(predicted_doublets=="False"))
dim(N27c_sub)

N27 <- merge(N27a_sub, y = c(N27b_sub,N27c_sub), add.cell.ids = c("N27a","N27b","N27c"), project = "N27")
SaveLoom(N27, filename = "/data/AMLscRNA/N27.loom", verbose = FALSE)

N28a <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/28_1")
dim(N28a)
N28a <- CreateSeuratObject(counts = N28a, project = "N28a", min.cells = 3, min.features = 200)
N28a[["percent.mt"]] <- PercentageFeatureSet(N28a, pattern = "^MT-")
VlnPlot(N28a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N28a <- subset(N28a, subset = nFeature_RNA >200 & nFeature_RNA < 6000 & percent.mt < 10&nCount_RNA > 1000 & nCount_RNA < 40000)
N28a@meta.data$state<-"a"
N28a@meta.data$sample<-"N28"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/28_1_doublet.txt",header = T,sep = ",")
a<-N28a@meta.data
a$barcode<-rownames(N28a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N28a@meta.data<-c
N28a_sub=subset(N28a,subset=(predicted_doublets=="False"))
dim(N28a_sub)

N28b <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/28_2")
dim(N28b)
N28b <- CreateSeuratObject(counts = N28b, project = "N28b", min.cells = 3, min.features = 200)
N28b[["percent.mt"]] <- PercentageFeatureSet(N28b, pattern = "^MT-")
VlnPlot(N28b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N28b <- subset(N28b, subset = nFeature_RNA >200 & nFeature_RNA < 4500 & percent.mt < 10&nCount_RNA > 1000 & nCount_RNA < 25000)
N28b@meta.data$state<-"b"
N28b@meta.data$sample<-"N28"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/28_2_doublet.txt",header = T,sep = ",")
a<-N28b@meta.data
a$barcode<-rownames(N28b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N28b@meta.data<-c
N28b_sub=subset(N28b,subset=(predicted_doublets=="False"))
dim(N28b_sub)

N28c <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/28_3")
dim(N28c)
N28c <- CreateSeuratObject(counts = N28c, project = "N28c", min.cells = 3, min.features = 200)
N28c[["percent.mt"]] <- PercentageFeatureSet(N28c, pattern = "^MT-")
VlnPlot(N28c, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N28c <- subset(N28c, subset = nFeature_RNA >200 & nFeature_RNA < 4000 & percent.mt < 20&nCount_RNA > 1000 & nCount_RNA <20000)
N28c@meta.data$state<-"c"
N28c@meta.data$sample<-"N28"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/28_3_doublet.txt",header = T,sep = ",")
a<-N28c@meta.data
a$barcode<-rownames(N28c@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N28c@meta.data<-c
N28c_sub=subset(N28c,subset=(predicted_doublets=="False"))
dim(N28c_sub)

N28 <- merge(N28a_sub, y = c(N28b_sub,N28c_sub), add.cell.ids = c("N28a","N28b","N28c"), project = "N28")
SaveLoom(N28, filename = "/data/AMLscRNA/N28.loom", verbose = FALSE)

N29a <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/29_1")
dim(N29a)
N29a <- CreateSeuratObject(counts = N29a, project = "N29a", min.cells = 3, min.features = 200)
N29a[["percent.mt"]] <- PercentageFeatureSet(N29a, pattern = "^MT-")
VlnPlot(N29a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N29a <- subset(N29a, subset = nFeature_RNA >200 & nFeature_RNA < 7000 & percent.mt < 15&nCount_RNA > 1000 & nCount_RNA < 50000)
N29a@meta.data$state<-"a"
N29a@meta.data$sample<-"N29"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/29_1_doublet.txt",header = T,sep = ",")
a<-N29a@meta.data
a$barcode<-rownames(N29a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N29a@meta.data<-c
N29a_sub=subset(N29a,subset=(predicted_doublets=="False"))
dim(N29a_sub)

N29b <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/29_2")
dim(N29b)
N29b <- CreateSeuratObject(counts = N29b, project = "N29b", min.cells = 3, min.features = 200)
N29b[["percent.mt"]] <- PercentageFeatureSet(N29b, pattern = "^MT-")
VlnPlot(N29b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N29b <- subset(N29b, subset = nFeature_RNA >200 & nFeature_RNA < 6000 & percent.mt < 7.5&nCount_RNA > 1000 & nCount_RNA < 30000)
N29b@meta.data$state<-"b"
N29b@meta.data$sample<-"N29"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/29_2_doublet.txt",header = T,sep = ",")
a<-N29b@meta.data
a$barcode<-rownames(N29b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N29b@meta.data<-c
N29b_sub=subset(N29b,subset=(predicted_doublets=="False"))
dim(N29b_sub)

N29c <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/29_3")
dim(N29c)
N29c <- CreateSeuratObject(counts = N29c, project = "N29c", min.cells = 3, min.features = 200)
N29c[["percent.mt"]] <- PercentageFeatureSet(N29c, pattern = "^MT-")
VlnPlot(N29c, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N29c <- subset(N29c, subset = nFeature_RNA >200 & nFeature_RNA <6000 & percent.mt <10&nCount_RNA > 1000 & nCount_RNA < 30000)
N29c@meta.data$state<-"c"
N29c@meta.data$sample<-"N29"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/29_3_doublet.txt",header = T,sep = ",")
a<-N29c@meta.data
a$barcode<-rownames(N29c@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N29c@meta.data<-c
N29c_sub=subset(N29c,subset=(predicted_doublets=="False"))
dim(N29c_sub)

N29 <- merge(N29a_sub, y = c(N29b_sub,N29c_sub), add.cell.ids = c("N29a","N29b","N29c"), project = "N29")
SaveLoom(N29, filename = "/data/AMLscRNA/N29.loom", verbose = FALSE)

N30a <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/30_1")
dim(N30a)
N30a <- CreateSeuratObject(counts = N30a, project = "N30a", min.cells = 3, min.features = 200)
N30a[["percent.mt"]] <- PercentageFeatureSet(N30a, pattern = "^MT-")
VlnPlot(N30a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N30a <- subset(N30a, subset = nFeature_RNA >200 & nFeature_RNA <6000 & percent.mt < 15&nCount_RNA > 1000 & nCount_RNA < 40000)
N30a@meta.data$state<-"a"
N30a@meta.data$sample<-"N30"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/30_1_doublet.txt",header = T,sep = ",")
a<-N30a@meta.data
a$barcode<-rownames(N30a@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N30a@meta.data<-c
N30a_sub=subset(N30a,subset=(predicted_doublets=="False"))
dim(N30a_sub)

N30b <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/30_2")
dim(N30b)
N30b <- CreateSeuratObject(counts = N30b, project = "N30b", min.cells = 3, min.features = 200)
N30b[["percent.mt"]] <- PercentageFeatureSet(N30b, pattern = "^MT-")
VlnPlot(N30b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N30b <- subset(N30b, subset = nFeature_RNA >200 & nFeature_RNA < 5000 & percent.mt < 25&nCount_RNA > 1000 & nCount_RNA < 30000)
N30b@meta.data$state<-"b"
N30b@meta.data$sample<-"N30"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/30_2_doublet.txt",header = T,sep = ",")
a<-N30b@meta.data
a$barcode<-rownames(N30b@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N30b@meta.data<-c
N30b_sub=subset(N30b,subset=(predicted_doublets=="False"))
dim(N30b_sub)

N30c <- Read10X(data.dir = "/data/AMLscRNA/BM/blood/30_3")
dim(N30c)
N30c <- CreateSeuratObject(counts = N30c, project = "N30c", min.cells = 3, min.features = 200)
N30c[["percent.mt"]] <- PercentageFeatureSet(N30c, pattern = "^MT-")
VlnPlot(N30c, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
N30c <- subset(N30c, subset = nFeature_RNA >200 & nFeature_RNA < 4000 & percent.mt <15&nCount_RNA > 1000 & nCount_RNA < 20000)
N30c@meta.data$state<-"c"
N30c@meta.data$sample<-"N30"
KA<-read.table("/data/AMLscRNA/BM/blood/scrublet/30_3_doublet.txt",header = T,sep = ",")
a<-N30c@meta.data
a$barcode<-rownames(N30c@meta.data)
c<-merge(a,KA,by="barcode")
rownames(c)=c[,1] 
c=c[,-1]
N30c@meta.data<-c
N30c_sub=subset(N30c,subset=(predicted_doublets=="False"))
dim(N30c_sub)

AML <- merge(P1b_sub, y = c(P1a_sub,P2b_sub,P2a_sub,P3b_sub,P3a_sub,P4b_sub,P4a_sub,P5b_sub,P5a_sub,P6b_sub,P6a_sub,P7b_sub,P7a_sub,P8b_sub,P8a_sub,P9b_sub,P9a_sub,P10b_sub,P10a_sub), add.cell.ids = c("P1b","P1t","P2b","P2t","P3b","P3t","P4b","P4t","P5b","P5t","P6b","P6t","P7b","P7t","P8b","P8t","P9b","P9t","P10b","P10t"), project = "AML")
healt <- merge(N26a_sub, y = c(N26b_sub,N26c_sub,N27a_sub,N27b_sub,N27c_sub,N28a_sub,N28b_sub,N28c_sub,N29a_sub,N29b_sub,N29c_sub,N30a_sub,N30b_sub,N30c_sub), add.cell.ids = c("N26a","N26b","N26c","N27a","N27b","N27c","N28a","N28b","N28c","N29a","N29b","N29c","N30a","N30b","N30c"), project = "health")

library(SeuratDisk)
SaveLoom(AML, filename = "/data/AMLscRNA/AML_10p.loom", verbose = FALSE)
SaveLoom(healt, filename = "/data/AMLscRNA/BM/blood/blood5n.loom", verbose = FALSE)
