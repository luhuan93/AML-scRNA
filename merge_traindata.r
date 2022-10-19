setwd("/data/yifan_AMLdata/dataset/Cell2019/raw/")
library(Seurat)
AML1012_D0<- read.table(gzfile("GSM3587923_AML1012_D0.dem.txt.gz"),header = T,row.names="Gene")
AML1012_D0<- CreateSeuratObject(counts = AML1012_D0, min.cells = 3, project = "AML1012_D0")
AML1012_D0_anno<- read.table(gzfile("GSM3587924_AML1012_D0.anno.txt.gz"),header = T,sep="\t")
AML1012_D0_anno$Cell<- gsub("-",".",AML1012_D0_anno$Cell)
a<-AML1012_D0@meta.data
a$Cell<-rownames(AML1012_D0@meta.data)
c<-merge(a,AML1012_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML1012_D0@meta.data<-c
AML210A_D0<- read.table(gzfile("GSM3587925_AML210A_D0.dem.txt.gz"),header = T,row.names="Gene")
AML210A_D0<- CreateSeuratObject(counts = AML210A_D0, min.cells = 3, project = "AML210A_D0")
AML210A_D0_anno<- read.table(gzfile("GSM3587926_AML210A_D0.anno.txt.gz"),header = T,sep="\t")
AML210A_D0_anno$Cell<- gsub("-",".",AML210A_D0_anno$Cell)
a<-AML210A_D0@meta.data
a$Cell<-rownames(AML210A_D0@meta.data)
c<-merge(a,AML210A_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML210A_D0@meta.data<-c
AML314_D0<- read.table(gzfile("GSM3587927_AML314_D0.dem.txt.gz"),header = T,row.names="Gene")
AML314_D0<- CreateSeuratObject(counts = AML314_D0, min.cells = 3, project = "AML314_D0")
AML314_D0_anno<- read.table(gzfile("GSM3587928_AML314_D0.anno.txt.gz"),header = T,sep="\t")
AML314_D0_anno$Cell<- gsub("-",".",AML314_D0_anno$Cell)
a<-AML314_D0@meta.data
a$Cell<-rownames(AML314_D0@meta.data)
c<-merge(a,AML314_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML314_D0@meta.data<-c
AML314_D31<- read.table(gzfile("GSM3587929_AML314_D31.dem.txt.gz"),header = T,row.names="Gene")
AML314_D31<- CreateSeuratObject(counts = AML314_D31, min.cells = 3, project = "AML314_D31")
AML314_D31_anno<- read.table(gzfile("GSM3587930_AML314_D31.anno.txt.gz"),header = T,sep="\t")
AML314_D31_anno$Cell<- gsub("-",".",AML314_D31_anno$Cell)
a<-AML314_D31@meta.data
a$Cell<-rownames(AML314_D31@meta.data)
c<-merge(a,AML314_D31_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML314_D31@meta.data<-c
AML328_D0<- read.table(gzfile("GSM3587931_AML328_D0.dem.txt.gz"),header = T,row.names="Gene")
AML328_D0<- CreateSeuratObject(counts = AML328_D0, min.cells = 3, project = "AML328_D0")
AML328_D0_anno<- read.table(gzfile("GSM3587932_AML328_D0.anno.txt.gz"),header = T,sep="\t")
AML328_D0_anno$Cell<- gsub("-",".",AML328_D0_anno$Cell)
a<-AML328_D0@meta.data
a$Cell<-rownames(AML328_D0@meta.data)
c<-merge(a,AML328_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML328_D0@meta.data<-c
AML328_D113<- read.table(gzfile("GSM3587933_AML328_D113.dem.txt.gz"),header = T,row.names="Gene")
AML328_D113<- CreateSeuratObject(counts = AML328_D113, min.cells = 3, project = "AML328_D113")
AML328_D113_anno<- read.table(gzfile("GSM3587934_AML328_D113.anno.txt.gz"),header = T,sep="\t")
AML328_D113_anno$Cell<- gsub("-",".",AML328_D113_anno$Cell)
a<-AML328_D113@meta.data
a$Cell<-rownames(AML328_D113@meta.data)
c<-merge(a,AML328_D113_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML328_D113@meta.data<-c
AML328_D171<- read.table(gzfile("GSM3587935_AML328_D171.dem.txt.gz"),header = T,row.names="Gene")
AML328_D171<- CreateSeuratObject(counts = AML328_D171, min.cells = 3, project = "AML328_D171")
AML328_D171_anno<- read.table(gzfile("GSM3587936_AML328_D171.anno.txt.gz"),header = T,sep="\t")
AML328_D171_anno$Cell<- gsub("-",".",AML328_D171_anno$Cell)
a<-AML328_D171@meta.data
a$Cell<-rownames(AML328_D171@meta.data)
c<-merge(a,AML328_D171_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML328_D171@meta.data<-c
AML328_D29<- read.table(gzfile("GSM3587937_AML328_D29.dem.txt.gz"),header = T,row.names="Gene")
AML328_D29<- CreateSeuratObject(counts = AML328_D29, min.cells = 3, project = "AML328_D29")
AML328_D29_anno<- read.table(gzfile("GSM3587938_AML328_D29.anno.txt.gz"),header = T,sep="\t")
AML328_D29_anno$Cell<- gsub("-",".",AML328_D29_anno$Cell)
a<-AML328_D29@meta.data
a$Cell<-rownames(AML328_D29@meta.data)
c<-merge(a,AML328_D29_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML328_D29@meta.data<-c
AML329_D0<- read.table(gzfile("GSM3587940_AML329_D0.dem.txt.gz"),header = T,row.names="Gene")
AML329_D0<- CreateSeuratObject(counts = AML329_D0, min.cells = 3, project = "AML329_D0")
AML329_D0_anno<- read.table(gzfile("GSM3587941_AML329_D0.anno.txt.gz"),header = T,sep="\t")
AML329_D0_anno$Cell<- gsub("-",".",AML329_D0_anno$Cell)
a<-AML329_D0@meta.data
a$Cell<-rownames(AML329_D0@meta.data)
c<-merge(a,AML329_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML329_D0@meta.data<-c
AML329_D20<- read.table(gzfile("GSM3587942_AML329_D20.dem.txt.gz"),header = T,row.names="Gene")
AML329_D20<- CreateSeuratObject(counts = AML329_D20, min.cells = 3, project = "AML329_D20")
AML329_D20_anno<- read.table(gzfile("GSM3587943_AML329_D20.anno.txt.gz"),header = T,sep="\t")
AML329_D20_anno$Cell<- gsub("-",".",AML329_D20_anno$Cell)
a<-AML329_D20@meta.data
a$Cell<-rownames(AML329_D20@meta.data)
c<-merge(a,AML329_D20_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML329_D20@meta.data<-c
AML329_D37<- read.table(gzfile("GSM3587944_AML329_D37.dem.txt.gz"),header = T,row.names="Gene")
AML329_D37<- CreateSeuratObject(counts = AML329_D37, min.cells = 3, project = "AML329_D37")
AML329_D37_anno<- read.table(gzfile("GSM3587945_AML329_D37.anno.txt.gz"),header = T,sep="\t")
AML329_D37_anno$Cell<- gsub("-",".",AML329_D37_anno$Cell)
a<-AML329_D37@meta.data
a$Cell<-rownames(AML329_D37@meta.data)
c<-merge(a,AML329_D37_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML329_D37@meta.data<-c
AML371_D0<- read.table(gzfile("GSM3587946_AML371_D0.dem.txt.gz"),header = T,row.names="Gene")
AML371_D0<- CreateSeuratObject(counts = AML371_D0, min.cells = 3, project = "AML371_D0")
AML371_D0_anno<- read.table(gzfile("GSM3587947_AML371_D0.anno.txt.gz"),header = T,sep="\t")
AML371_D0_anno$Cell<- gsub("-",".",AML371_D0_anno$Cell)
a<-AML371_D0@meta.data
a$Cell<-rownames(AML371_D0@meta.data)
c<-merge(a,AML371_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML371_D0@meta.data<-c
AML371_D34<- read.table(gzfile("GSM3587948_AML371_D34.dem.txt.gz"),header = T,row.names="Gene")
AML371_D34<- CreateSeuratObject(counts = AML371_D34, min.cells = 3, project = "AML371_D34")
AML371_D34_anno<- read.table(gzfile("GSM3587949_AML371_D34.anno.txt.gz"),header = T,sep="\t")
AML371_D34_anno$Cell<- gsub("-",".",AML371_D34_anno$Cell)
a<-AML371_D34@meta.data
a$Cell<-rownames(AML371_D34@meta.data)
c<-merge(a,AML371_D34_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML371_D34@meta.data<-c
AML419A_D0<- read.table(gzfile("GSM3587950_AML419A_D0.dem.txt.gz"),header = T,row.names="Gene")
AML419A_D0<- CreateSeuratObject(counts = AML419A_D0, min.cells = 3, project = "AML419A_D0")
AML419A_D0_anno<- read.table(gzfile("GSM3587951_AML419A_D0.anno.txt.gz"),header = T,sep="\t")
AML419A_D0_anno$Cell<- gsub("-",".",AML419A_D0_anno$Cell)
a<-AML419A_D0@meta.data
a$Cell<-rownames(AML419A_D0@meta.data)
c<-merge(a,AML419A_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML419A_D0@meta.data<-c
AML420B_D0<- read.table(gzfile("GSM3587953_AML420B_D0.dem.txt.gz"),header = T,row.names="Gene")
AML420B_D0<- CreateSeuratObject(counts = AML420B_D0, min.cells = 3, project = "AML420B_D0")
AML420B_D0_anno<- read.table(gzfile("GSM3587954_AML420B_D0.anno.txt.gz"),header = T,sep="\t")
AML420B_D0_anno$Cell<- gsub("-",".",AML420B_D0_anno$Cell)
a<-AML420B_D0@meta.data
a$Cell<-rownames(AML420B_D0@meta.data)
c<-merge(a,AML420B_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML420B_D0@meta.data<-c
AML420B_D14<- read.table(gzfile("GSM3587955_AML420B_D14.dem.txt.gz"),header = T,row.names="Gene")
AML420B_D14<- CreateSeuratObject(counts = AML420B_D14, min.cells = 3, project = "AML420B_D14")
AML420B_D14_anno<- read.table(gzfile("GSM3587956_AML420B_D14.anno.txt.gz"),header = T,sep="\t")
AML420B_D14_anno$Cell<- gsub("-",".",AML420B_D14_anno$Cell)
a<-AML420B_D14@meta.data
a$Cell<-rownames(AML420B_D14@meta.data)
c<-merge(a,AML420B_D14_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML420B_D14@meta.data<-c
AML420B_D35<- read.table(gzfile("GSM3587957_AML420B_D35.dem.txt.gz"),header = T,row.names="Gene")
AML420B_D35<- CreateSeuratObject(counts = AML420B_D35, min.cells = 3, project = "AML420B_D35")
AML420B_D35_anno<- read.table(gzfile("GSM3587958_AML420B_D35.anno.txt.gz"),header = T,sep="\t")
AML420B_D35_anno$Cell<- gsub("-",".",AML420B_D35_anno$Cell)
a<-AML420B_D35@meta.data
a$Cell<-rownames(AML420B_D35@meta.data)
c<-merge(a,AML420B_D35_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML420B_D35@meta.data<-c
AML475_D0<- read.table(gzfile("GSM3587959_AML475_D0.dem.txt.gz"),header = T,row.names="Gene")
AML475_D0<- CreateSeuratObject(counts = AML475_D0, min.cells = 3, project = "AML475_D0")
AML475_D0_anno<- read.table(gzfile("GSM3587960_AML475_D0.anno.txt.gz"),header = T,sep="\t")
AML475_D0_anno$Cell<- gsub("-",".",AML475_D0_anno$Cell)
a<-AML475_D0@meta.data
a$Cell<-rownames(AML475_D0@meta.data)
c<-merge(a,AML475_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML475_D0@meta.data<-c
AML475_D29<- read.table(gzfile("GSM3587961_AML475_D29.dem.txt.gz"),header = T,row.names="Gene")
AML475_D29<- CreateSeuratObject(counts = AML475_D29, min.cells = 3, project = "AML475_D29")
AML475_D29_anno<- read.table(gzfile("GSM3587962_AML475_D29.anno.txt.gz"),header = T,sep="\t")
AML475_D29_anno$Cell<- gsub("-",".",AML475_D29_anno$Cell)
a<-AML475_D29@meta.data
a$Cell<-rownames(AML475_D29@meta.data)
c<-merge(a,AML475_D29_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML475_D29@meta.data<-c
AML556_D0<- read.table(gzfile("GSM3587963_AML556_D0.dem.txt.gz"),header = T,row.names="Gene")
AML556_D0<- CreateSeuratObject(counts = AML556_D0, min.cells = 3, project = "AML556_D0")
AML556_D0_anno<- read.table(gzfile("GSM3587964_AML556_D0.anno.txt.gz"),header = T,sep="\t")
AML556_D0_anno$Cell<- gsub("-",".",AML556_D0_anno$Cell)
a<-AML556_D0@meta.data
a$Cell<-rownames(AML556_D0@meta.data)
c<-merge(a,AML556_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML556_D0@meta.data<-c
AML556_D15<- read.table(gzfile("GSM3587965_AML556_D15.dem.txt.gz"),header = T,row.names="Gene")
AML556_D15<- CreateSeuratObject(counts = AML556_D15, min.cells = 3, project = "AML556_D15")
AML556_D15_anno<- read.table(gzfile("GSM3587966_AML556_D15.anno.txt.gz"),header = T,sep="\t")
AML556_D15_anno$Cell<- gsub("-",".",AML556_D15_anno$Cell)
a<-AML556_D15@meta.data
a$Cell<-rownames(AML556_D15@meta.data)
c<-merge(a,AML556_D15_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML556_D15@meta.data<-c
AML556_D31<- read.table(gzfile("GSM3587967_AML556_D31.dem.txt.gz"),header = T,row.names="Gene")
AML556_D31<- CreateSeuratObject(counts = AML556_D31, min.cells = 3, project = "AML556_D31")
AML556_D31_anno<- read.table(gzfile("GSM3587968_AML556_D31.anno.txt.gz"),header = T,sep="\t")
AML556_D31_anno$Cell<- gsub("-",".",AML556_D31_anno$Cell)
a<-AML556_D31@meta.data
a$Cell<-rownames(AML556_D31@meta.data)
c<-merge(a,AML556_D31_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML556_D31@meta.data<-c
AML707B_D0<- read.table(gzfile("GSM3587969_AML707B_D0.dem.txt.gz"),header = T,row.names="Gene")
AML707B_D0<- CreateSeuratObject(counts = AML707B_D0, min.cells = 3, project = "AML707B_D0")
AML707B_D0_anno<- read.table(gzfile("GSM3587970_AML707B_D0.anno.txt.gz"),header = T,sep="\t")
AML707B_D0_anno$Cell<- gsub("-",".",AML707B_D0_anno$Cell)
a<-AML707B_D0@meta.data
a$Cell<-rownames(AML707B_D0@meta.data)
c<-merge(a,AML707B_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML707B_D0@meta.data<-c
AML707B_D113<- read.table(gzfile("GSM3587971_AML707B_D113.dem.txt.gz"),header = T,row.names="Gene")
AML707B_D113<- CreateSeuratObject(counts = AML707B_D113, min.cells = 3, project = "AML707B_D113")
AML707B_D113_anno<- read.table(gzfile("GSM3587972_AML707B_D113.anno.txt.gz"),header = T,sep="\t")
AML707B_D113_anno$Cell<- gsub("-",".",AML707B_D113_anno$Cell)
a<-AML707B_D113@meta.data
a$Cell<-rownames(AML707B_D113@meta.data)
c<-merge(a,AML707B_D113_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML707B_D113@meta.data<-c
AML707B_D18<- read.table(gzfile("GSM3587973_AML707B_D18.dem.txt.gz"),header = T,row.names="Gene")
AML707B_D18<- CreateSeuratObject(counts = AML707B_D18, min.cells = 3, project = "AML707B_D18")
AML707B_D18_anno<- read.table(gzfile("GSM3587974_AML707B_D18.anno.txt.gz"),header = T,sep="\t")
AML707B_D18_anno$Cell<- gsub("-",".",AML707B_D18_anno$Cell)
a<-AML707B_D18@meta.data
a$Cell<-rownames(AML707B_D18@meta.data)
c<-merge(a,AML707B_D18_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML707B_D18@meta.data<-c
AML707B_D41<- read.table(gzfile("GSM3587975_AML707B_D41.dem.txt.gz"),header = T,row.names="Gene")
AML707B_D41<- CreateSeuratObject(counts = AML707B_D41, min.cells = 3, project = "AML707B_D41")
AML707B_D41_anno<- read.table(gzfile("GSM3587976_AML707B_D41.anno.txt.gz"),header = T,sep="\t")
AML707B_D41_anno$Cell<- gsub("-",".",AML707B_D41_anno$Cell)
a<-AML707B_D41@meta.data
a$Cell<-rownames(AML707B_D41@meta.data)
c<-merge(a,AML707B_D41_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML707B_D41@meta.data<-c
AML707B_D97<- read.table(gzfile("GSM3587977_AML707B_D97.dem.txt.gz"),header = T,row.names="Gene")
AML707B_D97<- CreateSeuratObject(counts = AML707B_D97, min.cells = 3, project = "AML707B_D97")
AML707B_D97_anno<- read.table(gzfile("GSM3587978_AML707B_D97.anno.txt.gz"),header = T,sep="\t")
AML707B_D97_anno$Cell<- gsub("-",".",AML707B_D97_anno$Cell)
a<-AML707B_D97@meta.data
a$Cell<-rownames(AML707B_D97@meta.data)
c<-merge(a,AML707B_D97_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML707B_D97@meta.data<-c
AML722B_D0<- read.table(gzfile("GSM3587980_AML722B_D0.dem.txt.gz"),header = T,row.names="Gene")
AML722B_D0<- CreateSeuratObject(counts = AML722B_D0, min.cells = 3, project = "AML722B_D0")
AML722B_D0_anno<- read.table(gzfile("GSM3587981_AML722B_D0.anno.txt.gz"),header = T,sep="\t")
AML722B_D0_anno$Cell<- gsub("-",".",AML722B_D0_anno$Cell)
a<-AML722B_D0@meta.data
a$Cell<-rownames(AML722B_D0@meta.data)
c<-merge(a,AML722B_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML722B_D0@meta.data<-c
AML722B_D49<- read.table(gzfile("GSM3587982_AML722B_D49.dem.txt.gz"),header = T,row.names="Gene")
AML722B_D49<- CreateSeuratObject(counts = AML722B_D49, min.cells = 3, project = "AML722B_D49")
AML722B_D49_anno<- read.table(gzfile("GSM3587983_AML722B_D49.anno.txt.gz"),header = T,sep="\t")
AML722B_D49_anno$Cell<- gsub("-",".",AML722B_D49_anno$Cell)
a<-AML722B_D49@meta.data
a$Cell<-rownames(AML722B_D49@meta.data)
c<-merge(a,AML722B_D49_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML722B_D49@meta.data<-c
AML870_D0<- read.table(gzfile("GSM3587984_AML870_D0.dem.txt.gz"),header = T,row.names="Gene")
AML870_D0<- CreateSeuratObject(counts = AML870_D0, min.cells = 3, project = "AML870_D0")
AML870_D0_anno<- read.table(gzfile("GSM3587985_AML870_D0.anno.txt.gz"),header = T,sep="\t")
AML870_D0_anno$Cell<- gsub("-",".",AML870_D0_anno$Cell)
a<-AML870_D0@meta.data
a$Cell<-rownames(AML870_D0@meta.data)
c<-merge(a,AML870_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML870_D0@meta.data<-c
AML870_D14<- read.table(gzfile("GSM3587986_AML870_D14.dem.txt.gz"),header = T,row.names="Gene")
AML870_D14<- CreateSeuratObject(counts = AML870_D14, min.cells = 3, project = "AML870_D14")
AML870_D14_anno<- read.table(gzfile("GSM3587987_AML870_D14.anno.txt.gz"),header = T,sep="\t")
AML870_D14_anno$Cell<- gsub("-",".",AML870_D14_anno$Cell)
a<-AML870_D14@meta.data
a$Cell<-rownames(AML870_D14@meta.data)
c<-merge(a,AML870_D14_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML870_D14@meta.data<-c
AML916_D0<- read.table(gzfile("GSM3587988_AML916_D0.dem.txt.gz"),header = T,row.names="Gene")
AML916_D0<- CreateSeuratObject(counts = AML916_D0, min.cells = 3, project = "AML916_D0")
AML916_D0_anno<- read.table(gzfile("GSM3587989_AML916_D0.anno.txt.gz"),header = T,sep="\t")
AML916_D0_anno$Cell<- gsub("-",".",AML916_D0_anno$Cell)
a<-AML916_D0@meta.data
a$Cell<-rownames(AML916_D0@meta.data)
c<-merge(a,AML916_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML916_D0@meta.data<-c
AML921A_D0<- read.table(gzfile("GSM3587990_AML921A_D0.dem.txt.gz"),header = T,row.names="Gene")
AML921A_D0<- CreateSeuratObject(counts = AML921A_D0, min.cells = 3, project = "AML921A_D0")
AML921A_D0_anno<- read.table(gzfile("GSM3587991_AML921A_D0.anno.txt.gz"),header = T,sep="\t")
AML921A_D0_anno$Cell<- gsub("-",".",AML921A_D0_anno$Cell)
a<-AML921A_D0@meta.data
a$Cell<-rownames(AML921A_D0@meta.data)
c<-merge(a,AML921A_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML921A_D0@meta.data<-c
AML997_D0<- read.table(gzfile("GSM3587992_AML997_D0.dem.txt.gz"),header = T,row.names="Gene")
AML997_D0<- CreateSeuratObject(counts = AML997_D0, min.cells = 3, project = "AML997_D0")
AML997_D0_anno<- read.table(gzfile("GSM3587993_AML997_D0.anno.txt.gz"),header = T,sep="\t")
AML997_D0_anno$Cell<- gsub("-",".",AML997_D0_anno$Cell)
a<-AML997_D0@meta.data
a$Cell<-rownames(AML997_D0@meta.data)
c<-merge(a,AML997_D0_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML997_D0@meta.data<-c
AML997_D35<- read.table(gzfile("GSM3587994_AML997_D35.dem.txt.gz"),header = T,row.names="Gene")
AML997_D35<- CreateSeuratObject(counts = AML997_D35, min.cells = 3, project = "AML997_D35")
AML997_D35_anno<- read.table(gzfile("GSM3587995_AML997_D35.anno.txt.gz"),header = T,sep="\t")
AML997_D35_anno$Cell<- gsub("-",".",AML997_D35_anno$Cell)
a<-AML997_D35@meta.data
a$Cell<-rownames(AML997_D35@meta.data)
c<-merge(a,AML997_D35_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
AML997_D35@meta.data<-c
BM1<- read.table(gzfile("GSM3587996_BM1.dem.txt.gz"),header = T,row.names="Gene")
BM1<- CreateSeuratObject(counts = BM1, min.cells = 3, project = "BM1")
BM1_anno<- read.table(gzfile("GSM3587996_BM1.anno.txt.gz"),header = T,sep="\t")
BM1_anno$Cell<- gsub("-",".",BM1_anno$Cell)
a<-BM1@meta.data
a$Cell<-rownames(BM1@meta.data)
c<-merge(a,BM1_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
BM1@meta.data<-c
BM2<- read.table(gzfile("GSM3587997_BM2.dem.txt.gz"),header = T,row.names="Gene")
BM2<- CreateSeuratObject(counts = BM2, min.cells = 3, project = "BM2")
BM2_anno<- read.table(gzfile("GSM3587997_BM2.anno.txt.gz"),header = T,sep="\t")
BM2_anno$Cell<- gsub("-",".",BM2_anno$Cell)
a<-BM2@meta.data
a$Cell<-rownames(BM2@meta.data)
c<-merge(a,BM2_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
BM2@meta.data<-c
BM3<- read.table(gzfile("GSM3587998_BM3.dem.txt.gz"),header = T,row.names="Gene")
BM3<- CreateSeuratObject(counts = BM3, min.cells = 3, project = "BM3")
BM3_anno<- read.table(gzfile("GSM3587999_BM3.anno.txt.gz"),header = T,sep="\t")
BM3_anno$Cell<- gsub("-",".",BM3_anno$Cell)
a<-BM3@meta.data
a$Cell<-rownames(BM3@meta.data)
c<-merge(a,BM3_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
BM3@meta.data<-c
BM4<- read.table(gzfile("GSM3588000_BM4.dem.txt.gz"),header = T,row.names="Gene")
BM4<- CreateSeuratObject(counts = BM4, min.cells = 3, project = "BM4")
BM4_anno<- read.table(gzfile("GSM3588001_BM4.anno.txt.gz"),header = T,sep="\t")
BM4_anno$Cell<- gsub("-",".",BM4_anno$Cell)
a<-BM4@meta.data
a$Cell<-rownames(BM4@meta.data)
c<-merge(a,BM4_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
BM4@meta.data<-c
BM5_34p<- read.table(gzfile("GSM3588002_BM5_34p.dem.txt.gz"),header = T,row.names="Gene")
BM5_34p<- CreateSeuratObject(counts = BM5_34p, min.cells = 3, project = "BM5_34p")
BM5_34p_anno<- read.table(gzfile("GSM3588002_BM5_34p.anno.txt.gz"),header = T,sep="\t")
BM5_34p_anno$Cell<- gsub("-",".",BM5_34p_anno$Cell)
a<-BM5_34p@meta.data
a$Cell<-rownames(BM5_34p@meta.data)
c<-merge(a,BM5_34p_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
BM5_34p@meta.data<-c
BM5_34p38n<- read.table(gzfile("GSM3588003_BM5_34p38n.dem.txt.gz"),header = T,row.names="Gene")
BM5_34p38n<- CreateSeuratObject(counts = BM5_34p38n, min.cells = 3, project = "BM5_34p38n")
BM5_34p38n_anno<- read.table(gzfile("GSM3588003_BM5_34p38n.anno.txt.gz"),header = T,sep="\t")
BM5_34p38n_anno$Cell<- gsub("-",".",BM5_34p38n_anno$Cell)
a<-BM5_34p38n@meta.data
a$Cell<-rownames(BM5_34p38n@meta.data)
c<-merge(a,BM5_34p38n_anno,by="Cell")
rownames(c)=c[,1]
c=c[,-1]
BM5_34p38n@meta.data<-c
AML <- merge(AML1012_D0, y = c(AML210A_D0,AML314_D0,AML314_D31,AML328_D0,AML328_D113,AML328_D171,AML328_D29,AML329_D0,AML329_D20,AML329_D37,AML371_D0,AML371_D34,AML419A_D0,AML420B_D0,AML420B_D14,AML420B_D35,AML475_D0,AML475_D29,AML556_D0,AML556_D15,AML556_D31,AML707B_D0,AML707B_D113,AML707B_D18,AML707B_D41,AML707B_D97,AML722B_D0,AML722B_D49,AML870_D0,AML870_D14,AML916_D0,AML921A_D0,AML997_D0,AML997_D35,BM1,BM2,BM3,BM4,BM5_34p,BM5_34p38n), add.cell.ids = c("AML1012_D0","AML210A_D0","AML314_D0","AML314_D31","AML328_D0","AML328_D113","AML328_D171","AML328_D29","AML329_D0","AML329_D20","AML329_D37","AML371_D0","AML371_D34","AML419A_D0","AML420B_D0","AML420B_D14","AML420B_D35","AML475_D0","AML475_D29","AML556_D0","AML556_D15","AML556_D31","AML707B_D0","AML707B_D113","AML707B_D18","AML707B_D41","AML707B_D97","AML722B_D0","AML722B_D49","AML870_D0","AML870_D14","AML916_D0","AML921A_D0","AML997_D0","AML997_D35","BM1","BM2","BM3","BM4","BM5_34p","BM5_34p38n"), project = "AML")
library(SeuratDisk)
SaveLoom(AML, filename = "cell2019AML.loom", verbose = FALSE)
