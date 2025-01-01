if(!require(Seurat))BiocManager::install("Seurat")
if(!require(ggplot2))BiocManager::install("ggplot2")
if(!require(dplyr))BiocManager::install("dplyr")
if(!require(patchwork))BiocManager::install("patchwork")
if(!require(cowplot))BiocManager::install("cowplot")

setwd("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/1-QC-8.8")
load('sce.all.filt_qc_no_doublet.Rdata')
sce.all.list <- SplitObject(sce.all.filt, split.by = "orig.ident")
sce.all.list[[2]]
#一般双细胞是细胞量大于6000做
#检测doublets 
# 筛选高变基因
sce.all.list[[2]] <- FindVariableFeatures(sce.all.list[[2]], verbose = F)
#对数据进行缩放，缩放的参数vars.to.regress按照自己的目的决定，一般选择percent.mt，nFeature
sce.all.list[[2]] <- ScaleData(sce.all.list[[2]], 
                          vars.to.regress = c("nFeature_RNA", "percent_mito"))

#PC 20
sce.all.list[[2]] <- RunPCA(sce.all.list[[2]], npcs = 20)#提取20个主成分
sce.all.list[[2]] <- RunTSNE(sce.all.list[[2]], npcs = 20)
sce.all.list[[2]] <- RunUMAP(sce.all.list[[2]], dims = 1:20)#保留前 10 个主成分
sce.all.list[[2]]
### define the expected number of doublet cellscells.
# 按每增加1000个细胞，双细胞比率增加千分之8来计算（推荐）
nExp <- round(ncol(sce.all.list[[2]]) * 0.067)  ### expect 4% doublets也可以直接按4%去除
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
#这个DoubletFinder包的输入是经过预处理（包括归一化、降维，但不一定要聚类）的 Seurat 对象
sce.all.list[[2]] <- doubletFinder_v3(sce.all.list[[2]], pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:20)
###### 找Doublets  
### DF的名字是不固定的，因此从sce.all.filt@meta.data列名中提取比较保险
DF.name <- colnames(sce.all.list[[2]]@meta.data)[grepl("DF.classification", colnames(sce.all.list[[2]]@meta.data))]
p5.dimplot <- cowplot::plot_grid(ncol = 2, DimPlot(sce.all.list[[2]], group.by = "orig.ident") + NoAxes(), 
                                 DimPlot(sce.all.list[[2]], group.by = DF.name) + NoAxes())
ggsave(filename="doublet_dimplot2.pdf", p5.dimplot, width = 10, height = 6, dpi = 600)
#(图左)红色部分上的黑点就是doublets
p5.vlnplot <- VlnPlot(sce.all.list[[2]], features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
ggsave(filename="doublet_vlnplot1.pdf", p5.vlnplot, width = 8, height = 6, dpi = 600)
###### 过滤doublet
sce.all.list[[2]] <- sce.all.list[[2]][, sce.all.list[[2]]@meta.data[, DF.name] == "Singlet"]
table(sce.all.list[[2]]$orig.ident)
table(sce.all$orig.ident)
#i=1,nExp <- round(ncol(sce.all.list[[i]]) * 0.054);
#i=2,nExp <- round(ncol(sce.all.list[[i]]) * 0.067);
#i=3,nExp <- round(ncol(sce.all.list[[i]]) * 0.1);
#i=4,nExp <- round(ncol(sce.all.list[[i]]) * 0.03);
#i=5,nExp <- round(ncol(sce.all.list[[i]]) * 0.035);
#i=6,nExp <- round(ncol(sce.all.list[[i]]) * 0.043);
#i=7,nExp <- round(ncol(sce.all.list[[i]]) * 0.032);
#i=8,nExp <- round(ncol(sce.all.list[[i]]) * 0.033);
#i=9,nExp <- round(ncol(sce.all.list[[i]]) * 0.032);

for (i in 9) {
  sce.all.list[[i]] <- FindVariableFeatures(sce.all.list[[i]], verbose = F)
  sce.all.list[[i]] <- ScaleData(sce.all.list[[i]], 
                                 vars.to.regress = c("nFeature_RNA", "percent_mito"))
  sce.all.list[[i]] <- RunPCA(sce.all.list[[i]], npcs = 20)
  sce.all.list[[i]] <- RunTSNE(sce.all.list[[i]], npcs = 20)
  sce.all.list[[i]] <- RunUMAP(sce.all.list[[i]], dims = 1:20)
  nExp <- round(ncol(sce.all.list[[i]]) * 0.032)
  sce.all.list[[i]] <- doubletFinder_v3(sce.all.list[[i]], pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:20)
  DF.name <- colnames(sce.all.list[[i]]@meta.data)[grepl("DF.classification", colnames(sce.all.list[[i]]@meta.data))]
  p5.dimplot <- cowplot::plot_grid(ncol = 2, DimPlot(sce.all.list[[i]], group.by = "orig.ident") + NoAxes(), 
                                   DimPlot(sce.all.list[[i]], group.by = DF.name) + NoAxes())
  ggsave(filename = paste0("doublet_dimplot_", i, ".pdf"), p5.dimplot, width = 10, height = 6, dpi = 600)
  p5.vlnplot <- VlnPlot(sce.all.list[[i]], features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
  ggsave(filename = paste0("doublet_vlnplot_", i, ".pdf"), p5.vlnplot, width = 8, height = 6, dpi = 600)
  sce.all.list[[i]] <- sce.all.list[[i]][, sce.all.list[[i]]@meta.data[, DF.name] == "Singlet"]
}

table(sce.all.list[[9]]$orig.ident)
table(sce.all$orig.ident)


sce.all <- sce.all.filt
sce.all#28658 features across 53609 samples
table(sce.all@meta.data$orig.ident)

#CCA合并数据
#寻找锚点细胞
alldata.anchors <- FindIntegrationAnchors(object.list = sce.all.list, dims = 1:20, 
                                          reduction = "cca")
##整合高变基因找出两两数据集间的锚点细胞
sce.all.int <- IntegrateData(anchorset = alldata.anchors, dims = 1:20, new.assay.name = "CCA")
names(sce.all.int@assays)
sce.all.int@active.assay

sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize", scale.factor = 10000)
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst", nfeatures = 2000, verbose = F)
sce.all <- ScaleData(sce.all)

sce.all <- RunPCA(sce.all, npcs = 20)
sce.all <- RunTSNE(sce.all, dims = 1:20)
sce.all <- RunUMAP(sce.all, dims = 1:20)

sce.all.int <- ScaleData(sce.all.int)
sce.all.int <- RunPCA(sce.all.int, npcs = 20)
sce.all.int <- RunTSNE(sce.all.int, dims = 1:20)
sce.all.int <- RunUMAP(sce.all.int, dims = 1:20)
#比较整合前后的降维结果
p1.compare <- plot_grid(ncol = 3,
                        DimPlot(sce.all, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA raw_data"),
                        DimPlot(sce.all, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE raw_data"),
                        DimPlot(sce.all, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP raw_data"),
                        DimPlot(sce.all.int, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA integrated"),
                        DimPlot(sce.all.int, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE integrated"),
                        DimPlot(sce.all.int, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP integrated")
)

ggsave(filename="Before&After_int_cca.pdf", p1.compare, width = 15, height = 8, dpi = 600)
save(sce.all,sce.all.int, file = 'sce.all.int_CCA.Rdata')
load('sce.all.int_CCA.Rdata')
p1.compare <- plot_grid(ncol = 2,
                        DimPlot(sce.all, reduction = "umap", group.by = "orig.ident")+ggtitle("raw_data"),
                        DimPlot(sce.all.int, reduction = "umap", group.by = "orig.ident")+ggtitle("CCA_integrated")
)
ggsave(filename="Before&After_int_cca_umap.pdf", p1.compare, width = 13, height = 6, dpi = 600)
