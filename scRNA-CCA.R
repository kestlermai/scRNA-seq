if(!require(Seurat))BiocManager::install("Seurat")
if(!require(ggplot2))BiocManager::install("ggplot2")
if(!require(dplyr))BiocManager::install("dplyr")
if(!require(patchwork))BiocManager::install("patchwork")

getwd()
setwd("D:/scRNA")
setwd("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV")

samples <- list.files('raw/')
samples
scList <- lapply(samples,function(pro){ 
  folder <- file.path('raw/',pro) 
  print(pro)
  print(folder)
  print(list.files(folder))
  sc <- CreateSeuratObject(counts = Read10X(folder), project = pro, min.cells = 3, min.features = 200)
  return(sc)
})

names(scList) <- samples
names(scList)
scList
#merge数据
sce.all <- merge(scList[[1]], y = scList[-1], add.cell.ids = samples) 

#QC
dir.create("./1-QC-7.28")
setwd("./1-QC-8.8")
#计算线粒体基因比例
mito_genes <- rownames(sce.all)[grep("^MT-", rownames(sce.all))]
mito_genes #13个线粒体基因
sce.all <- PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
fivenum(sce.all@meta.data$percent_mito)
#计算核糖体基因比例
ribo_genes <- rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
ribo_genes
sce.all <- PercentageFeatureSet(sce.all, "^RP[SL]", col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)#最小值、下四分位数、中位数、上四分位数、最大值
#计算红血细胞基因比例
rownames(sce.all)[grep("^HB[^(p)]", rownames(sce.all))]
sce.all <- PercentageFeatureSet(sce.all, "^HB[^(p)]", col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)
#可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA")
p1 <- VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
ggsave(filename="Vlnplot1.pdf", p1, width = length(unique(sce.all$orig.ident))/3+5, height = 5, dpi = 600)
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2 <- VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims = T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
ggsave(filename="Vlnplot2.pdf", p2, width = length(unique(sce.all$orig.ident))/3+5, height = 5, dpi = 600)
#
feats <- c("nFeature_RNA", "nCount_RNA","percent_mito")
p1 <- VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
ggsave(filename="Vlnplot_raw.pdf", p1, width = length(unique(sce.all$orig.ident))/3+5, height = 5, dpi = 600)


p3.feature <- FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p3.mito <- FeatureScatter(sce.all, "nCount_RNA", "percent_mito", group.by = "orig.ident", pt.size = 0.5)
p3.rb <- FeatureScatter(sce.all, "nCount_RNA", "percent_ribo", group.by = "orig.ident", pt.size = 0.5)
p3.rb.mito <- FeatureScatter(sce.all, "percent_ribo", "percent_mito", group.by = "orig.ident", pt.size = 0.5)
p3 <- p3.feature + p3.mito #+ p3.rb + p3.rb.mito
#p31 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent_mito") 
#p32 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent_ribo") 
#p33 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent_hb")
#p34 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#p3 <- p31 + p32 + p33 + p34
ggsave(filename="Scatterplot.pdf", p3, width = 10, height = 6, dpi = 600)

#常见的标准是，保留线粒体基因比率<10%，200 < nFeature_RNA < 5000 or 6000
#如果质控更加严格的话，可以保留线粒体基因比率<5%, 200 < nFeature_RNA < 2500
#过滤指标:最少表达基因数的细胞&最少表达细胞数的基因
selected_c <- WhichCells(sce.all, expression = nFeature_RNA > 200 & nFeature_RNA < 5000)#600
selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA@counts > 0 ) > 3]
sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)
dim(sce.all)
#过滤前 30163 56866
dim(sce.all.filt) 
#过滤后 28778 56623
#可以看到，主要是过滤了基因，其次才是细胞

#展现top50高表达基因
C <- sce.all.filt@assays$RNA@counts
dim(C)
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]
pdf("TOP50_most_expressed_gene.pdf",width = 14,height = 12)
boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
        cex = 0.1, las = 1, 
        xlab = "% total count per cell", 
        col = (scales::hue_pal())(50)[50:1], 
        horizontal = TRUE)
dev.off()
rm(C)

#过滤线粒体/核糖体基因比例(根据上面的violin图)
selected_mito <- WhichCells(sce.all.filt, expression = percent_mito < 15)
#selected_ribo <- WhichCells(sce.all.filt, expression = percent_ribo > 3)
#selected_hb <- WhichCells(sce.all.filt, expression = percent_hb < 0.1)
length(selected_hb)
length(selected_ribo)
length(selected_mito)#53609

sce.all.filt <- subset(sce.all.filt, cells = selected_mito)#只过滤线粒体
dim(sce.all.filt)#28778 53609
sce.all.filt <- subset(sce.all.filt, cells = selected_ribo)
dim(sce.all.filt)#18628 21547
sce.all.filt <- subset(sce.all.filt, cells = selected_hb)
dim(sce.all.filt)#18628 21344
table(sce.all.filt$orig.ident)
#
#HC_1  HC_2  HC_3 HIV_1 HIV_2 HIV_3 HIV_4 HIV_5 HIV_6 
#6831  8377 12576  3787  4454  5391  4050  4149  3994

#可视化过滤后的情况
feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered <- VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
ggsave(filename="Vlnplot1_filtered.pdf", p1_filtered, width = length(unique(sce.all.filt$orig.ident))/3+5, height = 5, dpi = 600)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered <- VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
ggsave(filename="Vlnplot2_filtered.pdf", p2_filtered, width = length(unique(sce.all.filt$orig.ident))/3+5, height = 5, dpi = 600)

feats <- c("nFeature_RNA", "nCount_RNA","percent_mito")
p1_filtered <- VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
ggsave(filename="Vlnplot1_filtered.pdf", p1_filtered, width = length(unique(sce.all.filt$orig.ident))/3+5, height = 5, dpi = 600)

print(dim(sce.all))#30163 56866
print(dim(sce.all.filt))#28778 53609

#过滤特定基因
### Filter MALAT1 管家基因
sce.all.filt <- sce.all.filt[!grepl("MALAT1", rownames(sce.all.filt),ignore.case = T), ]
dim(sce.all.filt)#28777 53609
### Filter Mitocondrial 线粒体基因
sce.all.filt <- sce.all.filt[!grepl("^MT-", rownames(sce.all.filt),ignore.case = T), ]
dim(sce.all.filt)#28764 53609
### Filter 核糖体基因
sce.all.filt <- sce.all.filt[!grepl("^RP[SL]", rownames(sce.all.filt),ignore.case = T), ]
dim(sce.all.filt)#28658 53609

# 在分析单细胞数据时，同一类型的细胞往往来自于不同的细胞周期阶段，这可能对下游聚类分析，细胞类型注释产生混淆；
# 由于细胞周期也是通过cell cycle related protein 调控，即每个阶段有显著的marker基因；
# 通过分析细胞周期有关基因的表达情况，可以对细胞所处周期阶段进行注释；
# 在单细胞周期分析时，通常只考虑三个阶段：G1、S、G2M。
# 细胞周期评分
sce.all.filt <- NormalizeData(sce.all.filt, normalization.method = "LogNormalize", scale.factor = 10000)
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
sce.all.filt <- CellCycleScoring(object = sce.all.filt, 
                              s.features = s.genes, 
                              g2m.features = g2m.genes, 
                              set.ident = TRUE)
p4 <- VlnPlot(sce.all.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
           ncol = 2, pt.size = 0)
ggsave(filename="Vlnplot4_cycle.pdf", p4, width = 8, height = 5, dpi = 600)

sce.all.filt@meta.data %>% ggplot(aes(S.Score,G2M.Score))+ geom_point(aes(color = Phase))+
  theme_minimal()
ggsave(filename="cycle_details.pdf", width = 8, height = 5, dpi = 600)
# S.Score较高的为S期，G2M.Score较高的为G2M期，都比较低的为G1期

library(ggsci)
#寻找高变基因
sce.all.filt <- FindVariableFeatures(sce.all.filt, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sce.all.filt), 10)
top10
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sce.all.filt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
top10.Variable <- plot1 + plot2
ggsave(filename="top10_VariableFeatures.pdf", top10.Variable, width = 12, height = 8, dpi = 600)
#需要先归一化
sce.all.filt <- ScaleData(sce.all.filt, features = rownames(sce.all.filt))

sce.all.filt <- RunPCA(sce.all.filt, features = c(s.genes, g2m.genes))
p.ccs <- DimPlot(sce.all.filt, cols = pal_npg("nrc", alpha = 0.7)(3))
#明显聚在一起，没有明显的区分
ggsave(filename="cycle_pca.pdf", p.ccs, width = 8, height = 6, dpi = 600)
save(sce.all,sce.all.filt,file = 'sce.all.filt_qc_no_doublet.Rdata')
#一般双细胞是细胞量大于6000做
#检测doublets 
# 筛选高变基因
sce.all.filt <- FindVariableFeatures(sce.all.filt, selection.method = "vst", nfeatures = 2000, verbose = F)
#对数据进行缩放，缩放的参数vars.to.regress按照自己的目的决定，一般选择percent.mt，nFeature
sce.all.filt <- ScaleData(sce.all.filt, 
                         vars.to.regress = c("nFeature_RNA", "percent_mito"))
#降维
#确定降维的PC数了，具体选择多少比较合适，这个需要不断的尝试，没有标准，达到自己理想的效果即可。
#一般10都可以解释数据90%的信息
sce.all.filt <- RunPCA(sce.all.filt, features = VariableFeatures(object = sce.all.filt))
VizDimLoadings(sce.all.filt , dims = 1:2, reduction = "pca")
DimHeatmap(sce.all.filt, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(object = sce.all.filt, ndims = 20)#使用JackStraw和Elbow图方法确定了下游聚类数据集的最佳维数。
#文章里面用的是PC 10
sce.all.filt <- RunPCA(sce.all.filt, npcs = 10)#提取20个主成分
sce.all.filt <- RunTSNE(sce.all.filt, npcs = 10)
sce.all.filt <- RunUMAP(sce.all.filt, dims = 1:10)#保留前 10 个主成分
sce.all.filt
save(sce.all.filt,file = 'sce.all.filt_qc_no_doublet.Rdata')
#大多数此类软件包都需要对数据集中预期双细胞数的数量/比例进行假设。我们这里使用的数据是经过二次抽样的
#但是原始数据集中每个样本包含了大约5000个细胞，因此我们可以假设它们加载了大约9000个细胞，并且其双细胞率约为4％。
# 5000细胞对应的doublets rate是3.9%
#8个样共31732个细胞，平均4000个细胞
table(sce.all.filt$orig.ident)
sce.all.list <- SplitObject(sce.all.filt, split.by = "orig.ident")
sce.all.list[[3]]
### define the expected number of doublet cellscells.
nExp <- round(ncol(sce.all.filt) * 0.04)  ### expect 4% doublets也可以直接按4%去除
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
#这个DoubletFinder包的输入是经过预处理（包括归一化、降维，但不一定要聚类）的 Seurat 对象
sce.all.filt <- doubletFinder_v3(sce.all.filt, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:20)
###### 找Doublets  
### DF的名字是不固定的，因此从sce.all.filt@meta.data列名中提取比较保险
DF.name <- colnames(sce.all.filt@meta.data)[grepl("DF.classification", colnames(sce.all.filt@meta.data))]
p5.dimplot <- cowplot::plot_grid(ncol = 2, DimPlot(sce.all.filt, group.by = "orig.ident") + NoAxes(), 
                              DimPlot(sce.all.filt, group.by = DF.name) + NoAxes())
ggsave(filename="doublet_dimplot.pdf", p5.dimplot, width = 10, height = 6, dpi = 600)
#(图左)红色部分上的黑点就是doublets
p5.vlnplot <- VlnPlot(sce.all.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
ggsave(filename="doublet_vlnplot.pdf", p5.vlnplot, width = 8, height = 6, dpi = 600)

###### 过滤doublet
sce.all.filt <- sce.all.filt[, sce.all.filt@meta.data[, DF.name] == "Singlet"]
table(sce.all.filt$orig.ident)
table(sce.all$orig.ident)
save(sce.all.filt,file = 'sce.all.filt_qc.Rdata')

#sce.all.filt <- FindNeighbors(sce.all.filt, dims = 1:30)
#sce.all.filt <- FindClusters(sce.all.filt, resolution = 0.5)#一般0.5都可以分的很细
#table(sce.all.filt@meta.data$RNA_snn_res.0.5)


#CCA
#CCA去批次,对内存有要求，没有服务器用harmony跑吧
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
rm(list = ls())
setwd("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/1-QC-7.28")
load('sce.all.filt_qc.Rdata')

sce.all <- sce.all.filt
sce.all#19314 features across 30780 samples
table(sce.all@meta.data$orig.ident)
#拆分为 个seurat子对象
sce.all.list <- SplitObject(sce.all, split.by = "orig.ident")
##对每项数据集进行归一化和寻找高变基因
for (i in 1:length(sce.all.list)) {
  print(i)
  sce.all.list[[i]] <- NormalizeData(sce.all.list[[i]], verbose = FALSE)
  sce.all.list[[i]] <- FindVariableFeatures(sce.all.list[[i]], selection.method = "vst", 
                                            nfeatures = 2000, verbose = FALSE)
}
#寻找锚点细胞
alldata.anchors <- FindIntegrationAnchors(object.list = sce.all.list, dims = 1:30, 
                                          reduction = "cca")
##整合高变基因找出两两数据集间的锚点细胞
sce.all.int <- IntegrateData(anchorset = alldata.anchors, dims = 1:30, new.assay.name = "CCA")
#运行IntegrateData之后，Seurat对象将包含一个具有整合（或“批次校正”）表达矩阵的新Assay。
#请注意，原始值（未校正的值）仍存储在“RNA”分析的对象中，因此可以来回切换。然后可以使用这个新的整合矩阵进行下游分析和可视化
names(sce.all.int@assays)
#[1] "RNA" "CCA"
sce.all.int@active.assay#查看所用的assay
#[1] "CCA"

##执行将为算法
sce.all <- RunPCA(sce.all, npcs = 10)
sce.all <- RunTSNE(sce.all, dims = 1:10)
sce.all <- RunUMAP(sce.all, dims = 1:10)

sce.all.int <- ScaleData(sce.all.int)
sce.all.int <- RunPCA(sce.all.int, npcs = 10)
sce.all.int <- RunTSNE(sce.all.int, dims = 1:10)
sce.all.int <- RunUMAP(sce.all.int, dims = 1:10)
names(sce.all.int@reductions)
names(sce.all@reductions) 
colnames(sce.all@meta.data)

#比较整合前后的降维结果
p1.compare <- plot_grid(ncol = 3,
                     DimPlot(sce.all, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA raw_data"),
                     DimPlot(sce.all, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE raw_data"),
                     DimPlot(sce.all, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP raw_data"),
                     DimPlot(sce.all.int, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA integrated"),
                     DimPlot(sce.all.int, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE integrated"),
                     DimPlot(sce.all.int, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP integrated")
)
ggsave(filename="Before&After_int_cca1.pdf", p1.compare, width = 15, height = 8, dpi = 600) 

save(sce.all,sce.all.int, file = 'sce.all.int_CCA.Rdata')

#还是有些区别的，后面还是用已去除DF的数据
#sce.all.int <- sce.all.int.DF

load("sce.all.int_CCA.Rdata")
##分cluster
sce.all.int #After
sce.all #Before
sce.all <- sce.all.int
sce.all@active.assay
#使用Seurat包中的FindNeighbors函数计算构建SNN图。
sce.all <- FindNeighbors(sce.all, dims = 1:20, k.param = 60, prune.SNN = 1/15)
#设置不同的分辨率，观察分群效果(选择哪一个？)
for (res in c(0.01,0.05,0.1,0.2,0.3,0.5,0.8,1,1.2,1.5,2)) {
  sce.all <- FindClusters(sce.all, graph.name = "CCA_snn", resolution = res, algorithm = 1)}
apply(sce.all@meta.data[,grep("CCA_snn_res",colnames(sce.all@meta.data))],2,table)

p2_tree <- clustree(sce.all@meta.data, prefix = "CCA_snn_res.")
ggsave(p2_tree, filename="Tree_diff_resolution.pdf", width = 8, height = 10, dpi = 600)

#低分辨率的分群情况。（0.01，0.1，0.2）
p1_dim <- plot_grid(ncol = 3, DimPlot(sce.all, reduction = "umap", group.by = "CCA_snn_res.0.01") + 
                   ggtitle("louvain_0.01"), 
                   DimPlot(sce.all, reduction = "umap", group.by = "CCA_snn_res.0.1") + 
                   ggtitle("louvain_0.1"), 
                   DimPlot(sce.all, reduction = "umap", group.by = "CCA_snn_res.0.2") + 
                   ggtitle("louvain_0.2"))
ggsave(p1_dim, filename="Dimplot_diff_resolution_low.pdf", width = 14, height = 6, dpi = 600)


#高分辨率的分群情况。（0.3，0.5, 0.8，1）
p2_dim <- plot_grid(ncol = 2, DimPlot(sce.all, reduction = "umap", group.by = "CCA_snn_res.0.3") + 
                   ggtitle("louvain_0.3"), 
                   DimPlot(sce.all, reduction = "umap", group.by = "CCA_snn_res.0.5") + 
                     ggtitle("louvain_0.5"), 
                 DimPlot(sce.all, reduction = "umap", group.by = "CCA_snn_res.0.8") + 
                   ggtitle("louvain_0.8"), 
                 DimPlot(sce.all, reduction = "umap", group.by = "CCA_snn_res.1") + 
                   ggtitle("louvain_1"))
ggsave(p2_dim, filename="Dimplot_diff_resolution_high.pdf", width = 12, height = 8, dpi = 600)

#接下来分析，按照分辨率为0.8进行 
sel.clust <- "CCA_snn_res.0.8"
sce.all <- SetIdent(sce.all, value = sel.clust)
table(sce.all@active.ident) 
saveRDS(sce.all, "sce.all_int.rds")

##细胞类型注释(根据marker gene)
setwd("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8")
rm(list=ls())
sce.all <- readRDS("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/1-QC-8.8/sce.all_int.rds")
table(sce.all@meta.data$seurat_clusters)
table(sce.all@meta.data$CCA_snn_res.0.8)
#展示按照不同标准进行分群的结果，并没有什么区别
p1 <- DimPlot(sce.all, reduction = "umap", group.by = "seurat_clusters", label = T) 
ggsave("umap_by_seurat_clusters.pdf", p1, width = 8, height = 7, dpi = 600)
p2 <- DimPlot(sce.all, reduction = "umap", group.by = "CCA_snn_res.0.8", label = T) 
ggsave("umap_by_CCA_snn_res.0.8.pdf", p2, width = 8, height = 7, dpi = 600)
p3 <- DimPlot(sce.all, reduction = "tsne", group.by = "CCA_snn_res.0.8", label = T) 
p3ggsave("tsne_by_CCA_snn_res.0.8.pdf", p3, width = 8, height = 7, dpi = 600)
p4 <- DimPlot(sce.all, reduction = "umap", group.by = "CCA_snn_res.1", label = T) 
ggsave("umap_by_CCA_snn_res.1.pdf", p4, width = 8, height = 7, dpi = 600)
p5 <- DimPlot(sce.all, reduction = "tsne", group.by = "CCA_snn_res.1", label = T) 
ggsave("tsne_by_CCA_snn_res.1.pdf", p5, width = 8, height = 7, dpi = 600)


DefaultAssay(sce.all) <- "RNA"
markerGenes <- c("CD3D", "CD8B",#"TRAC",# 定位T细胞
                 #"TRAC",# 定位T细胞
                 "IL7R", # CD4 T cells
                 #"GZMA", # NK T /效应T
                 "GNLY",
                 "NKG7","GZMB", # NK cells
                 #"CD8B","GZMK",# CD8 T cells
                 "FCGR3A", # CD16(FCGR3A)+ Mono / NK / 效应 CD8+ T
                 "CD14", # CD14+ monocyte
                 "MS4A1", #B细胞
                 #"FCER1A","LILRA4","TPM2", #DC
                 "PPBP",#"GP1BB"# platelets
                 "CCR7","CD27","LYZ","S100A9",
                 "MZB1"#Plasma
)
#Nucleic Acids Res. 2018 Apr ，dropClust: efficient clustering of ultra-large scRNA-seq data.
markerGenes <- c("IL7R","CCR7","CD8A","CD8B",# Naive CD4 T cells
                 "IL7R","CD27","CCR7",#Memory CD4+
                 "ZNF683","CD8A","CD8B",#NK T Cell
                 "CD79A","CD37","MS4A1",#B
                 "GZMK","CD8A","CD8B",#CD8+ T 
                 "CD160","NKG7","GNLY","CD247","CCL3","GZMB",# NK
                 "CD68","CD16","CD14","S100A12","FCGR3A",#CD16+ and CD14+ monocytes
                 "CCR10","CD25","CD52","CMTM7","FOXP3",#Regulatory T
                 "CST3","CD1C","FCERA1",#Monocyte derived dendritic cells树突细胞
                 "PF4","PPBP","PLA2G12A",#Megakaryocyte progenitors巨核细胞
                 "ID2",#Progenitor-NK cells祖细胞NK细胞
                 "GZMB","CD123",#Plasmacytoid dendritic cells浆细胞样树突状细胞
                 "LYZ", "MARCO", "FCGR1A", "C1QB"#Macrophages巨噬细胞
)
p1 <- VlnPlot(sce.all, features = markerGenes, pt.size = 0, stack = T)+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 12))+
  NoLegend()
p1
ggsave("VlnPlot_markerGenes-stack.pdf", p1, width = 18, height = 12, dpi = 600)
p2 <- DotPlot(sce.all, features = markerGenes, dot.scale = 8) + RotatedAxis() + coord_flip()
ggsave("DotPlot_markerGenes.pdf", p2, width = 12, height = 8, dpi = 600)
p3 <- FeaturePlot(sce.all, features = markerGenes, label.size = 4, repel = T, label = T)+
  theme(plot.title = element_text(size = 8))
ggsave("FeaturePlot_markerGenes.pdf", p3, width = 15, height = 10, dpi = 600)
#识别每个cluster的marker
markers <- FindAllMarkers(sce.all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
DoHeatmap(sce.all, features = top10$gene) + NoLegend()

#在不同条件下识别保守marker
#由于使用的数据集中有不同条件的样本，因此最好的选择是找到保守标记
#在开始标记识别之前，我们将明确设置使用原始计数而不是整合数据
DefaultAssay(sce.all) <- "RNA"
#提取所有样本名称
all_samples <- colnames(sce.all)
## 将包含"GSM4775588"和"pbmc3k"字符串的细胞样本分配为"control"组
control_samples <- grep("HC_", all_samples, value = TRUE)
sce.all@meta.data$group <- ifelse(all_samples %in% control_samples, "control", "HIV")
table(sce.all@meta.data$group)

if(!require(multtest))BiocManager::install('multtest')
if(!require(metap))BiocManager::install('metap')
conserved_markers <- FindConservedMarkers(sce.all, ident.1 = "Neutrophil", #此参数表示一次仅评估一个cluster
                                grouping.var = "group", #metadata中的变量（列标题），标记细胞对应的样品来源；
                                only.pos = TRUE, min.diff.pct = 0.25, 
                                min.pct = 0.25, logfc.threshold = 0.25)

celltype <- data.frame(ClusterID = 0:14, celltype = 'unkown')
celltype[celltype$ClusterID %in% c(0,7,13),2]='NK cells'
celltype[celltype$ClusterID %in% c(1,2,4,5,6,12),2]='T cells'
celltype[celltype$ClusterID %in% c(11),2]='Megakaryocytes'
#celltype[celltype$ClusterID %in% c(10),2]='Neutrophils'
celltype[celltype$ClusterID %in% c(14),2]='pDCs'
celltype[celltype$ClusterID %in% c(3,10),2]='Monocytes'
#celltype[celltype$ClusterID %in% c(9),2]='CD16+ Mono'
celltype[celltype$ClusterID %in% c(8,9),2]='B cells'
#celltype[celltype$ClusterID %in% c(12),2]='Mk'

head(celltype)
table(celltype$celltype)
sce.all@meta.data$celltype = "NA"
#先新增列celltype，值均为NA，然后利用下一行代码循环填充
for(i in 1:nrow(celltype)){
  sce.all@meta.data[which(sce.all@meta.data$CCA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.all@meta.data$celltype)
p1_tsne <- DimPlot(sce.all, reduction = "tsne", group.by = "celltype",label = TRUE, pt.size = 1)
p2_umap <- DimPlot(sce.all, reduction = "umap", group.by = "celltype",label = TRUE, pt.size = 1)
p3_group <- DimPlot(sce.all, reduction = "tsne", group.by = "group",label = F, pt.size = 1)
p4_split.tsne <- DimPlot(sce.all,reduction = "tsne", label = T, split.by ='group',pt.size = 1)
p4_split.umap <- DimPlot(sce.all,reduction = "umap", label = T, split.by ='group',pt.size = 1)
mycol <- c("#E15759","#59A14F","#76B7B2","#F28E2B","#EDC948","#4E79A7")
p5_split.umap.group <- DimPlot(sce.all,reduction = "umap", label = T, split.by ='group', cols = mycol,
                               group.by = "celltype",pt.size = 1)
ggsave(filename="tsne_marker.pdf", p1_tsne, width = 10, height = 8, dpi = 600) 
ggsave(filename="umap_marker.pdf", p2_umap, width = 10, height = 8, dpi = 600) 
ggsave(filename="group_marker.pdf", p3_group, width = 10, height = 8, dpi = 600) 
ggsave(filename="split_marker.tsne.pdf", p4_split.tsne, width = 16, height = 8, dpi = 600) 
ggsave(filename="split_marker.umap.pdf", p4_split.umap, width = 16, height = 8, dpi = 600) 
ggsave(filename="split_marker.umap_group.pdf", p5_split.umap.group, width = 18, height = 10, dpi = 600) 
saveRDS(sce.all,"sce.all_marker.rds")
#
setwd("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8")
sce.all <- readRDS("sce.all_marker.rds")

#用qs包进行rdata储存更加迅速，Quick Serialization of R Objects
# install.packages("qs")
# library(qs)
# 保存
# qsave(sce.all, "sce.all_marker.qs")
# 读取
# sce.all <- qread("sce.all_marker.qs")

#提取特定细胞,譬如B细胞
B_cells <- subset(sce.all, cells = rownames(sce.all@meta.data[sce.all@meta.data$celltype=="B cells",]))
dim(B_cells)
dim(sce.all)
#用seurat的ident进行celltype进行赋值，多一行代码，不推荐，不如上面
sce.all$celltype <- Idents(sce.all)
#提取特定细胞,譬如B细胞
B_cells <- subset(sce.all, idents = 'B cells')

#harmony去批次
if(!require(harmony))devtools::install_github("immunogenomics/harmony")
rm(list = ls())
setwd("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/1-QC-7.28")
load('sce.all.filt.Rdata')
sce.all <- sce.all.filt
table(sce.all.filt$orig.ident)

sce.all.filt <- NormalizeData(sce.all.filt) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)

#添加一个肘图参数，可以更好的确定pca阈值
sce.all.int <- RunHarmony(sce.all.filt, c( "orig.ident" ), plot_convergence = T)
ggsave(filename = "Harmony_plot_convergence.pdf", width = 10, height = 8, dpi = 600)

harmony_embeddings <- Embeddings(sce.all.int, 'harmony')

sce.all.int <- FindNeighbors(sce.all.int, reduction = "harmony", dims = 1:10) %>% FindClusters(resolution = 0.8)

#注意一下，这里执行降维算法的时候，是基于harmony，不指定的话，默认是PCA
sce.all.int <- RunTSNE(sce.all.int, reduction = "harmony", dims = 1:20)
sce.all.int <- RunUMAP(sce.all.int, reduction = "harmony", dims = 1:20)
names(sce.all.int@reductions)
names(sce.all@reductions) 
p.umap_after_harmony <- DimPlot(sce.all.int, reduction = "umap", group.by = 'orig.ident')
p.tsne_after_harmony <- DimPlot(sce.all.int, reduction = "tsne", group.by = 'orig.ident')
p.after_harmony <- p.umap_after_harmony + p.tsne_after_harmony
ggsave('by_orig.ident_after_harmony.pdf', p.after_harmony, width = 18, height = 8, dpi = 600)
#比对整合结果
sce.all <- sce.all.filt
p1.compare <- plot_grid(ncol = 3,
                     DimPlot(sce.all, reduction = "pca", group.by = "orig.ident")+ NoAxes()+ggtitle("PCA raw_data"),
                     DimPlot(sce.all, reduction = "tsne", group.by = "orig.ident")+ NoAxes()+ggtitle("tSNE raw_data"),
                     DimPlot(sce.all, reduction = "umap", group.by = "orig.ident")+ NoAxes()+ggtitle("UMAP raw_data"),
                     DimPlot(sce.all.int, reduction = "pca", group.by = "orig.ident")+ NoAxes()+ggtitle("PCA integrated"),
                     DimPlot(sce.all.int, reduction = "tsne", group.by = "orig.ident")+ NoAxes()+ggtitle("tSNE integrated"),
                     DimPlot(sce.all.int, reduction = "umap", group.by = "orig.ident")+ NoAxes()+ggtitle("UMAP integrated")
)
ggsave(filename="Before&After_int_harmony.pdf", p1.compare, width = 15, height = 8, dpi = 600) 

##分cluster
sce.all <- sce.all.int
sce.all <- FindNeighbors(sce.all, reduction = "harmony", dims = 1:20)
colnames(sce.all@meta.data)
sce.all@graphs$RNA_snn
#设置不同的分辨率，观察分群效果(选择哪一个？)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1)){
  sce.all <- FindClusters(sce.all, resolution = res, algorithm = 1)
}
# 提取以 "RNA_snn_res" 开头的全部列名
apply(sce.all@meta.data[,grep("RNA_snn_res", colnames(sce.all@meta.data))],2,table)

p2_tree <- clustree(sce.all@meta.data, prefix = "RNA_snn_res.")
ggsave(p2_tree, filename="Tree_diff_resolution.pdf", width = 7, height = 10, dpi = 600)
#可以看到从 res=0.3开始，亚群间的细胞交换逐渐混乱，所以后续的分析选 0.3。
table(sce.all@meta.data$RNA_snn_res.0.5)

#主要是通过ROC分析来对marker基因的表达特异性进行评估，也就是将marker基因作为分类器来对细胞进行分类
#使用FindAllMarkers()函数计算roc值
#marker gene AUC 值
sce.all %>% 
  SetIdent(value = 'RNA_snn_res.0.5') %>% 
  FindAllMarkers(test.use = 'roc') %>% 
  filter(myAUC > 0.6) %>% 
  count(cluster, name = 'number')

sce.all <- SetIdent(sce.all, value = sce.all$RNA_snn_res.0.5)
#注释细胞
DotPlot(sce.all, features = c("MS4A1", "GNLY", "CD3E", 
                              "CD14", "FCER1A", "FCGR3A", 
                              'C1QA',  'C1QB',  # mac
                              'S100A9', 'S100A8', 'MMP19',# monocyte
                              'LAMP3', 'IDO1','IDO2',## DC3 
                              'CD1E','CD1C', # DC2
                              "LYZ", "PPBP", "CD8A",
                              'CD203c','CD63', 'CD123',"CD68"), assay = "RNA") + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) 
ggsave("anno_marker_gene_check.pdf", width = 10, height = 8, dpi = 600)

DimPlot(sce.all, reduction = "umap", pt.size = 1.0, label = T)#+ NoLegend()
ggsave(filename = "UMAP_res_0.3.pdf", width = 8, height = 7, dpi = 600)

DimPlot(sce.all, reduction = "tsne", pt.size = 1.0, label = T)#+ NoLegend()
ggsave(filename = "TSNE_res_0.3.pdf", width = 8, height = 7, dpi = 600)

markers <- FindAllMarkers(sce.all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
DoHeatmap(sce.all, features = top10$gene, assay = "CCA") + NoLegend()
ggsave(filename = "top10_heatmap_res_0.8.pdf", width = 18, height = 10)

saveRDS(sce.all, "sce.all.rds") 

sce.all <- readRDS("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/1-QC/sce.all.rds")
