if(!require(Seurat))BiocManager::install("Seurat")# V5版本
if(!require(ggplot2))BiocManager::install("ggplot2")
if(!require(dplyr))BiocManager::install("dplyr")
if(!require(patchwork))BiocManager::install("patchwork")
if(!require(qs))BiocManager::install("qs")

# 来自GSE242889数据
setwd("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HCC/raw/raw_data")
# 由于解压的文件不是标准命名中的.gz格式
# 同时，这个应该是cellranger v2版本处理的10x，所以需要genes.tsv改为features.tsv（seurat v4好像能直接读取genes.tsv）

# 批量修改一下
folder_list <- list.dirs(path = "./", full.names = TRUE, recursive = FALSE)
# 遍历文件夹
for (folder in folder_list) {
  file_list <- list.files(path = folder, full.names = TRUE, recursive = FALSE)
  for (file in file_list) {
    # 检查文件扩展名
    extension <- tools::file_ext(file)
    file_name <- tools::file_path_sans_ext(basename(file))
    # 检查是否是 genes.tsv 文件
    if (file_name == "genes" && extension == "tsv") {
      new_file <- file.path(dirname(file), "features.tsv")
      file.rename(file, new_file)  # 将 genes.tsv 重命名为 features.tsv
      file <- new_file  # 更新 `file` 路径为重命名后的路径
      cat("文件重命名:", new_file, "\n")
    }
    # 压缩文件
    if (extension %in% c("tsv", "mtx")) {
      gz_file <- paste0(file, ".gz")
      cmd <- paste("gzip", file)
      system(cmd)# 运行系统命令压缩文件
      file.rename(paste0(file, ".gz"), gz_file)# 重命名压缩后的文件
      cat("压缩文件:", file, "\n")
    }
  }
}

# 读入数据，说实话我还挺喜欢scanpy的读入方式的，seurat感觉每次合并都有毒 
# 这里要排除第一个元素，因为第一个是根目录
samples <- list.dirs()[-1]
samples
# 给每个文件取名字
names(samples) <- gsub("^\\.\\/", "", list.files(recursive = F))
samples
# 批量读取
scList <- lapply(samples, function(pro){
  pro <- gsub("^\\.\\/", "", pro)
  folder <- file.path("./",pro)
  print(pro)
  print(folder)
  print(list.files(folder))
  sc <- CreateSeuratObject(counts = Read10X(folder), project = pro, min.cells = 5, min.features = 300)
  return(sc)
})
scList
# merge数据
sce.all <- merge(scList[[1]], y = scList[-1], add.cell.ids = samples)
dim(sce.all)
# 27371 59867
names(sce.all@assays$RNA@layers)
# # V5版本相较于V4在assays对象下面多出了layers的结构
# # 这里用merge函数合并后发现layers结构并没有合并，后续提取counts时数据不完整
# # 用JoinLayers函数对layers进行合并
# sce.all <- JoinLayers(sce.all)
# sce.all
# # 应该用harmony处理完后再进行合并，不能直接合并

# 熟悉的QC---补充一下UMI，Unique Molecular Identifier 每个细胞的唯一分子标识
# 在构建文库的时候，用于标记每个原始 RNA 分子
# 由于每个 RNA 分子都被赋予一个独特的随机标识符，UMI 可以帮助区分原始分子和扩增过程中产生的 PCR 重复
# 在传统 RNA 测序中，会通过 PCR 扩增操作将少量 cDNA 库扩增到足以满足测序需要的浓度。
# 但是，由于 PCR 过程中扩增效率的不均一性，会导致部分转录本被过度扩增，从而影响表达量的实际准确性。
# 引入 UMI 后，即使在 PCR 扩增过程中某个特定转录本扩增了很多倍，所有扩增的拷贝都带有相同的 UMI。
# 因此，测序数据中可以通过识别同一 UMI 并将其计为一个分子，从而恢复原始 RNA 分子的真实数量。
# 通常，每个细胞的UMI计数应该在500以上，低于这个阈值的细胞可能质量较差，可能需要更深的测序。
## 提一嘴：barcode是标记细胞的，UMI是标记RNA分子
dir.create("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HCC/1-QC")
setwd("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HCC/1-QC")
# 计算线粒体基因比例
mito_genes <- rownames(sce.all)[grep("^MT-", rownames(sce.all))]
mito_genes #13个线粒体基因
sce.all <- PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
fivenum(sce.all@meta.data$percent_mito)
# 0.06546645  7.76969405 11.36820926 20.08885978 94.31183342
# 计算核糖体基因比例
ribo_genes <- rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
ribo_genes
sce.all <- PercentageFeatureSet(sce.all, "^RP[SL]", col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)#最小值、下四分位数、中位数、上四分位数、最大值
# 0.188781  5.780844  9.314725 14.462306 57.947556
# 计算红血细胞基因比例
rownames(sce.all)[grep("^HB[^(p)]", rownames(sce.all))]
sce.all <- PercentageFeatureSet(sce.all, "^HB[^(p)]", col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)
# 0.000000000  0.000000000  0.004597701  0.023835772 93.257476743

# nFeature_RNA: 这是每个细胞中检测到的不同基因的数量
# 该指标反映了细胞的基因表达多样性。一个细胞中表达的基因数量越多，通常意味着这个细胞是健康的并且质量较高。
# 然而，太多的检测到的基因（远高于正常水平）可能表明细胞存在双重细胞（即两个细胞一起被测序），需要排除这种异常情况。

# nCount_RNA:这是每个细胞中测序到的所有 RNA 分子的总数=====UMI
# 该指标反映了测序深度和细胞的转录活性。总 RNA 分子数量越多，通常表示细胞活性越高。
# 然而，异常高的 RNA 分子总数也可能是双重细胞的一个标志。

# percent.mt
# 线粒体基因的高比例通常被认为是细胞压力或死亡的标志。细胞在死亡过程中，线粒体 RNA 通常会增加。
# 因此，一个细胞的 percent.mt 越高，通常表明该细胞的质量越低，可能需要被过滤掉。

# 可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA")
p1 <- VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2, raster = FALSE) + 
  NoLegend()
ggsave(filename = "Vlnplot1.pdf", p1, width = length(unique(sce.all$orig.ident))/3+5, height = 5, dpi = 600)
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2 <- VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims = T) + 
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  NoLegend()
ggsave(filename = "Vlnplot2.pdf", p2, width = length(unique(sce.all$orig.ident))/3+5, height = 5, dpi = 600)
#
feats <- c("nFeature_RNA", "nCount_RNA","percent_mito")
p1 <- VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
ggsave(filename = "Vlnplot_raw.pdf", p1, width = length(unique(sce.all$orig.ident))/3+5, height = 5, dpi = 600)
# 
p3.feature <- FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p3.mito <- FeatureScatter(sce.all, "nCount_RNA", "percent_mito", group.by = "orig.ident", pt.size = 0.5)
# p3.rb <- FeatureScatter(sce.all, "nCount_RNA", "percent_ribo", group.by = "orig.ident", pt.size = 0.5)
# p3.rb.mito <- FeatureScatter(sce.all, "percent_ribo", "percent_mito", group.by = "orig.ident", pt.size = 0.5)
p3 <- p3.feature + p3.mito #+ p3.rb + p3.rb.mito
# p31 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent_mito")
# p32 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent_ribo")
# p33 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent_hb")
# p34 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# p3 <- p31 + p32 + p33 + p34
ggsave(filename = "Scatterplot.pdf", p3, width = 16, height = 6, dpi = 600)

# nCount_RNA和nFeature_RNA应该是正相关，越高越好，否则说明存在低质量细胞或者双细胞
# percent_mito应该与nCount_RNA不相关，不然测的越多，比例越多，说明细胞都死了

### 过滤一下
# 常见的标准是，保留线粒体基因比率<10%，200 < nFeature_RNA < 5000 or 6000
# 如果质控更加严格的话，可以保留线粒体基因比率<5%, 200 < nFeature_RNA < 2500
# 过滤指标:最少表达基因数的细胞&最少表达细胞数的基因
# 说实话不过滤也可以
# 这里按照原作者来 http://links.lww.com/HEP/I90
sce.all.filt <- subset(sce.all, subset = nFeature_RNA > 500 & percent_mito < 30)
dim(sce.all.filt)
# 27371 52784

#可视化过滤后的情况
feats <- c("nFeature_RNA", "nCount_RNA","percent_mito")
p1_filtered <- VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
ggsave(filename = "Vlnplot1_filtered.pdf", p1_filtered, 
       width = length(unique(sce.all.filt$orig.ident))/3+5, height = 5, dpi = 600)
# 警告是因为Seurat v5：计数矩阵存储在pbmc[["RNA"]]$counts。

### 然后作者是用SCT = NormalizeData + FindVariableFeatures + ScaleData，SCT感觉一点都不好用。
# 对线粒体基因表达、UMI计数和检测到的基因的检测率进行回归
library(glmGamPoi)# 说是可以加速，去除变异的混淆源
sce.all.filt <- SCTransform(sce.all.filt, 
                            variable.features.n = 3000, 
                            vars.to.regress = c("percent_mito", "nCount_RNA", "nFeature_RNA"),
                            vst.flavor = 'v2')# 默认V2版本
# # top10 高变基因
# top10 <- head(VariableFeatures(sce.all.filt), 10)
#### 不知道为什么会报错
# Error in SCTResults(object = object, slot = "feature.attributes")[, vars] : 
#   量度数目不对
# p1_VF <- VariableFeaturePlot(sce.all.filt)
# LabelPoints(plot = p1_VF, points = top10, repel = TRUE)

## 
sce.all.filt <- RunPCA(sce.all.filt, features = VariableFeatures(object = sce.all.filt))# 用高变基因作为PCA
VizDimLoadings(sce.all.filt, dims = 1:2, reduction = "pca")
DimPlot(sce.all.filt, reduction = "pca", group.by = "orig.ident")
DimHeatmap(sce.all.filt, dims = 1:10, cells = 500, balanced = TRUE)
# 选PC
ElbowPlot(sce.all.filt, ndims = 30)

# 考虑用一下哈佛研究团队方法确定最佳PC数
# Determine percent of variation associated with each PC
pct <- sce.all.filt[["pca"]]@stdev / sum(sce.all.filt[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
ggsave(filename = "确定最佳PC数.pdf", width = 8, height = 6, dpi = 600)
## PC选18，那选20即可
qsave(sce.all.filt, 'sce.all.filt.qs')
################ 去除批次效应
sce.all.filt <- qread("sce.all.filt.qs")
##可视化
sce.all.filt <- RunUMAP(sce.all.filt, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(sce.all.filt, reduction = "umap.unintegrated", group.by = "orig.ident")+
  NoAxes()+ ggtitle("raw_data")

#### 去批次，踩坑了，如果我前面直接把layer直接合并了，后面就不能再用去除批次了
sce.all.int <- IntegrateLayers(object = sce.all.filt, method = HarmonyIntegration,
                               orig.reduction = "pca", new.reduction = "harmony",
                               normalization.method = "SCT", verbose = T)# 
Layers(sce.all.int)
sce.all.int <- FindNeighbors(sce.all.int, dims = 1:20, 
                              k.param = 20, # Defines k for the k-nearest neighbor algorithm
                              # 增大 k.param 可以引入更多的细胞节点进行聚类
                              prune.SNN = 1/15, reduction = "harmony")
# 在计算SNN构造的邻域重叠时，设置可接受的Jaccard索引的截止值。任何值小于或等于此值的边都将被设置为0并从SNN图中删除。
##以0.5分辨率对细胞聚类，计算构建KNN和SNN图
# 1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 
# 3 = SLM algorithm; 4 = Leiden algorithm；
sce.all.int <- FindClusters(sce.all.int, resolution = 0.3, algorithm = 1)# 作者说直接用0.3
# Leiden算法使聚类结构更加稳定和精细，Leiden requires the leidenalg python.
# 设置不同的分辨率，观察分群效果(resolution越大分群越多越细)
sce.all.int <- RunUMAP(sce.all.int, dims = 1:20, reduction = "harmony")
p.compare <- cowplot::plot_grid(ncol = 2,
                                DimPlot(sce.all.filt, reduction = "umap.unintegrated", group.by = "orig.ident")+
                                  NoAxes()+ ggtitle("raw data"),
                                DimPlot(sce.all.int, reduction = "umap", group.by = "orig.ident")+
                                  NoAxes()+ ggtitle("harmony integrated"))
ggsave(filename = "Before&After_int_harmony.pdf", p.compare, width = 14, height = 6, dpi = 600)
qsave(sce.all.int, "sce.all.int.qs")

############# 2-cluster
dir.create("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HCC/2-cluster")
setwd("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HCC/2-cluster")


