# UMAP可视化时候总是出现毛毛躁躁的边缘和大量散在细胞，还有很多细胞亚群之间有连续的细胞
# 可能时真是存在的过渡态细胞，不然需要考虑这个用SoupX或者DecontX预测和去除单细胞转录组的环境游离RNA污染：


# SoupX -----------------------------------------------------------------
# 所有基于液滴的单细胞 RNA-seq 实验还捕获了输入溶液中存在的环境 mRNA 以及感兴趣的细胞特异性 mRNA。
# 这种污染无处不在，并且在不同实验之间可能会有很大差异 （2% - 50%），尽管大约 10% 似乎相当普遍。
# 没有办法提前知道实验中的污染是什么，尽管实体瘤和低活力细胞往往会产生更高的污染分数。
# 由于污染 mRNA 的来源是输入溶液中的裂解细胞，因此污染的概况是特定于实验的，并产生批次效应。
# 即使您出于任何原因决定不想使用 SoupX 校正方法，您至少应该想知道您的数据受到的污染程度。

# 在基于液滴的scRNA-seq实验中，稀释液中始终存在一定量的背景mRNA，这些mRNA与细胞一起分布到液滴中，
# 并与它们一起测序。这样做的净效应是产生背景污染，该污染代表的表达不是来自液滴中包含的细胞，
# 而是来自包含细胞的溶液。
# 这种漂浮在输入溶液中的游离mRNA集合（作者这里称为soup）是由输入溶液中被裂解的细胞产生的。
# 因此，每个输入溶液的soup看起来都不同，并且与通过对所有单个细胞求和获得的表达模式非常相似。
# SoupX的目的是提供一种方法来估计这种soup的组成，每个液滴中从soup中得出的 UMI 分数，
# 并生成一个校正的计数表，其中去除了基于soup的表达式。

# 分为三部分
# 1.计算soup的profile：计算估计空液滴中RNA表达量
# 2.估计细胞特异性污染分数：每个细胞受到污染的比例有多少
# 3.推断校正后的表达式矩阵：矫正矩阵

# devtools::install_github("constantAmateur/SoupX", ref = 'devel')
library(SoupX)
library(Seurat)

## 这里demo还是用的是10X PBMC 4k
tmpDir <- tempdir(check = TRUE)
download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz", 
              destfile = file.path(tmpDir, "tod.tar.gz"))
download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz", 
              destfile = file.path(tmpDir, "toc.tar.gz"))
# 解压
untar(file.path(tmpDir, "tod.tar.gz"), exdir = tmpDir)
untar(file.path(tmpDir, "toc.tar.gz"), exdir = tmpDir)

## 读取数据
sc <- load10X(tmpDir)
# 读取过滤后的seurat对象--toc：table of counts
toc <- Seurat::Read10X(file.path(tmpDir, "filtered_gene_bc_matrices", "GRCh38")) # 路径
# 读取raw—seurat对象--tod：Table of droplets 滴液raw矩阵
tod <- Seurat::Read10X(file.path(tmpDir, "raw_gene_bc_matrices", "GRCh38"))
# 创建一个SoupChannel对象，包括raw和filter后
sc <- SoupChannel(tod = tod, toc = toc, calcSoupProfile = TRUE) # 使用默认值的estimateSaup计算
sc # Channel with 33694 genes and 4340 cells
# estimateSaup默认计算soupRange = c(0, 100)，当然也可以自己调整范围，后面再详细说
### 不想下载可以直接用包自带的PBMC数据
data(PBMC_sc)
sc <- PBMC_sc
sc # Channel with 33694 genes and 2170 cells
# 这里其实都是用的是PBMC 4k数据，但细胞量只有一半，源代码是这样写的：
# PBMC_sc <- SoupChannel(PBMC_sc$tod,PBMC_sc$toc[,sample(ncol(PBMC_sc$toc),round(ncol(PBMC_sc$toc)*0.5))])

### 如果要对数据进行插补，应该是先跑SoupX，先质控再插补
# 而且需要用adjustCounts函数设置：roundToInt = T

#### 如果只有filter后的矩阵，没有空滴液
# library(Matrix)
# toc <- sc$toc
# scNoDrops <- SoupChannel(toc, toc, calcSoupProfile = FALSE)
# # Calculate soup profile
# soupProf <- data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))
# scNoDrops <- setSoupProfile(scNoDrops, soupProf)

#### 向SoupChannel对象添加额外的metadata，因为SoupX可以用上cluster信息，通过cluster更好地自动估计污染分数，
# 作者在github说基于图形的聚类（应该说的是Louvain和Leiden聚类算法），效果也足够好，
# 加上celltype注释可能效果更好
# 这里sc用PBMC_sc，不能用完整的PBMC 4k，metadata对应不上（应该是没有处理过的）
data(PBMC_metaData)
# 导入细胞注释信息
sc <- setClusters(sc, setNames(PBMC_metaData$Cluster, rownames(PBMC_metaData)))
# 降维 DR: dimension reduction，降维好像默认是tsne
sc <- setDR(sc, DR = PBMC_metaData[colnames(sc$toc), c("RD1", "RD2")]) # 降维坐标

library(ggplot2)
dd <- PBMC_metaData[colnames(sc$toc), ]
mids <- aggregate(cbind(RD1, RD2) ~ Annotation, data = dd, FUN = mean)
gg <- ggplot(dd, aes(RD1, RD2)) + 
  geom_point(aes(colour = Annotation), size = 0.2) + 
  geom_label(data = mids, aes(label = Annotation)) + 
  ggtitle("PBMC 4k Annotation") + 
  guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

### 我们需要了解一组细胞中一个基因（或一组基因）的表达是否来自soup
# 假设我们对基因 IGKC 的表达感兴趣，IGKC 是 B 细胞高度表达的关键成分免疫球蛋白（即抗体）。
dd$IGKC <- sc$toc["IGKC", ]
gg <- ggplot(dd, aes(RD1, RD2)) + 
  geom_point(aes(colour = IGKC > 0))
plot(gg)
# 右下角很明显是B细胞，但B细胞边上有一部分T细胞，CD4⁺ T细胞激活B细胞分化成浆细胞从而产生抗体
# 检查一下这些分散细胞中 IGKC 的表达是否超出了我们从soup中随机预期的水平，

# plotMarkerMap计算观察到的计数与预期计数的绘图比，也就是预计的soup与实际观察的count相比较（泊松分布）
gg <- plotMarkerMap(sc, "IGKC")
plot(gg)
# 我们看到 B 细胞簇中的细胞呈现红色，表明它们的表达远超过我们随机期望的水平
# 我们改变范式、产生抗体的 T 细胞表现不佳。它们都呈现出明显的蓝色阴影，
# 这表明这些细胞中 IGKC 的表达完全有可能是由于混合液的污染
# 上面只是假设每个液滴只包含背景污染，但这显然不真实

## 估计细胞特异性污染分数：手动和自动，作者推荐自动估计
# 识别在我们的数据中某些细胞未表达的基因，以及我们观察到这些基因在这些细胞中的表达必须是由于污染

# 手动指定污染分数--作者推荐20%-可去除99%
sc <- setContaminationFraction(sc, 0.2)

# decontX -----------------------------------------------------------------

# 基于液滴的微流控装置（主要是10X Genomics）已被广泛用于单细胞 RNA 测序(scRNA-seq)。
# 然而，存在于细胞悬浮液中的环境 RNA 可以与细胞的天然 mRNA 一起被异常计数，
# 并导致不同细胞群之间转录物的交叉污染。DecontX 是一种新的贝叶斯方法，用于估计和去除单个细胞中的污染。

devtools::install_github("campbio/decontX")
