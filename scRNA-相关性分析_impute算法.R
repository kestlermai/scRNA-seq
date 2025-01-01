rm(list = ls())
library(Seurat)
library(qs)
macr <- qread("C:/Users/maihuanzhuo/Desktop/zhanghuilin/macrophage.qs")
macr
table(Idents(macr))
case_macr <- subset(macr, cells= rownames(macr@meta.data[macr@meta.data$group=="Case",]))
table(Idents(case_macr))
rm(macr)
#直接用seurat包中FeatureScatter函数直接比较，在质控的时候也能看到MT跟Count的相关性
FeatureScatter(sce,feature1 = 'GeneA',feature2 = 'GeneB')
#或者提取表达矩阵然后进行循环cor
exp <- LayerData(case_macr, assay = "RNA", layer = "counts")#V5版
exp <- LayerData(case_macr, assay = "SCT", layer = "counts")
exp <- GetAssayData(case_macr, assay = "RNA", slot = "data")#V4版
class(exp)
# #转换成
exprSet <- as.data.frame(exprSet)
target_gene = 'CTSK'
exp_CTSK <- exp['CTSK',]
# 但有个问题，单细胞表达矩阵是稀疏矩阵，这意味着有很多细胞基因表达量都为零，
# 那么都不用浪费时间去跑都知道结果肯定是不理想的
# 基于马尔可夫亲和力的细胞图插补 (MAGIC) 是一种对高维数据进行去噪的算法，最常应用于单细胞 RNA 测序数据。
# MAGIC 学习流形数据，使用生成的图来平滑特征并恢复数据的结构。

#install.packages("Rmagic"),这个包已经被cran删了，手动去官网找之前的版本装吧
# https://cran.r-project.org/src/contrib/Archive/Rmagic/
# 好像安装这个包要调用python，不懂弄好配置env没
# 教程：https://github.com/KrishnaswamyLab/magic
gc()#释放内存
library(Rmagic)
library(reticulate)

#激活magic的conda环境
reticulate::use_condaenv("magic", required = TRUE)
#标准化处理
exp_magic <- library.size.normalize(exp)
#对归一化后取平方根
exp_magic <- sqrt(exp_magic)


######
# scImpute: accurate and robust imputation of scRNA-seq data
# options(timeout = 1000000)
# devtools::install_github("Vivianstats/scImpute")
# install.packages("C:/Users/maihuanzhuo/Desktop/scImpute-master/",
#                  repos = NULL,type = "source")
# install.packages('penalized')
library(scImpute)

scimpute(# full path to raw count matrix
  count_path = system.file("extdata", "raw_count.csv", package = "scImpute"), 
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = "./",           # full path to output directory
  labeled = FALSE,          # cell type labels not available
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 2,             # 2 cell subpopulations
  ncores = 10)              # number of cores used in parallel computation
