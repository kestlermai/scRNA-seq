if(!require(Seurat))BiocManager::install("Seurat")
if(!require(harmony))devtools::install_github("immunogenomics/harmony")
if(!require(ggplot2))BiocManager::install("ggplot2")
if(!require(dplyr))BiocManager::install("dplyr")
if(!require(tidyverse))BiocManager::install("tidyverse")
if(!require(patchwork))BiocManager::install("patchwork")


setwd("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV")

#提取文件格式
#fileformat <-  unlist(lapply(fs,function(x){strsplit(x,'_')[[1]][3]}))
#这个数据只有一个_符号分隔，上面代码不适用
#批量对名字改名GSM4775588_C1barcodes.tsv.gz改成GSM4775588_C1_barcodes.tsv.gz
# 设置文件夹路径
folder_path <- "C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/raw"
# 获取文件夹中的所有文件名
file_names <- list.files(folder_path)
file_names
# 定义正则表达式模式
pattern <- "(C1|Q[1-7])(.*)"
# 批量改名
for (file_name in file_names) {
  new_file_name <- sub(pattern, "\\1_\\2", file_name)
  old_file_path <- file.path(folder_path, file_name)
  new_file_path <- file.path(folder_path, new_file_name)
  file.rename(old_file_path, new_file_path)
}
# 输出改名后的文件名
new_file_names <- list.files(folder_path)
print(new_file_names)

#提取文件名
filename <- list.files("./raw")
filename
#提取样品名
sample <- unlist(lapply(filename,function(x){strsplit(x,'_')[[1]][1]}))
sample
#提取文件格式
fileformat <-  unlist(lapply(filename,function(x){strsplit(x,'_')[[1]][3]}))
fileformat

for(i in 1: length(filename)){
  #创建每个样本自己的文件夹:
  dir.create(paste0('raw/',sample[i]))
  #把文件依次移动到新的文件夹中并命名:
  file.copy(paste0('raw/',filename[i]),
            paste0('raw/',sample[i],'/',fileformat[i]))
  #删除原文件:      
  unlink(paste0('raw/',filename[i]))
}

#构建路径
samplenames <- unique(sample)
filepath <- paste0("raw/", samplenames)
filepath
list.files(filepath)
names(filepath) <- samplenames
filepath
#CellRanger < 3.0 文件格式必须是barcodes.tsv, features.tsv, and matrix.mtx
#CellRanger ＞= 3.0，文件格式必须是barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
#批量读取
dir <- "C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/raw"
scList <- list()
for (sample in samplenames) {
  #读取文件
  scrna_data <- CreateSeuratObject(
    counts = Read10X(file.path(dir,sample)# ,gene.column=1 有些三联数据不符合规则
    ),
    project = sample,
    min.cells = 3, # 去除小于3个细胞中表达的基因
    min.features = 200  # 去除只有200以下基因表达的细胞
  )
  # metadata 添加一列
  scrna_data@meta.data$sample <- sample
  # 将 seurat 对象存入列表
  scList[[sample]] = scrna_data 
}

#合并数据
sc_merge <- merge(x = scList[[1]],
              y = scList[-1],
              add.cell.ids = names(scList), project = "scRNA-HIV")
sc_merge
dim(sc_merge)
table(sc_merge@meta.data$orig.ident)
# 统计细胞表达分布
hist(colSums(sc_merge),
     breaks = 150, main = "UMI count per cell",
     xlab = "UMI count per cell")

#QC
### 主要PercentageFeatureSet函数计算线粒体含量
### 人类使用pattern = "^MT-"，小鼠使用pattern = "^mt-"
#线粒体
sc_merge[["percent.mt"]] <- PercentageFeatureSet(sc_merge, pattern = "^MT-")
#核糖体 
#sc_merge[["percent.rp"]] <- PercentageFeatureSet(sc_merge, pattern = "^RP")
#提取细胞信息 meta.data
head(sc_merge@meta.data)
load()
### 质控数据可视化，使用VlnPlot函数
### nFeature_RNA,  每个细胞中有多少个基因
### nCount_RNA, 每个细胞中有多少个counts
### percent.mt, 每个细胞中线粒体基因的比例
VlnPlot(sc_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)
nFeature_RNA <- sc_merge@meta.data$nFeature_RNA
nCount_RNA <- sc_merge@meta.data$nCount_RNA
percent.mt <- sc_merge@meta.data$percent.mt
p <- ggplot(df,aes(x = x, y = y))+
  geom_jjviolin(width = 0.1,
                trim = F,
                type = 'right',
                shift = 0.005,
                position = position_dodge(width = 0.8)) +
  geom_jjboxplot(width = 0.3,
                 type = 'left',
                 shift = 0.005,
                 position = position_dodge(width = 0.8)) +
  theme_prism()
p
### 正式筛选，筛选的是细胞，最终细胞减少
### nFeature_RNA > 200
### nFeature_RNA < 2500 or 3000
### percent.mt < 5
sc_merge <- subset(sc_merge, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# 将QC指标可视化为QC后的小提琴曲线图
VlnPlot(sc_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
## Normalizing the data
sc_merge <- NormalizeData(sc_merge,normalization.method = "LogNormalize", scale.factor = 10000)
#Identification of highly variable features (feature selection)
sc_merge <- FindVariableFeatures(sc_merge, selection.method = "vst", nfeatures = 2000)#3000
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sc_merge), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sc_merge)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#线性缩放，PCA之前需要标准预处理步骤，归一化
sc_merge <- ScaleData(sc_merge, features = rownames(sc_merge))
#线性降维
sc_merge <- RunPCA(sc_merge, features = VariableFeatures(object = sc_merge), reduction.name = "pca")
print(sc_merge[["pca"]], dims = 1:5, nfeatures = 5)
#可视化
DimPlot(sc_merge, reduction = "pca")
DimHeatmap(sc_merge, dims = 1:6, cells = 500, balanced = TRUE)

#维数选择
sc_merge <- JackStraw(sc_merge, num.replicate = 100)
sc_merge <- ScoreJackStraw(sc_merge, dims = 1:20)

JackStrawPlot(sc_merge, dims = 1:20)
ElbowPlot(sc_merge)

#细胞聚类--基于图的聚类算法
#FindNeighbors中的dims指定聚类使用的维数; 
#FindClusters 中的 resolution 指定类别的精度，越大则分出越多的类；越小则类别越少
sc_merge <- FindNeighbors(sc_merge, dims = 1:10)
sc_merge <- FindClusters(sc_merge, resolution = 0.5)

##非线性降维(UMAP/tSNE)
### 不进行批次矫正
sc_merge <- RunUMAP(sc_merge, dims = 1:10)
DimPlot(sc_merge, reduction = "umap", label = TRUE)
sc_merge <- RunTSNE(sc_merge, dims = 1:10)
DimPlot(sc_merge, reduction = "tsne", label = TRUE)

sc_merge <- RunUMAP(sc_merge,reduction = "pca", dims = 1:10, reduction.name = "umap_naive")
DimPlot(sc_merge, reduction = "umap")

### 使用harmony进行批次矫正
### 默认reduction = "pca",group.by.vars参数输入的是批次信息,reduction.save 是结果保存的名称
sc_merge <- RunHarmony(sc_merge,reduction = "pca",group.by.vars = "group",reduction.save = "harmony")
sc_merge <- RunUMAP(sc_merge,reduction = "harmony",dims = 1:30,reduction.name = "umap")

p1 <- DimPlot(sc_merge, reduction = "umap_naive",group.by = "group")
p2 <- DimPlot(sc_merge, reduction = "umap",group.by = "group")
p1+p2
