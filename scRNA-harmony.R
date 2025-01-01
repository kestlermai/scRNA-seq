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

#批量计算线粒体和红细胞比例
for(i in 1:length(scList)){
  sc <- scList[[i]]
  # 计算线粒体比例
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  #核糖体 
  sc[["rp_percent"]] <- PercentageFeatureSet(sc, pattern = "^RP-")
  # 计算红细胞比例
  HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB_genes, rownames(sc@assays$RNA))
  HB_genes <- rownames(sc@assays$RNA)[HB_m] 
  HB_genes <- HB_genes[!is.na(HB_genes)] 
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, features = HB_genes) 
  # 将sc赋值给scRNAlist[[i]]
  scList[[i]] <- sc
  # 删除sc
  rm(sc)
}
#查看
scList[[1]]@meta.data[1:5,1:5]

#nFeature_RNA：每个细胞检测表达的基因数目大于300，小于7000；
#mt_percent:每个细胞的线粒体基因表达量占总体基因的比例小于10%；
#批量过滤细胞和线粒体
scList <- lapply(X = scList, FUN = function(x){
  x <- subset(x,subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & mt_percent < 10)})

#合并Seurat对象
scList <- merge(x = scList[[1]],
                y = scList[-1],
                add.cell.ids = names(scList))
## 统计细胞数
table(scList[[]]$orig.ident)

#归一化
#筛选高变基因与降维
scList <- NormalizeData(scList) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)

#harmony整合
scRNA_harmony <- RunHarmony(scList, group.by.vars = "orig.ident")
scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]

# 聚类
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 1)
# umap/tsne降维
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)

# 绘图
umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "orig.ident") 
tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)
# 合并图片
umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,
                                          umap_integrated1,umap_integrated2), ncol = 2)
# 保存图片
ggsave("umap_tsne_integrated.pdf",umap_tsne_integrated, width = 20, height = 15, dpi = 1200)

#保存数据
save(scRNA_harmony,scList,file = "sc-rawdata-harmony.Rdata")
