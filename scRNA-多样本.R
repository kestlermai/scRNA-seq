if(!require(Seurat))BiocManager::install("Seurat")
if(!require(harmony))devtools::install_github("immunogenomics/harmony")
if(!require(ggplot2))BiocManager::install("ggplot2")
if(!require(dplyr))BiocManager::install("dplyr")
if(!require(patchwork))BiocManager::install("patchwork")

dirs = dir(pattern = "^R")
f = "dat.Rdata"
if(!file.exists(f)){
  scelist = list()
  for(i in 1:length(dirs)){
    x = Read10X(data.dir = dirs[[i]])
    scelist[[i]] <- CreateSeuratObject(counts = x, 
                                       project = paste0("R",i))
    scelist[[i]][["percent.mt"]] <- PercentageFeatureSet(scelist[[i]], pattern = "^MT-")##计算线粒体含量
    scelist[[i]][["percent.rp"]] <- PercentageFeatureSet(scelist[[i]], pattern = "^RP-")##计算核糖体含量
    scelist[[i]] <- subset(scelist[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
  }
  names(scelist)  = paste0("R",1:4)
  sum(sapply(scelist, function(x)ncol(x@assays$RNA@counts)))
  # normalize and identify variable features for each dataset independently
  scelist <- lapply(X = scelist, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  #整合每个数据集中重复的高变基因
  features <- SelectIntegrationFeatures(object.list = scelist)
  #整合高变基因找出两两数据集间的锚点细胞
  immune.anchors <- FindIntegrationAnchors(object.list = scelist, anchor.features = features)
  immune.combined <- IntegrateData(anchorset = immune.anchors)
  DefaultAssay(immune.combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  immune.combined <- ScaleData(immune.combined, verbose = FALSE)
  immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
  immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
  immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
  immune.combined <- FindClusters(immune.combined, resolution = 0.5)
  save(immune.combined,file = f)
}
load(f)
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
