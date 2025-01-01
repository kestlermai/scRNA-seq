if(!require(Seurat))BiocManager::install("Seurat")
if(!require(ggplot2))BiocManager::install("ggplot2")
if(!require(dplyr))BiocManager::install("dplyr")
if(!require(patchwork))BiocManager::install("patchwork")

data_dir <- 'G:/maihuanzhuo/zhang'
list.files(data_dir)

sc_data <- Read10X(data.dir = data_dir, gene.column = 1)

sc <- CreateSeuratObject(sc_data, project = 'scRNA', min.features = 200)

#计算线粒体基因比例
mito_genes <- rownames(sc)[grep("^MT-", rownames(sc))]
mito_genes #13个线粒体基因
sc <- PercentageFeatureSet(sc, "^MT-", col.name = "percent_mito")
fivenum(sc@meta.data$percent_mito)

feats <- c("nFeature_RNA", "nCount_RNA","percent_mito")
p1 <- VlnPlot(sc, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
ggsave(filename="Vlnplot_raw.pdf", p1, width = length(unique(sc$orig.ident))/3+5, height = 5, dpi = 600)


p2 <- FeatureScatter(sc, "nFeature_RNA", "percent_mito", group.by = "orig.ident", pt.size = 0.5)
ggsave(filename="Scatterplot_mito.pdf", p2, width = 10, height = 6, dpi = 600)

#filter mito
sc_selected <- WhichCells(sc, expression = nFeature_RNA > 1000 & percent_mito < 25)
sc.filt <- subset(sc, cells = sc_selected)
dim(sc.filt)#33694 119438
dim(sc)#33694 212566

#标化
sc.filt <- SCTransform(sc.filt)
# #筛选top 3000高变基因
# sc.filt <- FindVariableFeatures(sc.filt, selection.method = "vst", nfeatures = 3000)
# #缩放数据
# sc.filt <- ScaleData(sc.filt, features = rownames(sc.filt))
#PCA降维
sc.filt <- RunPCA(sc.filt)

####作者run PCA出来的dims选了20，为了省时间就直接用20了，不过我看他elbow plot用的是30
#
sc.filt <- RunUMAP(sc.filt, dims = 1:20)
sc.filt <- FindNeighbors(sc.filt, dims = 1:20)
sc.filt <- FindClusters(sc.filt, resolution = 0.01)#分辨率这么低吗？

p1 <- DimPlot(sc.filt, reduction = "umap", group.by = "seurat_clusters", label = T) 
ggsave("umap_by_seurat_clusters.pdf", p1, width = 8, height = 7, dpi = 600)

markers <- FindAllMarkers(sc.filt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

markerGenes <- c('PTPRC',#immune cells
                 'EPCAM',#epithelial cells
                 'PECAM1'#endothelial cells
)
p2 <- VlnPlot(sc.filt, features = markerGenes, pt.size = 0, stack = T)+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 12))+
  NoLegend()
ggsave("VlnPlot_markerGenes-stack.pdf", p2, width = 18, height = 12, dpi = 600)

Idents(sc.filt, cells = WhichCells(sc.filt, idents = c(0, 3))) <- "Immune"
Idents(sc.filt, cells = WhichCells(sc.filt, idents = c(1, 2))) <- "Epithelial"
Idents(sc.filt, cells = WhichCells(sc.filt, idents = c(4))) <- "Endothelial"
Idents(sc.filt, cells = WhichCells(sc.filt, idents = c(5))) <- "Mesenchymal"
sc.filt@meta.data$population <- Idents(sc.filt)

p3 <- DimPlot(sc.filt, reduction = 'umap', group.by = 'population', label = T, pt.size = 2)
p3
ggsave("UMAP.pdf", p3, width = 18, height = 12, dpi = 600)
saveRDS(sc.filt,'G:/maihuanzhuo/zhang/zhang.rds')

#读取过滤后的数据
sc.filt <- readRDS('G:/maihuanzhuo/zhang/zhang.rds')
# Subset into 4 major populations
immune <- subset(sc.filt, idents = "Immune")
epi <- subset(sc.filt, idents = "Epithelial")
endo <- subset(sc.filt, idents = "Endothelial")
meso <- subset(sc.filt, idents = "Mesenchymal")
# Subset into 4 major populations
table(immune@meta.data$population)#63058
table(epi@meta.data$population)#39280
table(endo@meta.data$population)#11299
table(meso@meta.data$population)#5801

# Immune subset
immune <- FindVariableFeatures(immune, nfeatures = 3000)
immune <- ScaleData(immune, features = row.names(immune@assays$SCT@data))
immune <- RunPCA(immune)
#上面这一步其实就是确定最佳dims数，好让后面分群
#按照作者确定的dims=20跑
immune <- RunUMAP(immune, dims = 1:21)
immune <- FindNeighbors(immune, dims = 1:20)
immune <- FindClusters(immune, resolution = 1)

immune_marker_genes <- c('MARCO',#Macrophages
                         'CD3E','CD8A','IL7R','CD8B','CD27',#T cells
                         'MS4A1','CD79A',#B cells
                         'NKG7','CD160','NCR1',# NK cells
                         'LYZ','S100A12','CD14',# monocytes
                         "LILRA4", "CLEC4C", "JCHAIN", "TCF4",#pDCs
                         "FCER1A", "CD1C", "CLEC9A", "FSCN1",#cDCs
                         'CPA3','KIT',#Mast cells
                         "JCHAIN", "IGHG1"#Plasma Cells
)
p1 <- DimPlot(immune,reduction = "umap", group.by = "seurat_clusters", label = T)
p1
p2 <- VlnPlot(immune, features = immune_marker_genes, pt.size = 0, stack = T)+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 12))+
  NoLegend()
ggsave("G:/maihuanzhuo/zhang/vlnplot_immune_marker.pdf", p2, width = 18, height = 12, dpi = 600)


Idents(immune, cells = WhichCells(immune, idents = c(25))) <- "B Cells"
Idents(immune, cells = WhichCells(immune, idents = c(29))) <- "Plasma Cells"
Idents(immune, cells = WhichCells(immune, idents = c(22,10,15,8,20,12,30))) <- "Monocytes"
Idents(immune, cells = WhichCells(immune, idents = c(19))) <- "Monocytes"
Idents(immune, cells = WhichCells(immune, idents = c(3,9,23,32,31,35))) <- "T Cells"
Idents(immune, cells = WhichCells(immune, idents = c(13))) <- "NK Cells"
Idents(immune, cells = WhichCells(immune, idents = c(27))) <- "Mast Cells"
Idents(immune, cells = WhichCells(immune, idents = c(36))) <- "pDCs/cDCs"
Idents(immune, cells = WhichCells(immune, idents = c(21,26,7,24,28,4,5,33,0,16,
                                                     34,11,6,17,18,14,1,2))) <- "Macrophages"

immune$celltype <- Idents(immune)
p3 <- DimPlot(immune)
p3
