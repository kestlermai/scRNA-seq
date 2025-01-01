library(Seurat)
library(gplots)
library(ggplot2)

setwd("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8")
sce.all <- readRDS("sce.all_marker_cluster.rds")

target_genes <- c("MX1","IFIT3","GBP1","SERPING1","IFIT1","IFI44","ISG15","IFI27")
target_genes <- c("ISG15","IFI27")

p1 <- VlnPlot(sce.all, features = target_genes, group.by = "group", pt.size = 0.5, ncol = 2)+
  NoLegend()
p1
ggsave(filename="target_genes_group.pdf", p1, width = 8, height = 4, dpi = 600)

p2 <- VlnPlot(sce.all, features = target_genes, pt.size = 0.5, ncol = 2)+
  NoLegend()
p2
ggsave(filename="target_genes_cell.pdf", p2, width = 8, height = 6, dpi = 600)

p3 <- FeaturePlot(sce.all, reduction = "umap", features = target_genes, pt.size = 0.5, ncol = 4,
              label = T)+
  NoLegend()
p3
ggsave(filename="target_genes_umap.pdf", p3, width = 16, height = 8, dpi = 600)

p4 <- VlnPlot(sce.all, features = target_genes,
                 split.by = "group", group.by = "celltype",
                 pt.size = 0.5, combine = T, ncol = 2)+
  theme(legend.position = "top")
p4
ggsave(filename="target_genes_group_celltype.pdf", p4, width = 10, height = 6, dpi = 600)
