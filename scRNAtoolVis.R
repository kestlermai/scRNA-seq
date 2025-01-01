# devtools::install_github('junjunlab/scRNAtoolVis')
library(ggplot2)
library(Seurat)
library(scRNAtoolVis)
setwd("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8")
sce.all <- readRDS("sce.all_marker.rds")
table(Idents(sce.all))
# find markers for every cluster compared to all remaining cells
markers <- FindAllMarkers(sce.all,only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0)
#基于 marker 基因的火山图
# 绘制火山图            
p1 <- markerVocalno(
  markers = markers,            
  topn = 5,            
  labelCol = ggsci::pal_npg()(9)            
)
p1
ggsave("marker-Vocalno.pdf",p1,width = 14, height = 8,dpi = 600)
#基于 marker 基因的曼哈顿图
#默认绘制
p2 <- jjVolcano(diffData = markers, topGeneN = 2, col.type = "adjustP")
p2
ggsave("Manhattan-Vocalno.pdf",p2,width = 14, height = 8,dpi = 600)
#指定展现marker基因  
jjVolcano(        
  diffData = markers,
  myMarkers = c('LTB','CD79B','CCR7','GNLY'), 
  log2FC.cutoff = 0.5, # X 轴的细胞注释的高度      
  base_size = 16, # Y 轴的标签及标题的大小
  aesCol = c('purple','orange'),# 上下两边的颜色
  tile.col = corrplot::COL2('RdBu', 16)[1:16],# X 轴的填充色
  size = 3, # marker 基因标签的大小
  cluster.order = rev(unique(pbmc.markers$cluster))
)

#弦图
jjVolcano(diffData = pbmc.markers,              
          topGeneN = 3,              
          tile.col = corrplot::COL2('RdYlBu', 15)[4:12],              
          size  = 3.5,              
          fontface = 'italic',              
          polar = T) +              
  ylim(-8,10)