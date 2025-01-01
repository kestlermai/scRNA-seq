#scRNA seurat violin plot 优化
library(ggplot2)
library(Seurat)
library(ggprism)
library(ggpubr)
#
fibro <- readRDS("C:/Users/maihuanzhuo/Desktop/fibro.rds")
fibro
dim(fibro)
table(Idents(fibro))
table(fibro$celltype)
#
p1 <- VlnPlot(subset(fibro, CTSK > 0), features = 'CTSK', 
              adjust = 0.5, cols = c("#385989","#d22027"),
              split.by = "group")+
  theme(legend.position = 'top',legend.direction = "horizontal",legend.justification = "center")
  
p1
ggsave(filename="C:/Users/maihuanzhuo/Desktop/fibro_VlnPlot_CTSK.pdf", p1, width = 6, height = 5, dpi = 600)
#
feature_genes <- c('CILP','COL14A1','CPXM1','CTSK','GREM1','HMGCS2','LEMD1','MMP12','NPR3')
p2 <- VlnPlot(subset(fibro, CTSK > 0), features = feature_genes, 
              adjust = 0.5, cols = c("#385989","#d22027"),
              split.by = "group",
              combine = T, ncol = 3)+
  NoLegend()
ggsave(filename="C:/Users/maihuanzhuo/Desktop/fibro_VlnPlot.pdf", p2, width = 14, height = 12, dpi = 600)
#
macr <- readRDS("C:/Users/maihuanzhuo/Desktop/macrophage.rds")
#
p3 <- VlnPlot(subset(macr, CTSK > 0), features = 'CTSK', 
              adjust = 0.5, cols = c("#385989","#d22027"),
              split.by = "group")+
  theme(legend.position = 'top',legend.direction = "horizontal",legend.justification = "center")

p3
ggsave(filename="C:/Users/maihuanzhuo/Desktop/macr_VlnPlot_CTSK.pdf", p3, width = 6, height = 5, dpi = 600)
#
p4 <- VlnPlot(subset(macr, CTSK > 0), features = feature_genes, 
              adjust = 0.5, cols = c("#385989","#d22027"),
              split.by = "group",
              combine = T, ncol = 3)+
  NoLegend()
ggsave(filename="C:/Users/maihuanzhuo/Desktop/macr_VlnPlot.pdf", p4, width = 14, height = 12, dpi = 600)
#
data <- data.frame(p1[[1]][["data"]])
p1_box <- ggplot(data,aes(x = ident, y = CTSK, fill = split))+
  geom_boxplot(outlier.size = 1, outlier.colour = "black", 
              fill = c("#385989","#d22027"))+
  geom_pwc(aes(group = split), method = "t_test",
           label.size = 5, size = 1, family = "sans",
           #label = "italic(p)= {p}",
           tip.length = 0,bracket.nudge.y = 0.1#调整整体label的位置
  )+
  theme_prism()+
  labs(x = "Identity", y = "Expression Level", title = "CTSK")+
  theme(legend.position = 'top',legend.direction = "horizontal",legend.justification = "center",
        axis.text = element_text(color="black", family = "sans", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p1_box
ggsave(filename="C:/Users/maihuanzhuo/Desktop/fibro_BoxPlot.pdf", p1_box, width = 6, height = 5, dpi = 600)
#
data <- data.frame(p3[[1]][["data"]])
p3_box <- ggplot(data,aes(x = ident, y = CTSK, fill = split))+
  geom_boxplot(outlier.size = 1, outlier.colour = "black", 
               fill = c("#385989","#d22027"))+
  geom_pwc(aes(group = split), method = "t_test",
           label.size = 5, size = 1, family = "sans",
           #label = "italic(p)= {p}",
           tip.length = 0,bracket.nudge.y = 0.1#调整整体label的位置
           )+
  theme_prism()+
  labs(x = "Identity", y = "Expression Level", title = "CTSK")+
  theme(legend.position = 'top',legend.direction = "horizontal",legend.justification = "center",
        axis.text = element_text(color="black", family = "sans", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p3_box
ggsave(filename="C:/Users/maihuanzhuo/Desktop/macr_BoxPlot.pdf", p3_box, width = 6, height = 5, dpi = 600)
#
Myof_sc <- fibro[,fibro@meta.data$celltype %in% 'Myofibroblasts']
p_Myof <- VlnPlot(subset(Myof_sc, CTSK > 0), features = 'CTSK', 
              adjust = 0.5, cols = c("#385989","#d22027"),
              split.by = "group")+
  labs(title = 'Myofibroblasts_CTSK')+
  theme(legend.position = 'top',legend.direction = "horizontal",legend.justification = "center")

p_Myof
ggsave(filename="C:/Users/maihuanzhuo/Desktop/Myof_VlnPlot_CTSK.pdf", p_Myof, width = 6, height = 5, dpi = 600)
#
data <- data.frame(p_Myof[[1]][["data"]])
p_Myof_box <- ggplot(data,aes(x = ident, y = CTSK, fill = split))+
  geom_boxplot(outlier.size = 1, outlier.colour = "black", 
               fill = c("#385989","#d22027"))+
  geom_pwc(aes(group = split), method = "t_test",
           label.size = 5, size = 1, family = "sans",
           #label = "italic(p)= {p}",
           tip.length = 0,bracket.nudge.y = 0.1#调整整体label的位置
  )+
  theme_prism()+
  labs(x = "Identity", y = "Expression Level", title = "Myofibroblasts_CTSK")+
  theme(legend.position = 'top',legend.direction = "horizontal",legend.justification = "center",
        axis.text = element_text(color="black", family = "sans", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p_Myof_box
ggsave(filename="C:/Users/maihuanzhuo/Desktop/Myof_BoxPlot.pdf", p_Myof_box, width = 6, height = 5, dpi = 600)