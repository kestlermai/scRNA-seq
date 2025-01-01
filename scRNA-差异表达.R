#最终细胞类型的差异基因 
library(Seurat)
library(gplots)
library(ggplot2)

setwd("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker")
sce.all <- readRDS("sce.all_marker_cluster.rds")
tab.1 <- table(sce.all@meta.data$CCA_snn_res.0.8,sce.all@meta.data$celltype)

#可视化每种细胞类型的基因数和Count数
feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered <- VlnPlot(sce.all, group.by = "celltype", 
                    features = feats, pt.size = 0.1, ncol = 2) + 
  NoLegend()
p1_filtered
ggsave(filename="Vlnplot1_filtered_marker.pdf",p1_filtered, width = 10, height = 8, dpi = 600)
#分别可视化"percent_mito", "percent_ribo", "percent_hb"在不同类型细胞中的比例
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered <- VlnPlot(sce.all, group.by = "celltype", 
                    features = feats, pt.size = 0.1, ncol = 3) + 
  NoLegend()
p2_filtered
ggsave(filename="Vlnplot2_filtered_marker.pdf",p2_filtered, width = 12, height = 8, dpi = 600)

#不同细胞亚群比例变化
#可视化2：所有细胞注释分类的饼图
tb <- as.data.frame(table(phe$celltype),stringsAsFactors = F)
tb <- tb[order(tb$Freq,decreasing=F),]
colnames(tb) <- c('Celltype', 'Freq')
percentage <- scales::percent(tb$Freq / sum(tb$Freq))
percentage
labs <- paste0(tb$Celltype,'(', percentage, ')')#设置标签名
labs
tb$cell <- factor(tb$Celltype,levels=tb$Celltype)
p.pie <- ggplot(tb, aes(x = "", y = Freq, fill = cell)) + 
  geom_bar(stat = "identity",width = 1) + 
  scale_fill_discrete(breaks = tb$cell, labels = labs) +
  labs(x = "", y = "", title = "") +
  coord_polar(theta = "y")  +
  theme(axis.ticks = element_blank()) 
p.pie
ggsave(p.pie,filename="Total_pie.pdf",width = 14, height = 7,dpi = 600)

#柱状图
head(tb)
tb$percentage <- tb$Freq / sum(tb$Freq)
tb$Celltype <- factor(tb$Celltype,tb$Celltype)
p <- ggplot(tb, aes(x=Celltype, y=percentage, fill=Celltype)) +
  geom_bar(stat="identity") +theme_bw()  +  theme(legend.position ="none") +
  geom_text( aes(x=Celltype, y=percentage, label=labs ), #check_overlap = T,
             size=4, vjust = -1,fontface="bold")
p 
ggsave(p,filename="Total_pie-barplot.pdf",width = 14, height = 7,dpi = 600)


#绘制玫瑰花图：
p2 <- p+coord_polar()+labs(x = "", y = "", title = "Cell Types") + 
  theme(axis.text.y = element_blank()) +     
  theme(axis.ticks = element_blank()) +     
  theme(panel.border = element_blank()) + 
  theme(plot.title = element_text(hjust=0.5,size=14,face = "bold") )+
  theme(axis.text.x=  element_blank() )
p2
ggsave(plot=p2,filename="Total_pie2.pdf",width = 14, height = 7, dpi = 600)

#可视化3：每个sample/group的细胞类型组成柱状图
colnames(phe)
table(phe$orig.ident)
#test=phe[,c("orig.ident" , "patient","group","celltype")]
test <- phe[,c("orig.ident" ,"celltype")]
test$cell <- factor(test$celltype,levels=rev(tb$Celltype))
head(test)
p.bar.sample=ggplot(test, aes(x=orig.ident, fill=cell)) +
  geom_bar(position="fill") + 
  labs(title='Cell proportion of each sample') +theme_bw()  
p.bar.sample
ggsave(p.bar.sample,filename="Bar.sample.pdf",width = 14, height = 7)
### 
