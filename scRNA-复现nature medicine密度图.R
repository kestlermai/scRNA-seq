library(ggplot2)
library(Seurat)
# BiocManager::install('TENxPBMCData')
library(TENxPBMCData)#单细胞数据集
args(TENxPBMCData) 
#使用最经典的pbmc3k数据集
pbmc <- TENxPBMCData(dataset = "pbmc3k")
pbmc
class(pbmc)#"SingleCellExperiment"
#居然不是seurat格式，基因名还是ENSEMBL_ID，那还不如老老实实装seuratdata包
#####转换一下格式
#提取原始稀疏矩阵
ct <- as.data.frame(assay(pbmc))
#改一下列名，不然全是V1V2V3，改成1234
colnames(ct) <- 1:ncol(ct)
#将ENSEMBL-ID转换成Entrez
#BiocManager::install('AnnoProbe')
library(AnnoProbe)
ag <- annoGene(rownames(ct),
            ID_type = 'ENSEMBL',
            species = 'human'
            )
head(ag)
#scrna-seq居然测到lncRNA，说实话生物学真不懂，有大佬解释说lncRNA有polyA尾巴，所以10X平台可以测到
#去除重复基因
ag <- ag[!duplicated(ag$SYMBOL),]
ag <- ag[!duplicated(ag$ENSEMBL),]
head(ag)

#匹配一下行名
pos <- match(ag$ENSEMBL,rownames(ct))
#赋值到原始稀疏矩阵中
ct <- ct[pos,]
rownames(ct) <- ag$SYMBOL
#创建seurat对象
sce <- CreateSeuratObject(counts =  ct, project = 'pbmc3k')
sce
#后面就慢慢注释吧，我直接用SeuratData包算了

#devtools::install_github('satijalab/seurat-data')
#老子的github又犯病了
#本地装吧
# install.packages("C:/Users/maihuanzhuo/Desktop/SeuratData-main/SeuratData-main/",
#                  repos = NULL,type = "source")
library(SeuratData)
library(ggplot2)
library(Seurat)
library(dplyr)
library(grid)
#AvailableData()
options(timeout = 10000)
InstallData("pbmc3k")
#嫌慢就手动下
#https://seurat.nygenome.org/src/contrib/pbmc3k.SeuratData_3.1.4.tar.gz
detach("package:SeuratData")
install.packages("C:/Users/maihuanzhuo/Desktop/pbmc3k.SeuratData_3.1.4.tar.gz",
                 repos = getOption(x = "SeuratData.repo.use"),
                 type = "source")
library(pbmc3k.SeuratData)
data("pbmc3k.final")
pbmc3k <- UpdateSeuratObject(object = pbmc3k.final)

# install.packages("C:/Users/maihuanzhuo/Desktop/ggSCvis-master/ggSCvis-master/",
#                  repos = NULL,type = "source")
# #ERROR: dependency 'ggcirclize' is not available for package 'ggSCvis'
# install.packages("C:/Users/maihuanzhuo/Desktop/ggcirclize-master/",
#                  repos = NULL,type = "source")
# #ERROR: dependency 'DescTools' is not available for package 'ggcirclize'
# install.packages("DescTools")
library(ggSCvis)
library(ggh4x)
library(cowplot)
library(viridis)

genes <- c("LYZ", "GNLY", "CD3E")

p1 <- ggscplot(object = pbmc3k, features = genes) +
  geom_density2d(bins = 8, show.legend = F, color = "grey50") +
  geom_scPoint(aes(color = value)) +
  scale_color_gradient(low = "white", high = "#CC0033", name = "gene expression") +
  facet_wrap(~ gene_name, ncol = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = rel(1)),
        axis.text = element_blank())
p1
#填充颜色
p2 <- ggscplot(object = pbmc3k, features = genes) +
  geom_density2d_filled(bins = 10, show.legend = F) +
  geom_scPoint(aes(color = value), size = 0.2) +
  scale_color_gradient(low = "white", high = "#CC0033", name = "gene expression") +
  facet_wrap(~gene_name,ncol = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = rel(1)),
        axis.text = element_blank()) +
  scale_fill_manual(values = colorRampPalette(c("white","#336633"))(10))
p2

p3 <- ggscplot(object = pbmc3k) +
  geom_scPoint(aes(color = seurat_annotations,cluster = seurat_annotations))
p3

p4 <- ggscplot(object = pbmc3k,features = genes) +
  geom_scPoint(color = "grey80") +
  geom_density2d_filled(data = . %>%
                          filter(seurat_annotations %in% c("CD14+ Mono","NK","Naive CD4 T")),
                        bins = 6,show.legend = F,alpha = 0.5) +
  scale_fill_manual(values = colorRampPalette(c("white","#336633"))(6)) +
  geom_scPoint(data = . %>%
                 filter(seurat_annotations %in% c("CD14+ Mono","NK","Naive CD4 T")),
               aes(color = value)) +
  scale_color_gradient(low = "white",high = "#CC0033",name = "gene expression") +
  facet_wrap(~gene_name,ncol = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic",size = rel(1)),
        axis.text = element_blank())
p4

p5 <- ggscplot(object = pbmc3k,features = genes) +
  stat_density2d(geom = "raster",aes(fill = ..density..),
                 contour = F) +
  geom_scPoint(aes(color = value),size = 0.2) +
  scale_color_gradient(low = "white",high = "black",name = "gene expression") +
  facet_wrap(~gene_name,ncol = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic",size = rel(1)),
        axis.text = element_blank()) +
  scale_fill_viridis_c(option = "rocket",direction = -1,name = "cell density") +
  coord_cartesian(expand = F)
p5

p6 <- ggscplot(object = pbmc3k,features = genes) +
  stat_density2d(geom = "raster",aes(fill = ..density..),
                 contour = F) +
  geom_scPoint(aes(color = value),size = 0.2) +
  scale_color_gradient(low = "black",high = "white",name = "gene expression") +
  facet_wrap(~gene_name,ncol = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic",size = rel(1)),
        axis.text = element_blank()) +
  scale_fill_viridis_c(option = "rocket",direction = 1,name = "cell density") +
  coord_cartesian(expand = F)
p6

#按照数据集进行分组
pbmc3k$group <- sample(LETTERS[1:3], 2638, replace = T)

p7 <- ggscplot(object = pbmc3k) +
  stat_density2d(geom = "raster", #绘制栅格填充色样式
                 aes(fill = ..density..),#颜色基于密度值
                 contour = F, #不设置等高线
                 show.legend = F) +#颜色映射与密度相关，而不是图例中的离散变量，所以不显示
  geom_scPoint(color = "white",size = 0.2) +
  facet_wrap(~ group, ncol = 1) +
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.spacing.y = unit(0.1,'cm'),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.line.y = element_blank()) +
  scale_fill_viridis_c(option = "magma",direction = 1) +#viridis，magma，inferno，cividis，plasma
  coord_cartesian(expand = F)+
  xlab("")+
  ylab("")
p7

p8 <- ggscplot(object = pbmc3k) +
  geom_scPoint(aes(color = seurat_annotations,
                   cluster = seurat_annotations),
               show.legend = F,
               label.gp = gpar(fontsize = 8,fontface = "bold.italic")) +
  #facet_wrap与facet_grid的区别是能够自定义行列数
  facet_wrap2(~ group, ncol = 1, 
              strip.position = "left", 
              scales = "free",
              axes = 'x',
              strip = strip_nested(background_y = elem_list_rect(fill = c('#E15759','#85C17E','#7DA4D2')))
              )+
  #theme_bw() +
  theme_classic()+
  theme(panel.grid = element_blank(),
        #panel.background = element_rect(fill='white',color = 'white'),
        panel.spacing.y = unit(0.1,'cm'),
        axis.ticks = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(color = 'white'),
        strip.text.y.left = element_text(face = "bold.italic",
                                         size = rel(1), 
                                         angle = 0,
                                         color = 'white'),
        axis.line.y = element_line(color = "white"),
        axis.text = element_blank()
        ) +
  xlab("")+
  ylab("")
p8
#组合
plot_grid(plotlist = list(p8,NULL,p7), ncol = 3,rel_widths = c(1,0.05,1))
ggsave("C:/Users/maihuanzhuo/Desktop/scRNA-seq密度图.pdf", width = 12, height = 10, dpi = 300)
