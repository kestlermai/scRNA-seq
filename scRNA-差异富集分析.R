library(Seurat)
options(timeout = 100000)
# BiocManager::install("MAST")
library(MAST)
fibro <- readRDS("C:/Users/maihuanzhuo/Desktop/fibro.rds")
fibro
#提取表达矩阵
exp <- GetAssayData(fibro, assay = "RNA")

target_gene = 'CTSK'
#分组
group <- ifelse(exp[c(target_gene),] > median(exp[c(target_gene),]), "CTSK_High", "CTSK_Low")
group <- factor(group,levels = c("CTSK_High","CTSK_Low"))
#
Idents(fibro) <- group
#
table(Idents(fibro))
table(fibro$celltype)
#给FindAllMarkers加线程
library(future)
availableCores()#32 #查看几个核可用
nbrOfWorkers()#1 当前可用的核有多少个
plan()
plan(multisession, workers = 32)#怕不是要跑崩，用scanpy跑岂不是更好
#设置可用的内存
options(future.globals.maxSize = 10 * 1024^3)#10g

#进行亚群比较
degs <- lapply(unique(fibro$celltype), function(x){
  FindMarkers(fibro[,fibro$celltype==x], test.use = "wilcox",#MAST
              ident.1 = 'CTSK_High', ident.2 = 'CTSK_Low',#latent.vars = c("sex"),#MAST法矫正协变量
              logfc.threshold = 0.25, only.pos = F)
})
#
do.call(rbind,lapply(degs, function(x){
  table(x$p_val_adj < 0.05 )
}))
#循环命名
cellnames <- unique(fibro$celltype)
named_degs <- list()
for (i in seq_along(cellnames)){
    current_degs <- degs[[i]]
    filtered_degs <- current_degs[current_degs$p_val_adj < 0.05, ]#过滤P值，P值是Bonferroni校正的
    named_degs[[cellnames[i]]] <- filtered_degs
}
#导出文件
output_folder <- "C:/Users/maihuanzhuo/Desktop/zhanghuilin_output_degs/fibro/"
for (i in names(named_degs)) {
  data_frame <- named_degs[[i]]
  output_file <- file.path(output_folder, paste0(i, "_degs.csv"))
  write.csv(data_frame, file = output_file, row.names = T, col.names = !file.exists(output_file))
}

#整群细胞差异分析
allcluster_degs <- FindMarkers(fibro, test.use = "wilcox",#MAST
                    ident.1 = 'CTSK_High', ident.2 = 'CTSK_Low',
                    logfc.threshold = 0.25, only.pos = F, group.by = NULL)
#
filtered_allcluster_degs <- allcluster_degs[allcluster_degs$p_val_adj < 0.05, ]#过滤P值
#同时过滤P值和log2FC值
filtered_allcluster_degs <- allcluster_degs[allcluster_degs$p_val_adj < 0.05 & abs(allcluster_degs$avg_log2FC) > 1, ]
write.csv(filtered_allcluster_degs, paste(output_folder, "filtered_allcluster_degs.csv",sep = ""))
#
#FindMarkers跑太慢，试试cosg包，参数说是跟FindMarker一样
#remotes::install_github(repo = 'genecell/COSGR')
library(COSG)
degs_cosg <- cosg(fibro, groups = c('CTSK_High','CTSK_Low'),#默认groups='all'
                  assay = 'RNA',slot = 'data')
####
macr <- readRDS("C:/Users/maihuanzhuo/Desktop/macrophage.rds")
macr <- qread("C:/Users/maihuanzhuo/Desktop/zhanghuilin/macrophage.qs")
macr
#提取表达矩阵
exp <- GetAssayData(macr, assay = "RNA")

target_gene = 'CTSK'
#分组
group <- ifelse(exp[c(target_gene),]>median(exp[c(target_gene),]), "CTSK_High", "CTSK_Low")
group <- factor(group,levels = c("CTSK_High","CTSK_Low"))
#
Idents(macr) <- group
#
table(Idents(macr))
table(macr$celltype)
#整群细胞差异分析
allcluster_degs <- FindMarkers(macr, test.use = "wilcox",#MAST
                               ident.1 = 'CTSK_High', ident.2 = 'CTSK_Low',
                               logfc.threshold = 0.25, only.pos = F, group.by = NULL)
#FindMarkers跑太慢，试试cosg包，参数说是跟FindMarker一样
#remotes::install_github(repo = 'genecell/COSGR')
library(COSG)
degs_cosg <- cosg(macr, groups = c('CTSK_High','CTSK_Low'),#默认groups='all'
                  assay = 'RNA',slot = 'data')
#好像跟FindMarkers结果有点出入

filtered_allcluster_degs <- allcluster_degs[allcluster_degs$p_val_adj < 0.05, ]#过滤P值
#同时过滤P值和log2FC值
filtered_allcluster_degs <- allcluster_degs[allcluster_degs$p_val_adj < 0.05 & abs(allcluster_degs$avg_log2FC) > 1, ]
#导出文件
output_folder <- "C:/Users/maihuanzhuo/Desktop/zhanghuilin_output_degs/macrophage/"
write.csv(filtered_allcluster_degs, paste(output_folder, "filtered_allcluster_degs.csv",sep = ""))

#富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
#BiocManager::install("ReactomePA")
library(ReactomePA)
library(ggplot2)
library(enrichplot)
#这里用的是cosg跑出来的差异基因，其实就是一个list
head(degs_cosg[[1]])
head(degs_cosg[[1]]$CTSK_High)
#GO富集
#BP
BP <- compareCluster(degs_cosg[[1]], fun = 'enrichGO', 
                      OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL', ont = "BP")
#MF
MF <- compareCluster(degs_cosg[[1]], fun = 'enrichGO', 
                     OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL', ont = "MF")
#CC
CC <- compareCluster(degs_cosg[[1]], fun = 'enrichGO', 
                     OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL', ont = "CC")
#绘图
p_BP <- dotplot(BP, showCategory = 10)+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text.y = element_text(size = 12), panel.spacing = unit(5, "mm"))+
  scale_colour_gradientn(colours = c('#E15759','#7DA4D2'))
p_BP

#KEGG富集
#KEGG进行注释时，需要将全部的symbol需要转为entrezID
symbols_list <- degs_cosg[[1]]
entrezID_list <- lapply(symbols_list, function(y){
  out <- bitr(y,fromType = 'SYMBOL', 
              toType = 'ENTREZID', OrgDb = org.Hs.eg.db)[, 2]
})
# 找到最长向量长度 
max_len <- max(sapply(entrezID_list, length))
max_len
# 用NA补齐其他向量，因为不一定能够所有的ID都能完成转换，会有缺失
# 将list直接转为data.frame
my_list <- lapply(entrezID_list, function(x) {
  if(length(x) < max_len){
    #c函数将两个向量连接在一起，使两个list连接在一起
    c(x, rep(NA, max_len - length(x)))#用于将缺失值填充到每个向量的末尾，使它们的长度一致。
  }else{
    x
  }
})
my_list
# 转换为数据框
gene <- as.data.frame(my_list)
colnames(gene) <- names(entrezID_list)
gene
#KEGG富集
KEGG <- compareCluster(gene, fun = "enrichKEGG",
                     organism = 'hsa', pvalueCutoff = 0.05)
#绘图
KEGG_plot <- dotplot(KEGG, showCategory = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text.y = element_text(size = 12), panel.spacing = unit(5, "mm"))+
  scale_colour_gradientn(colours = c('#E15759','#7DA4D2'))
KEGG_plot

####ReactomePA富集，一样也是用entrezID
ReactomePA <- compareCluster(gene, fun = "enrichPathway",
                             organism = 'human', pvalueCutoff = 0.05)
#绘图
ReactomePA_plot <- dotplot(ReactomePA, showCategory = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text.y = element_text(size = 12), panel.spacing = unit(5, "mm"))+
  scale_colour_gradientn(colours = c('#E15759','#7DA4D2'))
ReactomePA_plot

#如果是用FindMarkers跑出来的结果，其实就是bulk rna分析那一套
#GO富集分析
ego_All <- enrichGO(gene = row.names(filtered_allcluster_degs),
                    OrgDb = 'org.Hs.eg.db',
                    keyType = 'SYMBOL',
                    ont = "ALL", #设置为ALL时BP, CC, MF都计算
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05)
#绘图
GO_all_plot <- barplot(ego_All,
                       showCategory = 10,
                       split = "ONTOLOGY",
                       label_format = 30,# 默认对名字超过30个字符的进行折叠
                       x = "Count",# X轴展示的变量，默认Count,也可以是GeneRatio
                       )+
  facet_grid(ONTOLOGY~., scale="free")
GO_all_plot

#KEGG富集分析
#同样也需要进行转换成ENTREZID
genelist_ENTREZID <- bitr(row.names(filtered_allcluster_degs),
                          fromType = "SYMBOL",toType = "ENTREZID", 
                          OrgDb = 'org.Hs.eg.db')
#pull函数将提取数据框的指定列
genelist_ENTREZID <- pull(genelist_ENTREZID, ENTREZID)
ekegg <- enrichKEGG(gene = genelist_ENTREZID, 
                    organism = 'hsa',
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    pAdjustMethod = "BH")
#绘图
KEGG_plot <- barplot(ekegg, showCategory = 10,
                     x = "Count",# X轴展示的变量，默认Count,也可以是GeneRatio
                     font.size = 12,
                     title = "KEGG enrichment barplot")
KEGG_plot

#说实话，GO的图有点丑，美化一下吧
library(ggplot2)
library(ggprism)
library(stringr)
library(ggh4x)
GO_plot_data <- ego_All@result
#首先我们要提取BP、CC、MF前十的通路的数据
ONTOLOGY <- c("BP", "CC", "MF")
top10_GO_data <- lapply(ONTOLOGY, function(x) {
  GO_plot_data[GO_plot_data$ONTOLOGY == x, ][1:10, ]
})
#合并数据
top10_GO_data <- do.call(rbind, top10_GO_data)
#进行因子排序
top10_GO_data$Description <- factor(top10_GO_data$Description, 
                                    levels = rev(top10_GO_data$Description))
colnames(top10_GO_data)
#文字颜色映射
colors <- c(rep('#279D77', 10), rep('#CF6611', 10), rep('#7974A1', 10))
#这里GeneRatio是字符型，无法转换成数值型，我们重新计算一下并*100
ratios <- strsplit(top10_GO_data$GeneRatio, "/")
numerators <- sapply(ratios, function(x) as.numeric(x[1]))
denominators <- sapply(ratios, function(x) as.numeric(x[2]))
top10_GO_data$GeneRatio_1 <- numerators / denominators * 100

#条形图
p <- ggplot(top10_GO_data, aes(x = Description, y = Count, fill = ONTOLOGY))+
  geom_bar(stat = "identity",width = 0.8)+
  scale_fill_manual(values = c("#279D77", "#CF6611", "#7974A1"), name = "ONTOLOGY") +
  theme_prism(border = T)+
  theme(axis.text.y = element_text(color = rev(colors)),
        legend.position = "right",
        legend.title = element_text(color = 'black', face = 'bold', size = 12),
        legend.text = element_text(color = 'black', face = 'bold', size = 12))+
  coord_flip()
p

#棒棒图
p <- ggplot(top10_GO_data, aes(x = Description, y = Count))+
  geom_segment(aes(y = 0, yend = Count, x = Description, xend = Description),
               color = 'grey', linewidth = 1) +#添加横线，即棒棒糖图
  geom_point(aes(color = -log10(p.adjust), size = GeneRatio_1)) +
  scale_color_gradient(low = 'blue', high = 'red',  breaks = seq(10, 20, 5), limits = c(5, 20),
                       name = 'Significance\n(-log10 FDR)') +
  scale_size_continuous(name = 'GeneRatio(%)', range = c(2, 6), limits = c(0.5, 4.5)) +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+ #设置间隔距离
  # scale_x_discrete(labels = function(x) str_wrap(x, width = 50))+ #调整label长度
  #facet_grid(ONTOLOGY ~ ., scales = 'free')+
  # facet_wrap2(.~ ONTOLOGY, ncol = 1, strip.position = "right", 
  #             scales = "free_y", 
#             strip = strip_nested(background_y = elem_list_rect(fill = c('#279D77','#CF6611','#7974A1')))
# )+
  theme_prism(border = T)+
  theme(axis.text.y = element_text(color = rev(colors)),
        legend.position = "right",
        legend.title = element_text(color = 'black', face = 'bold', size = 12),
        legend.text = element_text(color = 'black', face = 'bold', size = 12))+
  coord_flip()
p
#添加颜色背景
rect.data <- data.frame(
  xmin = c(0, 10.5, 20.5),
  xmax = c(10.5, 20.5, 30.6),
  colors = letters[1:3]
)
head(rect.data)
p1 <- p + 
  geom_rect(data = rect.data, aes(ymin = -Inf, ymax = Inf, xmin = xmin, xmax = xmax, fill = colors), 
            inherit.aes = FALSE, alpha = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = c('#7974A1', '#CF6611', '#279D77'))+#颜色是从下往上映射的
  guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))#调整图例顺序
p1

###GSEA
#需要根据logFC对基因进行排序
#ENTREZID转换
filtered_allcluster_degs$symbol <- rownames(filtered_allcluster_degs)
df <- bitr(unique(filtered_allcluster_degs$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
#
DEG <- filtered_allcluster_degs
DEG <- merge(DEG,df,by.y='SYMBOL',by.x='symbol')
colnames(DEG)
#
data_all_sort <- DEG %>% 
  arrange(desc(avg_log2FC))
#把foldchange按照从大到小提取出来
geneList <- data_all_sort$avg_log2FC
#给上面提取的foldchange加上对应上ENTREZID
names(geneList) <- data_all_sort$ENTREZID
head(geneList)
#GSEA-KEGG
kegg_ges <- gseKEGG(geneList = geneList,
                    organism = 'hsa',
                    nPerm = 1000,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    eps = 0)
#
kegg_result <- as.data.frame(kegg_ges)
rownames(kegg_ges@result)[head(order(kegg_ges@result$enrichmentScore))]
af <- as.data.frame(kegg_ges@result)
#write.table(af,file=paste0("2.","all_GSEA.xls"),sep="\t",quote=F,col.names=T)

#排序后分别取GSEA结果enrichmentScore负向关前10个和正相关前10个
pdf(paste0("top10_down_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kegg_ges, geneSetID = rownames(kegg_ges@result)[head(order(kegg_ges@result$enrichmentScore),10)])
dev.off()
pdf(paste0("top10_up_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kegg_ges, geneSetID = rownames(kegg_ges@result)[tail(order(kegg_ges@result$enrichmentScore),10)])
dev.off()

#GO BP enrich
BP_ges <- gseGO(geneList = geneList,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                eps = 0)
#
BP_ges_result <- BP_ges@result

#试试hallmark geneset
library(msigdbr)
H.gmt <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
#这里用的是gene symbol不是ENTREZID
gene_list <- data_all_sort$avg_log2FC
names(gene_list) <- data_all_sort$symbol
head(gene_list)
#
Hallmark_gse <- GSEA(gene_list, TERM2GENE = H.gmt)
Hallmark_gse_result <- Hallmark_gse@result
#
gseaplot2(Hallmark_gse,
          geneSetID = rownames(Hallmark_gse@result)[tail(order(Hallmark_gse@result$NES),5)],
          base_size = 14)



#参考NC富集条形图美化，这里用的是自己分析单细胞数据HIV-control的PBMC
library(qs)
library(Seurat)
library(MAST)
#读取数据
#sce.all <- readRDS('C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8/sce.all_marker_cluster.rds')
sce.all <- qread("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8/sce.all_marker_cluster.qs")
#qsave(sce.all, "C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8/sce.all_marker_cluster.qs")
sce.all
table(sce.all@active.ident)
table(sce.all$group)
##提取HIV患者的细胞出来
# sce.hiv <- subset(sce.all, cells = rownames(sce.all@meta.data[sce.all@meta.data$group == 'HIV',]))
# dim(sce.hiv)
#赋值一下分组信息
Idents(sce.all) <- sce.all$group
Idents(sce.all)
#然后就是每个cluster进行差异分析了
degs <- lapply(unique(sce.all$celltype), function(x){
  FindMarkers(sce.all[,sce.all$celltype == x], test.use = "wilcox",#MAST,wilcox
              ident.1 = 'control', ident.2 = 'HIV',#latent.vars = c("sex"),#MAST法矫正协变量
              logfc.threshold = 0.25, only.pos = F)
})

#FindMarker跑的MAST巨慢，可以考虑用cosg包
# library(COSG)
# degs_cosg <- cosg(sce.all, groups = c('control','HIV'),#默认groups='all'
#                   assay = 'RNA', slot = 'data')

#看看总体筛选标准情况
do.call(rbind,lapply(degs, function(x){
  table(x$p_val_adj < 0.01 & abs(x$avg_log2FC) > 2)
}))
# FALSE TRUE
# [1,]  7752 1556
# [2,]  7562 1994
# [3,]  8850 1177
# [4,]  7040 2689
# [5,] 10752  228
# [6,]  9902  926
#还是筛出很多个基因
#保存一下degs结果
qsave(degs,'C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8/sce.all_degs.qs')
# degs <- qread('C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8/sce.all_degs.qs')
#循环命名
cellnames <- unique(sce.all$celltype)
cellnames <- gsub(" ", "_", cellnames)
named_degs <- list()
for (i in seq_along(cellnames)) {  
  current_degs <- degs[[i]]  
  filtered_degs <- current_degs[current_degs$p_val_adj < 0.01 & abs(current_degs$avg_log2FC) > 2, ]  
  named_degs[[cellnames[i]]] <- filtered_degs  
}

#####富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggprism)

BP_results <- list()
#GO富集分析
for (i in names(named_degs)){
  go_BP <- enrichGO(gene = row.names(named_degs[[i]]), OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL',
                    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
  BP_results[[i]] <- go_BP@result[1:3,]
}
#生成数据框
BP_merged <- data.frame(  
  Name = rep(names(BP_results), sapply(BP_results, nrow)),  
  do.call(rbind, BP_results)  
)
BP_merged <- BP_merged[,c(1,3,7)]
colnames(BP_merged)
#对每个细胞行下方插入blank行
blank <- data.frame(Name = paste("blank", 1:5, sep = ""), 
                    Description = paste("blank", 1:5, sep = ""), 
                    p.adjust = rep("9.079409e-08", 5))
blank$p.adjust <- as.numeric(blank$p.adjust)
#合并数据
df <- rbind(BP_merged, blank)
# #定义数据顺序
df$Name <- factor(df$Name, levels = c("T_cells", "blank1", "NK_cells", "blank2",
                                      "B_cells", "blank3", "Monocytes", "blank4",
                                      "Megakaryocytes", "blank5","Plasma"))
# #根据df$Name定义的顺序进行排列
df_sorted <- df[order(df$Name), ]
#发现出现重复值，正常因子顺序无法进行定义
# df$Description <- factor(df$Description, levels = rev(df$Description))
#还是在excel里面改吧，在R太麻烦了
write.csv(df_sorted,'C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8/df_sorted.csv')
df_sorted <- read.csv('C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8/df_sorted.csv',row.names = 1)
#重新定义一下
cell_levels <- c("T_cells", "blank1", "NK_cells", "blank2","B_cells", "blank3", 
                "Monocytes", "blank4", "Megakaryocytes", "blank5", "Plasma")
df_sorted$Name <- factor(df_sorted$Name, levels = cell_levels)
df_sorted <- df_sorted[order(df_sorted$Name), ]
df_sorted$Description <- factor(df_sorted$Description, levels = rev(df_sorted$Description))
#配色工具cols4all
c4a_gui()
c4a('superfishel_stone',7)
colors <- c("#6388B4","white","#FFAE34","white","#EF6F6A","white",
            "#55AD89","white","#C3BC3F","white","#BB7693")
#自定义函数修改标签颜色：
col_function <- function(x){
  col <- rep("black", length(x))
  blank <- which(x %in% paste0(rep('blank', times = 6), 1:6))
  col[blank] <- 'white' #将blank开头的文本标签设定为白色，实现"隐藏"
  col
}
ycol <- col_function(df_sorted$Description)
#根据目标位置创建一个新的图例数据集：
legend <- data.frame(x = rep(10, 6),
                     y = seq(2, 22, 4),
                     label = rev(c("T_cells",  "NK_cells", "B_cells",  
                               "Monocytes",  "Megakaryocytes", "Plasma")))
legend
#绘图
p <- df_sorted %>% ggplot() +
  geom_bar(aes(x = -log10(p.adjust), y = Description, fill = Name),
           stat = "identity")+
  geom_point(data = legend, aes(x = x, y = y),
             col = rev(c("#6388B4","#FFAE34","#EF6F6A","#55AD89","#C3BC3F","#BB7693")),
             shape = 15, size = 6)+
  geom_text(data = legend,aes(x = x + 0.3, y = y, label = label), 
            hjust = 0, label.size = 5,fontface = "bold") +
  scale_x_continuous(limits = c(0, 12.5), breaks = seq(0, 9, by = 3))+ #坐标轴自定义(为自定义图例留出足够空隙)
  scale_fill_manual(values = colors)+
  labs(title = 'Top ontologies (biological process)')+
  theme_prism(border = T)+
  theme(legend.position = 'none',
        axis.text.y = element_text(color = ycol),
        axis.ticks.y = element_blank())
p
ggsave('C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8/富集分析条形图美化.pdf',
       p, width = 12, height = 8, dpi = 300)

########兴致勃勃来也画个Circular barplot看看
# 首先明确一下思路，这个barplot类似一个堆积图，代表不同cluster所在的比例
# 然后就是中间内圈的分组，这个分组你可以是类似通路（代谢、细胞凋亡）分成一组ABCD
# 这里我就进行GO、KEGG、ReactomePA进行分组了
library(qs)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(tidyverse)
library(ggplot2)
library(Seurat)
#
sce.all <- qread("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8/sce.all_marker_cluster.qs")
#赋值一下分组信息
Idents(sce.all) <- sce.all$group
Idents(sce.all)
#
degs <- qread('C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8/sce.all_degs.qs')
#循环命名
cellnames <- unique(sce.all$celltype)
cellnames <- gsub(" ", "_", cellnames)
named_degs <- list()
for (i in seq_along(cellnames)) {
  current_degs <- degs[[i]]
  filtered_degs <- current_degs[current_degs$p_val_adj < 0.01 & abs(current_degs$avg_log2FC) > 2, ]
  named_degs[[cellnames[i]]] <- filtered_degs
}
# qsave(named_degs, 'C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8/sce.all_degs_filted.qs')
# named_degs <- qread('C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8/sce.all_degs_filted.qs')
#GO富集分析
BP_results <- list()
MF_results <- list()
CC_results <- list()
for (i in names(named_degs)){
  go_BP <- enrichGO(gene = row.names(named_degs[[i]]), OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL',
                    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
  go_MF <- enrichGO(gene = row.names(named_degs[[i]]), OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL',
                    ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
  go_CC <- enrichGO(gene = row.names(named_degs[[i]]), OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL',
                    ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
  BP_results[[i]] <- go_BP@result
  MF_results[[i]] <- go_MF@result
  CC_results[[i]] <- go_CC@result
}
# KEGG富集分析之前需要转换一下，把SYMBOL改为ENTREZID
kegg_result <- list()
for (i in names(named_degs)){
  genelist_ENTREZID <- bitr(row.names(named_degs[[i]]), fromType = "SYMBOL",
                            toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
  genelist_ENTREZID <- pull(genelist_ENTREZID, ENTREZID)
  ekegg <- enrichKEGG(gene = genelist_ENTREZID, organism = 'hsa', 
                      pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
  kegg_result[[i]] <- ekegg@result
}
#####其实上面的函数可以换成多组富集分析compareCluster
# 把named_degs的rownames全部提取出来
degs_genelist <- lapply(named_degs, function(x) rownames(x))
#转ID 
for (i in names(named_degs)){
  degs_genelist[[i]] <- as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                                  keys = degs_genelist[[i]],
                                                                  columns = 'ENTREZID',
                                                                  keytype = 'SYMBOL')[,2]))
}
lapply(degs_genelist, head)
#GO富集分析
compareCluster.go <- compareCluster(degs_genelist, fun = "enrichGO", ont = 'BP',
                                    OrgDb = 'org.Hs.eg.db', pvalueCutoff = 0.05) 
table(compareCluster.go@compareClusterResult$Cluster) #每个基因集富集个数
df_go_diff <- as.data.frame(compareCluster.go)
## 选择富集到两组样本以上的通路可视化
df <- as.data.frame(table(df_go_diff$Description))
#提取进行降序
df <- df[order(df$Freq,decreasing = T),][1:10,]
#根据进行过滤
df_go_diff <- df_go_diff[df_go_diff$Description %in% df$Var1,]
#重新排序
df_go_diff <- df_go_diff[order(df_go_diff$Description,decreasing = F),]

#KEGG富集分析
compareCluster.kegg <- compareCluster(degs_genelist, fun = "enrichKEGG",
                                      organism = 'hsa', pvalueCutoff = 0.05) 
table(compareCluster.kegg@compareClusterResult$Cluster) #每个基因集富集个数
df_kegg_diff <- as.data.frame(compareCluster.kegg)
## 选择富集到两组样本以上的通路可视化
df <- as.data.frame(table(df_kegg_diff$Description))
#提取进行降序
df <- df[order(df$Freq,decreasing = T),][1:7,]
#根据进行过滤
df_kegg_diff <- df_kegg_diff[df_kegg_diff$Description %in% df$Var1,]
#重新排序
df_kegg_diff <- df_kegg_diff[order(df_kegg_diff$Description,decreasing = F),]

#pathway富集分析
compareCluster.pathway <- compareCluster(degs_genelist, fun = "enrichPathway",
                                         organism = 'human', pvalueCutoff = 0.05) 
table(compareCluster.pathway@compareClusterResult$Cluster) #每个基因集富集个数
df_pathway_diff <- as.data.frame(compareCluster.pathway)
## 选择富集到两组样本以上的通路可视化
df <- as.data.frame(table(df_pathway_diff$Description))
#提取进行降序
df <- df[order(df$Freq,decreasing = T),][1:8,]
#根据进行过滤
df_pathway_diff <- df_pathway_diff[df_pathway_diff$Description %in% df$Var1,]
#重新排序
df_pathway_diff <- df_pathway_diff[order(df_pathway_diff$Description,decreasing = F),]

######目前以上富集分析了GO、KEGG、pathway（ReactomePA）
merged_df <- bind_rows(df_go_diff, df_kegg_diff, df_pathway_diff)
## 分割弦图需添加NA Set a number of 'empty bar'
to_add <- matrix(NA, 3, ncol(merged_df))
colnames(to_add) <- colnames(merged_df)
merged_df <- rbind(merged_df, to_add)
merged_df$group <- factor(c(rep('GO',50),
                            rep('KEGG',29),
                            rep('ReactomePA',38),
                            'GO','KEGG','ReactomePA'
))
#
merged_df_sorted <- merged_df[order(merged_df$group), ]
# #给相同通路添加序号
# merged_df_sorted <- merged_df_sorted %>%
#   mutate(id = group_indices(., Description))
#没有按照我想要的顺序，还是手动添加吧，excel大法也行
#统计一下每个通路的个数
tmp.data <- data.frame(table(merged_df_sorted$Description))
#
merged_df_sorted$id <- factor(c(rep(1:10, each = 5), #BP
                                '11',#NA
                                rep(12:15, each = 4),#KEGG
                                rep(16, each = 5),#KEGG
                                rep(17:18, each = 4),#KEGG
                                '19',#NA
                                rep(20:21, each = 5),#pathway
                                rep(22:23, each = 4),#pathway
                                rep(24:27, each = 5),#pathway
                                '28'#NA
))
merged_df_sorted$id <- factor(merged_df_sorted$id, levels = c(1:28))
#######
# 确定标签及倾斜角度 Get the name and the y position of each label
label_data <- merged_df_sorted %>% 
  group_by(id, Description) %>% 
  summarize(tot = sum(-log10(p.adjust)))
#
number_of_bar <- nrow(label_data)
#转换一下label_data$id
label_data$id <- as.numeric(label_data$id)
# I substract 0.5 because the letter must have the angle of the center of the bars. 
#Not extreme right(1) or extreme left (0)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
label_data
#####
## 添加分组 prepare a data frame for base lines
merged_df_sorted$id <- as.numeric(merged_df_sorted$id)
base_data <- merged_df_sorted %>% 
  group_by(group) %>% 
  summarize(start = min(id), end = max(id) - 1) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
base_data
#
colors <- c("#6388B4","#FFAE34","#EF6F6A","#55AD89","#C3BC3F","#BB7693")
#
p <- ggplot(merged_df_sorted,aes(x = as.factor(id), y = -log10(p.adjust), fill = Cluster))+
  geom_bar(stat = "identity", alpha = 0.9)+
  geom_text(data = label_data, aes(x = id, y = tot + 0.3, label = Description, hjust = hjust), 
            color = "black", fontface = "bold", size = 3, angle = label_data$angle, inherit.aes = F)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50))+ #调整label长度
  geom_segment(data = base_data, aes(x = start, y = -0.5, xend = end, yend = -0.5), 
               colour = "black", linewidth = 0.8, inherit.aes = F)+
  geom_text(data = base_data, aes(x = title, y = -1, label = group), #hjust = c(1,1,0,0),#会重叠,AI调整吧
            colour = "black", size = 3, fontface = "bold", inherit.aes = F)+
  annotate("text", x = c(rep(max(merged_df_sorted$id),5),24.5), #定义标签水平位置
           y = c(0,20,40,60,80,100), 
           label = c("0", "5", "10","15","20","-log10(p.adjust)") , 
           color = "grey", size = 3 , 
           angle = c(rep(0,5),0), #定义文本标签的旋转角度
           fontface = "bold", hjust = 1)+
  scale_fill_manual(values = colors)+
  ylim(-30,100)+
  theme_minimal() +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        panel.grid = element_blank())+
  coord_polar(start = 0)
p
ggsave("C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-HIV/2-marker-8.8/Circular barplot.pdf", 
       p, width = 12, height = 12, dpi = 300)
