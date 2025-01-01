plot.clusters.group <- function(data = seurat_data,clusters = seurat_clusters,
                                group = orig.ident,widths = c(3,1),log =TRUE,
                                legend.title = "Group",color = 1,xlab = ""){
  ## take an integrated Seurat object, plot distributions over orig.ident
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  library(paletteer)
  library(ggsci)
  mytheme = theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                  axis.title = element_text(size = 10,color ="black"),
                  axis.text = element_text(size=10,color = "black"),
                  #axis.line = element_line(color = "black"),
                  #axis.ticks = element_line(color = "black"),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  # panel.grid=element_blank(), # 去网格线
                  # legend.position = "none",
                  legend.text = element_text(size=8),
                  legend.title= element_text(size= 8),
                  # axis.text.x = element_text(angle = 45, hjust=1, vjust=1)
  )
  
  count_table <- table(data@meta.data[,clusters], data@meta.data[,group])
  count_mtx <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)
  
  cluster_size <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
  
  if("0" %in% cluster_size$cluster){
    sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
  }else{
    sorted_labels <- paste(cluster_size$cluster[order(cluster_size$value)])
  }
  
  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
  
  colnames(melt_mtx)[2] <- "dataset"
  
  ################### p1  
  if(log){
    p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",
                                                                     fill = "#F28E2B") +#"grey60"
      theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("") + mytheme
  }else{
    p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "#F28E2B") +
      theme_bw() + xlab("Cells per cluster") + ylab("") + mytheme
  }
  ################### p2    
  ########################### color 1  
  if(color==1){
    if(length(unique(melt_mtx$dataset)) < 21){
      p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
        geom_bar(position="fill", stat="identity",) + theme_bw() + coord_flip() +
        #scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.9))+
        #scale_fill_aaas(palette = c("default"), alpha = 0.8)+
        scale_fill_manual(values = c("#4E79A7","#E15759"))+
        ylab(paste0("Fraction of cells in each ",tolower(legend.title))) + xlab(xlab) +
        theme(legend.position="top") + guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme
    }else{
      warning("The color limit is <21")
      p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
        geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() +
        ylab(paste0("Fraction of cells in each ",tolower(legend.title))) + xlab(xlab) +
        theme(legend.position="top") + guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme
    }
  }
  ########################### color 2
  if(color==2){
    if(length(unique(melt_mtx$dataset)) < 9){
      p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
        geom_bar(position="fill", stat="identity",) + theme_bw() + coord_flip() +
        scale_fill_brewer(palette = "Set2")+
        ylab(paste0("Fraction of cells in each ",tolower(legend.title))) + xlab(xlab) +
        theme(legend.position="top") + guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme
    }else{
      warning("The color limit is <9")
      p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
        geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() +
        ylab(paste0("Fraction of cells in each ",tolower(legend.title))) + xlab(xlab) +
        theme(legend.position="top") + guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme
    }
  }
  wrap_plots(ncol = 2,p2,p1,widths = widths)
}

p <- plot.clusters.group(data = sce.all, clusters =  "celltype", 
                         xlab = "Cluster number", log = TRUE, group = "group",
                         legend.title = "",widths = c(3,2),color = 1)
p
ggsave("group_celltype_percent.pdf", p, units = "cm", width = 28,height = 15, dpi = 600)
