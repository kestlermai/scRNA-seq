# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 23:52:38 2024

@author: maihuanzhuo
"""
# Leiden clustering package
# conda install -c conda-forge leidenalg
# pip install bbknn 这里annoy包安装不上，原因是我的VS C++版本太久，升级完就好了https://my.visualstudio.com/Downloads/Featured?mkt=en-us
# pip install cellhint (cell整合算法)
# pip install celltypist （细胞注释）
# pip install scrublet 
# scanpy用命名函数主要包括pp（preprocessing）预处理（主要是QC）
# tl（tools）工具（包括聚类、降维、轨迹推断）
# pl（plot）绘图

# •adata.X 存储 count matrix，数据类型为稀疏矩阵 scipy.sparse.csr.csr_matrix
# •adata.obs 存储关于 obervations(cells) 的 metadata，数据类型为 dataframe
# •adata.var 存储关于 variables(genes) 的 metadata，数据类型为 dataframe
# •AnnData.uns 存储后续附加的其他非结构化信息
# •adata.obs_names 和 adata.var_names index

# import numpy as np
# import pandas as pd
# import scanpy as sc
# import anndata
# import os
# os.chdir('C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-CRC/2-cluster-7.28') ##修改路径
# # 
# adata = anndata.read_h5ad('sce.all.int.h5ad')
# adata

# # 找marker gene
# sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon') # louvain  leiden
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
# # 提取每个cluster的差异基因，这里leiden算法在0.5分辨率下算出32个cluster
# result = adata.uns["rank_genes_groups"]
# groups = result["names"].dtype.names
# # 展示每个cluster前十的基因
# pd.DataFrame(
#     {
#         group + "_" + key[:1]: result[key][group]
#         for group in groups
#         for key in ["names", "logfoldchanges", "pvals"]
#     }
# ).head(10)
# # 提取top10 marker gene
# marker_genes_top10 = pd.DataFrame(
#     {   
#     f"{group}_{key[:1]}": result[key][group][:10] # 只提取前10个基因
#     for group in groups
#     for key in ["names", "logfoldchanges", "pvals"]
#     }
# )
# #
# marker_genes = pd.DataFrame(
#     {   
#     f"{group}_{key[:1]}": result[key][group]
#     for group in groups
#     for key in ["names", "logfoldchanges", "pvals"]
#     }
# )
# marker_genes.to_csv("marker_genes.csv", index=False)# 注意某些基因会被excel改成日期

# 重新跑吧
import numpy as np
import pandas as pd
import scanpy as sc
import leidenalg
import anndata
import re # 处理正则表达式的库
import matplotlib.pyplot as plt
import os
os.chdir('C:/Users/maihuanzhuo/Desktop/scRNA/scRNA-CRC') ##修改路径
os.makedirs("scanpy")
# 
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
# 输出版本号
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")# set_figure_params 设置图片的分辨率/大小以及其他样式
# facecolor="white"即背景颜色为白色
plt.rcParams['figure.dpi'] = 80 # 将全局图形分辨率设置为300 dpi
# 定义保存图函数
def pic(pdf):
    searchObj = re.search( r'(.*).pdf', pdf)
    png = f"{searchObj.group(1)}.png"
    plt.savefig(pdf, bbox_inches="tight")
    plt.savefig(png, bbox_inches="tight", dpi=100)
    plt.close()

# 获取所有目录（排除根目录）
samples = [d for d in next(os.walk("./raw/raw_data"))[1]]
print(samples)

# # 读取矩阵
# adata = sc.read_10x_mtx(
#     "data/filtered_gene_bc_matrices/hg19/",  # the directory with the `.mtx` file
#     var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
#     cache=True,  # write a cache file for faster subsequent reading
# )
# 索引去重
# adata.var_names_make_unique()
# adata.obs_names_make_unique()

# 这里踩坑了，由于Cell Ranger V2是输出genes而不是features
# 所以使用scanpy的read_10x_mtx()函数进行读取的话，会读取失败
# 读取整合在一起的GSE178318，GSE178318是把所有数据都合在一起了，但是metadata是没有给出具体信息，考虑从cellranger重新下载
from scipy.io import mmread
path = './raw/GSE178318/'
genes = None
cell = None
mtx = None
for name in os.listdir(path):
    file= os.path.join(path, name)
    if 'barcodes.tsv' in name:
        cell = pd.read_csv(file, header=None)
    elif 'genes.tsv' in name:
        genes = pd.read_csv(file, sep='\t', header=None)
    elif 'matrix.mtx' in name:
        mtx = mmread(file)
#创建scanpy对象
adata = anndata.AnnData(mtx)
adata = adata.T
adata.obs.index = pd.Index(cell[0])
adata.var.index = pd.Index(genes[1])
adata.var['gene_ids'] = genes[0].to_list()
adata
print(adata.obs.columns)
# 批量读取数据并创建 Anndata 对象
sc_list = []
for sample in samples:
    folder = os.path.join('./raw/raw_data', sample)
    print(sample)
    print(folder)
    print(os.listdir(folder))
    # 读取数据，10x的三个标准文件 mtx, tsv, tsv
    adata = sc.read_10x_mtx(folder, var_names='gene_symbols', 
                            cache=True)# 在第一次读取后缓存读取结果，加速后续读取
    # 去重
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obs['Dataset'] = sample # 设置项目名称
    sc_list.append(adata)
# 合并数据
adata = sc.concat(sc_list, join='outer', # 外连接
                  label='Dataset', keys=samples)
# 查看每个数据集的细胞数
adata.obs['Dataset'].value_counts()

# QC
os.chdir('./scanpy') ##修改路径
os.makedirs("1-QC-7.31")
os.chdir('./1-QC-7.31')

# 展示前二十高表达基因
sc.pl.highest_expr_genes(adata, n_top=20, save='_raw.png')# 自动命名好了，还自动创建figures文件夹
# plt.savefig('highest_expr_genes.pdf', format='pdf', bbox_inches='tight', dpi=300)# bbox_inches='tight' 确保图像边界紧凑

# QC
# 保留至少在三个细胞中表达的基因，保留至少包含 200 个基因的细胞
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata
# 124700 × 30751

# 计算每个细胞中线粒体、核糖体、血红蛋白基因
# mitochondrial genes
adata.var["mito"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
# 
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mito", "ribo", "hb"], percent_top=None, log1p=True,
    inplace=True #直接在原始数据上去置换，节省内存inplace=True
)
adata
# 小提琴图可视化
sc.pl.violin(adata, 
             ['n_genes_by_counts', 'total_counts', 'pct_counts_mito'],
             jitter=0.4, multi_panel=True, save='_QC_raw.png')
# 可视化特征之间的关系
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mito', color="pct_counts_mito", show=False)
plt.savefig('./figures/scatter_mito.png')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color="n_genes_by_counts", show=False)
plt.savefig('./figures/scatter_feature.png')
# # 组合
# fig, axs = plt.subplots(1, 2, figsize=(12, 5))# 创建子图布局
# # 在子图上绘制散点图
# sc.pl.scatter(adata, x='total_counts', y='pct_counts_mito', ax=axs[0], show=False)
# sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=axs[1], show=False)
# plt.tight_layout()# 调整布局并显示图形
# plt.savefig('./figures/scatter_raw.png')

# 过滤
adata = adata[adata.obs.n_genes_by_counts < 5000, :] # 一般来说是2500（seurat中pbmc教程）
# 121459 × 30751
adata = adata[adata.obs.pct_counts_mito < 25, :]# 一般10%-15%，严格点就5%
adata
# 92719 × 30751
# 过滤后
sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4)
pic('genes_filter.pdf')
sc.pl.violin(adata, ['total_counts'], jitter=0.4)
pic('counts_filter.pdf')
sc.pl.violin(adata, ['pct_counts_mito'], jitter=0.4)
pic('percent_mito_filter.pdf')

# 处理双细胞检测
# scrublet是单个样本量进行处理的
# 拆分数据集并进行Scrublet去除双细胞
import scrublet as scr
adata_list = []
for Dataset in adata.obs['Dataset'].unique():
    adata_single = adata[adata.obs['Dataset'] == Dataset]
    # 计算双细胞评分
    scrub = scr.Scrublet(adata_single.X, expected_doublet_rate=0.05)
    # 预测的双细胞
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, 
                                                              min_gene_variability_pctl=85, n_prin_comps=30)
    # scrub.call_doublets(threshold=0.25) 如果自动阈值检测效果不佳，则可以使用call_doublets()函数调整阈值
    # 绘制双细胞评分直方图
    # scrub.plot_histogram() 
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, n_neighbors=10, 
                                             min_dist=0.3)) # min_dist默认0.1, 值越小会使数据点在低维空间中更紧密地聚集在一起
    scrub.plot_embedding('UMAP', order_points=True)
    save_path = './scrublet' # 设置保存路径
    pic(os.path.join(save_path, f'UMAP_{Dataset}.pdf'))
    # doublets占比
    print(scrub.detected_doublet_rate_)
    adata_single.obs['doublet_scores'] = doublet_scores
    adata_single.obs['predicted_doublets'] = predicted_doublets
    # 创建DataFrame保存双细胞信息
    doublet_info = pd.DataFrame({
        'cell_barcodes': adata_single.obs_names,
        'doublet_score': doublet_scores,
        'predicted_doublet': predicted_doublets
    })
    doublet_info.to_csv(os.path.join(save_path, f'doublet_info_{Dataset}.csv'),index=False)# 保存DataFrame为CSV文件
    # 去除预测为双细胞的细胞
    adata_single = adata_single[adata_single.obs['predicted_doublets'] == False]
    adata_list.append(adata_single)
    
# 检查结果
for adata_single in adata_list:
    print(f"Sample: {adata_single.obs['Dataset'][0]}")
    print(f"Number of cells after doublet removal: {adata_single.shape[0]}")

# 重新合并数据
adata = sc.concat(adata_list, join='outer', label='Dataset', keys=samples)
adata
# 92414 × 30751
adata.obs.Dataset.value_counts()

# 对数正态化基因表达（标准化至每个细胞的10,000计数）存储在 .X 中，而原始Count数据存储在 .raw 中。
# CellHint不依赖于后者，但为了确保单细胞分析的完整性，我们仍然从原始Count数据开始进行操作：
adata.raw = adata
# adata = adata.raw.to_adata()

# 标准化数据
sc.pp.normalize_total(adata, target_sum=1e4)
# 自然对数转换
sc.pp.log1p(adata)
# 识别高变的基因
sc.pp.highly_variable_genes(adata, # n_top_genes=2000, 
                            batch_key = 'Dataset', # 对不同批次的数据分别计算变异性
                            subset = True, # 过滤基因，标记高变异性基因
                            #flavor='seurat' # 使用 Seurat 包的算法来计算变异性
                            )
# 可视化
sc.pl.highly_variable_genes(adata, save='_filt.png')
# 将 AnnData 对象的 .raw 属性设置为经归一化和对数化的原始基因表达值，供之后的可视化分析使用
adata.raw = adata
# # 取出高度差异的基因
# adata = adata[:, adata.var.highly_variable]
# # 校正细胞基因计数和线粒体基因比例的影响
# sc.pp.regress_out(adata, ["total_counts", "pct_counts_mito"])
# 数据缩放
sc.pp.scale(adata, max_value=10)
#
sc.tl.pca(adata, n_comps=30)# 默认是25个PC
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
# 根据数据集进行分组可视化umap结果
sc.pl.umap(adata, color='Dataset', wspace=0.5)
# 
import cellhint
cellhint.integrate(adata, 'Dataset')
# 
sc.tl.umap(adata)
sc.pl.umap(adata, color ='Dataset', wspace = 0.5)
# # PCA降维
# sc.tl.pca(adata, svd_solver='arpack')# svd_solver 指定奇异值分解 SVD 的方法
# # 画图
# sc.pl.pca(adata, color="CST3", save='.png')
# # 碎石图
# sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save='.png')

# cluster
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
# leiden算法
sc.tl.leiden(
    adata,
    resolution=1, # 感觉一开始分辨率1就足够了吧
    random_state=0,
    flavor="igraph",
    n_iterations=2,# Leiden算法的迭代次数
    directed=False,# 设置为无向图，代表细胞之间关系没有方向
)
# louvain算法
sc.tl.louvain(
    adata,
    resolution=0.5,
    random_state=0,
    directed=False,# 设置为无向图，代表细胞之间关系没有方向
)
# UMAP
sc.tl.umap(adata)
# 绘制聚类图
sc.pl.umap(adata, color=['leiden'])
pic('umap_leiden_res_1.pdf')

# 用CellTypist算法自动注释
import celltypist
from celltypist import models
# 下载官方model https://www.celltypist.org/models
models.download_models(force_update = True)
models.models_path # 模型会下载到家目录'~/.celltypist/data/models'
models.models_description() # model description
# 去除重复的索引
adata = adata[~adata.obs.index.duplicated()].copy()
# 自动注释
adata = celltypist.annotate(adata, 
                            model = 'Human_Colorectal_Cancer.pkl', # 这里用结直肠癌患者的结肠组织的细胞类型
                            majority_voting = True).to_adata()# 使用多数投票法来决定细胞类型
# 
adata.obs[['predicted_labels', 'majority_voting', 'conf_score']]
# 牛逼，速度挺快的，比那GPT自动注释快多了。
sc.pl.umap(adata, color = 'majority_voting', save='_celltypist_test_annotation.pdf')
# 找 marker gene
# t 检验
# sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
# sc.settings.verbosity = 2  # reduce the verbosity
# Wilcoxon 秩和检验
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, 
                        n_genes=25, # 选择每个组中表达差异最大的25个基因来展示
                        sharey=False) # sharey 参数用于控制是否所有子图共享同一个y轴
# # 使用逻辑回归对基因进行排名
# sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
# 提取每个cluster的差异基因
result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names
# 展示每个cluster前十的基因
pd.DataFrame(
    {
        group + "_" + key[:1]: result[key][group]
        for group in groups
        for key in ["names", "logfoldchanges", "pvals"]
    }
).head(10)
# 提取top10 marker gene
marker_genes_top10 = pd.DataFrame(
    {   
    f"{group}_{key[:1]}": result[key][group][:10] # 只提取前10个基因
    for group in groups
    for key in ["names", "logfoldchanges", "pvals"]
    }
)
marker_genes_top10.to_csv("marker_genes_top10.csv", index=False)# 注意某些基因会被excel改成日期
# 保存结果
adata.write('sce.all.int.h5ad')
