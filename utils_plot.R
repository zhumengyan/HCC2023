# 
 #Author: Zhu Mengyan
 #Date: 2022-09-28 11:10:38
 #LastEditors: Zhu Mengyan
 #LastEditTime: 2022-10-09 15:04:32
 #FilePath: /CRC2Liver_script/utils/utils_plot.R
 #Description: 描述
 #
# 
 #Author: Zhu Mengyan
 #Date: 2022-08-23 21:10:30
 #LastEditors: Zhu Mengyan
 #LastEditTime: 2022-09-15 14:12:44
 #FilePath: /script3/utils/utils_plot.R
 #Description: 绘制图形的函数
 #

########################################## 
# 包的加载
##########################################
{
    library(Seurat)
    library(forcats)
    library(dplyr)
    library(stringr)
    library(ggplot2)
    library(RColorBrewer)
    library(cowplot)
    library(ggsci)
    library(infercnv)
    library(readr)
    library(tidyr)
    library(tidyverse)
}

########################################## 
# 绘图的颜色配置
########################################## 
{
    colors_ggsci_npg <- c(
        "#4DBBD5", "#7E6148","#E64B35", "#8491B4", "#F39B7F", "#00A087", "#91D1C2",  "#B09C85","#DC0000","#3C5488"
    )
    mycolors = c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
        '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999',
        "#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51", "#457b9d", "#a8dadc", "#81b29a", "#f2cc8f")    

    mycolors = c(
        '#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC',
        '#54B0E4','#A65628','#B3DE69','#222F75','#1B9E77', "#81b29a",'#B2DF8A',
        '#00CDD1','#E3BE00','#FB9A99',
        '#E7298A','#910241','#A6CEE3','#CE1261','#5E4FA2','#8CA77B',
        '#00441B','#DEDC00','#8DD3C7','#999999',"#264653", "#2a9d8f", 
        "#e9c46a", "#f4a261", "#e76f51", "#457b9d", "#a8dadc", "#f2cc8f"
    )    

    selected_colors=c("#e63946", "#a8dadc", "#457b9d", "#1d3557", "#f2cc8f", "#81b29a", 
        "#e07a5f", "#ffcb77", "#17c3b2", "#0081a7", "#fed9b7", "#c1121f", "#f4a261", "#2a9d8f",
        "#adc178", "#dde5b6", "#006d77", "#80ed99", "#f2542d", "#edddd4")




    #' code from Cai Yun
    #' Determine the color scheme. Can be specified for samples or clusters to avoid confusion. 
    #' @param type Type of scheme ("samples" or "clusters").
    #' @return A vector of colors.
    #' @importFrom RColorBrewer brewer.pal
    #' @importFrom ggsci pal_d3 pal_igv
    get_color_scheme <- function(type = "clusters") {
        if (type == "samples") {
            color_scheme <- c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
        }
        if (type == "clusters") {
            color_scheme <- unique(c(pal_d3("category20")(20), rev(pal_d3("category20b")(20)), 
                                    brewer.pal(9, "Set1"), brewer.pal(12, "Set3"),brewer.pal(8, "Set2"),
                                    brewer.pal(8, "Dark2"),brewer.pal(8, "Accent"),brewer.pal(12, "Paired"), 
                                    brewer.pal(9, "Pastel1"),brewer.pal(8, "Pastel2"),
                                    pal_igv("default")(51)))
        }
        return(color_scheme)
    }
}
colors_selected <- c("#1f946d","#eb8f42", "#b15f2f", "#e61f2f","#2579b0", "#8f9192","#2eaa54")

#theme_set(theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1))
theme_global <- theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)

########################################## 
# 绘制seurat object的常见umap图， 包括分组、样本、单个样本、细胞类别
########################################## 
#' Determine the point size for UMAP plots (smaller for larger datasets).
#' @param num_cells Number of cells (points on the plot).
#' @return Numeric point size.
get_dr_pt_size <- function(num_cells) {

    pt_size <- 1.8
    if (num_cells > 1000) pt_size <- 1.2
    if (num_cells > 5000) pt_size <- 1.0
    if (num_cells > 10000) pt_size <- 0.8
    if (num_cells > 25000) pt_size <- 0.6
    if (num_cells > 50000) pt_size <- 0.4
    if (num_cells > 100000) pt_size <- 0.2
    return(pt_size)
}



plot_seurat_summary <- function(seuratObject, outdir, prefix="test", colors="cols", cluster="seurat_clusters", group="groups", phase="Phase", sample="orig.ident", batch="batch") {
    #' 绘制seurat object的常见umap图
    #' @param seuratObject: Seurat object
    #' @param outdir: 输出目录
    #' @param prefix: 文件名前缀
    #' @param colors: 颜色表
    #' @param cluster: cluster
    #' @param group: 分组信息
    #' @param phase: cell cycle phase
    #' @param sample: 样本
    #' @param batch: batch effect
    library(Seurat)
    library(stringr)
    filepath = str_glue(outdir, prefix, "_")
    pt.size = get_dr_pt_size(ncol(seuratObject))
    pdf(str_glue(filepath, "clusters.pdf"), useDingbats=F)
    p <- DimPlot(seuratObject, reduction="umap", group.by=cluster, pt.size=pt.size, cols=colors)  +
        theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
    print(p)
    dev.off()    
    pdf(str_glue(filepath, "groups.pdf"), useDingbats=F)
    p <- DimPlot(seuratObject, reduction="umap", group.by=group, pt.size=pt.size)  + scale_color_manual("Groups", values = c("#5b8e7d", "#e09f3e", "#ef233c")) +
        theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
    print(p)
    dev.off()
    pdf(str_glue(filepath, "groups2.pdf"), useDingbats=F)
    for (i in unique(seuratObject@meta.data[[group]])) {
        p <- DimPlot(seuratObject, reduction="umap", cells.highlight = colnames(seuratObject)[seuratObject@meta.data[[group]] == i], pt.size=pt.size)  + 
            scale_color_manual("Groups", values = c("grey", "red")) + ggtitle(i) +
            theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
        print(p) 
    }
    dev.off()  
    pdf(str_glue(filepath, "phase.pdf"), useDingbats=F)
    p <- DimPlot(seuratObject, reduction="umap", group.by=phase, pt.size=pt.size) +
        theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
    print(p)
    dev.off()
    pdf(str_glue(filepath, "sample.pdf"), useDingbats=F)
    p <- DimPlot(seuratObject, reduction="umap", group.by=sample, pt.size=pt.size) +
        theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
    print(p)
    dev.off()    
    pdf(str_glue(filepath, "sample2.pdf"), useDingbats=F)
    for (i in unique(seuratObject@meta.data[[sample]])) {
        p <- DimPlot(seuratObject, reduction="umap", cells.highlight = colnames(seuratObject)[seuratObject@meta.data[[sample]] == i], pt.size=pt.size)  + 
            scale_color_manual("samples", values = c("grey", "red")) + ggtitle(i) +
            theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
        print(p) 
    }
    dev.off()  
    pdf(str_glue(filepath, "batch.pdf"), useDingbats=F)
    p <- DimPlot(seuratObject, reduction="umap", group.by=batch, pt.size=pt.size) +
        theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
    print(p)
    dev.off()      

}




 
    # pdf("./result/integration/integrated_cca_celltype.pdf", useDingbats=F)
    # colors_ggsci_npg <- c(
    #     "#4DBBD5", "#7E6148","#E64B35", "#F39B7F", "#00A087",  "#8491B4","#B09C85","#91D1C2", "#DC0000","#3C5488"
    # )
    # #scales::show_col(colors_ggsci_npg)
    # DimPlot(integrated, reduction="umap", group.by="celltype_markers", pt.size=0.75, label=T, cols=colors_ggsci_npg) +
    #     theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
    # dev.off()   






