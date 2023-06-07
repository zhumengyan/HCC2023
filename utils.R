# 
 #Author: Zhu Mengyan
 #Date: 2022-09-28 11:10:28
 #LastEditors: Zhu Mengyan
 #LastEditTime: 2022-10-21 15:29:31
 #FilePath: /CRC2Liver_script/utils/utils.R
 #Description: 描述
 #
# 
 #Author: Zhu Mengyan
 #Date: 2021-06-19 00:02:09
 #LastEditors: Zhu Mengyan
 #LastEditTime: 2022-09-27 18:46:59
 #FilePath: /script3/utils/utils.R
 #Description: 描述
 #
#### 主要记录单细胞分析的一些函数，方便重复调用
########################################## 
# 包的加载
##########################################
set.seed(12345)
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
    #source("~/CRC2Liver/ref/ssgseaMOD.r")
}


################ 利用Seurat整合单细胞转录组样本 ####################
seurat_integration <- function(dat_list, rpca = FALSE, k.filter = 200, k.weight = 100, vars = NULL) {
    #' @param dat_list: list containing seurat object
    #' @param rpca: integration methods in Seurat, rpca or cca
    #' @param k.filter: Default 200
    #' @param k.weight: Default 100 
    #' @param vars: vars need to be regressed out, such as nCount_RNA
    library(Seurat)
    library(foreach)
    dat_list <- lapply(X = dat_list, FUN = function(x) {
        x <- Seurat::NormalizeData(x, verbose = FALSE)
        x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
    })
    features <- Seurat::SelectIntegrationFeatures(object.list = dat_list)
    
    if(rpca) {
        dat_list <- lapply(X = dat_list, FUN = function(x) {
            x <- ScaleData(x, features = features, verbose = FALSE, vars.to.regress = vars)#c("S.Score", "G2M.Score")
            x <- RunPCA(x, features = features, verbose = FALSE)
        })
        anchors <- FindIntegrationAnchors(object.list = dat_list, k.filter = k.filter, anchor.features = features, reduction = "rpca")
    } else {
        anchors <- FindIntegrationAnchors(object.list = dat_list, k.filter = k.filter, anchor.features = features)
    }
    integrated <- IntegrateData(anchorset = anchors, k.weight = k.weight)
    DefaultAssay(integrated) <- "integrated"
    integrated <- ScaleData(integrated, verbose = FALSE, vars.to.regress = vars)
    integrated <- RunPCA(integrated, verbose = FALSE, npcs = 50)
    integrated <- FindNeighbors(integrated, dims = 1:20)
    integrated <- RunTSNE(integrated, dims = 1:20)
    integrated <- RunUMAP(integrated, dims = 1:20)
    integrated
}


################ Seurat数据的差异基因 ####################
run_FindMarkers <- function(seuratObject, groups, logfc.threshold = 0.25, padj = 0.05, cores = 10) {
    #' @param seuratObject: list containing seurat object
    #' @param groups: integration methods in Seurat, rpca or cca
    #' @param logfc.threshold: Default 200
    #' @param cores: Default 100 
    #' @param vars: vars need to be regressed out, such as nCount_RNA
    set.seed(12345)
    library(Seurat)
    library(dplyr)
    DefaultAssay(seuratObject) = "RNA"
    seuratObject@active.ident = factor(seuratObject@meta.data[[groups]])
    names(seuratObject@active.ident) = colnames(seuratObject)
    n_clust = unique(seuratObject@meta.data[[groups]])
    mcFindMarkers = function(i){
        ident1 = i
        ident2 = n_clust[n_clust != i]
        #table1 <- FindMarkers(seuratObject, ident.1 = ident1, ident.2 = ident2, only.pos = TRUE,
        #    test.use = "MAST", latent.vars = c("percent.mt", "percent.ribo"))
        table1 = FindMarkers(seuratObject, ident.1 = ident1, ident.2 = ident2, only.pos = TRUE, logfc.threshold = logfc.threshold)
        table1$gene = rownames(table1)
        table1$cluster = rep(i, nrow(table1))
        return(table1)
    }
    marker_results = list()[n_clust]
    marker_results = parallel::mclapply(n_clust, mcFindMarkers, mc.cores = cores)
    markers = purrr::reduce(marker_results, rbind)
    markers = markers %>% filter(p_val < 0.05)
    mean_res = apply(seuratObject@assays$RNA@data[unique(as.character(markers$gene)),], 1, function(x) {
        tmp = data.frame(value = as.numeric(x), group = as.character(seuratObject@meta.data[[groups]]))
        tmp$group = factor(tmp$group, levels = n_clust)
        res = tmp %>% dplyr::group_by(group) %>% dplyr::summarise(mean_ = mean(value))
        res$mean_
    })
    mean_res1 = as.data.frame(t(mean_res))
    colnames(mean_res1) = str_c("AvgExp_", as.character(n_clust))
    mean_res1$gene = rownames(mean_res1)
    markers1 = dplyr::inner_join(markers, mean_res1, by = c("gene"))
    return(markers1)
}

run_FindMarkers_subsample <- function(seuratObject, groups, sample_numbers, sample_ratio, cut = 0.7, logfc.threshold = 0.25, cores = 10) {
    ## 主要目的：得到更加稳定的特征基因，直接得到的是基因列表
    #' @param seuratObject: seurat object containing different samples
    #' @param groups: groups name in meta data
    #' @param sample_numbers: 抽样次数
    #' @param sample_ratio: sample_ratio：抽样比例
    #' @param logfc.threshold: logFC, default 0.25
    #' @param cores: threads
    library(stringr)
    library(Seurat)
    DefaultAssay(seuratObject) <- "RNA"
    seuratObject@active.ident <- factor(as.character(seuratObject@meta.data[[groups]]))
    names(seuratObject@active.ident) <- colnames(seuratObject)
    n_clust <- unique(as.character(seuratObject@meta.data[[groups]]))
    DEG_list = lapply(1:sample_numbers, function(x) {
        set.seed(x)
        cells_subset = sample(colnames(seuratObject), size = ncol(seuratObject) * sample_ratio, replace=F)
        seuratObject2 = subset(seuratObject, cells = cells_subset)
        mcFindMarkers <- function(i){
            ident1 <- i
            ident2 <- n_clust[n_clust != i]
            table1 <- FindMarkers(seuratObject2, ident.1 = ident1, ident.2 = ident2, only.pos = TRUE, logfc.threshold = logfc.threshold)
            table1$gene <- rownames(table1)
            table1$cluster <- rep(i, nrow(table1))
            return(table1)
        }
        marker_results <- list()[n_clust]
        marker_results <- parallel::mclapply(n_clust, mcFindMarkers, mc.cores = cores)
        markers <- purrr::reduce(marker_results, rbind)
        markers = markers %>% filter(p_val < 0.05)
        markers$sample_number <- x
        return(markers)
    })
    markers_all <- purrr::reduce(DEG_list, rbind)
    geneList = lapply(n_clust, function(x) {
        tmp = markers_all %>% filter(cluster == x)
        genes = names(table(tmp$gene))[table(tmp$gene) > sample_numbers * cut]
        return(genes)
    })
    names(geneList) = n_clust
    return(geneList)
}

################ 利用VISION计算单细胞得分 ####################
run_vision <- function(seuratobject, geneList, latentSpace=NULL, projection_methods="UMAP", min_signature_genes = 5) {
    #' @param seuratobject: Seurat object
    #' @param geneList: gene list
    #' @param latentSpace: NULL/PCA/scVI and so on 
    #' @param projection_methods: UMAP/TSNE  
    #' @param min_signature_genes: min gene number in geneList
    #DefaultAssay(seuratobject) <- "RNA"
    library(VISION)
    sig_list <- sapply(names(geneList), function(x) {
        sigData = rep(1, length(geneList[[x]]))
        names(sigData) = geneList[[x]]
        sig <- createGeneSignature(name = x, sigData = sigData)
        sig
    })
    vis <- Vision(seuratobject, signatures=sig_list, latentSpace=latentSpace, min_signature_genes = min_signature_genes)
    options(mc.cores = 6)
    vis <- analyze(vis)
    return(vis)
}

################ 利用infercnv计算拷贝数 ####################
run_infercnv_v1.0 <- function(outdir, SeuratObject, sampleName, ct_column, aberrant_ct, ref_ct, HMM=TRUE) {
    #' @param outdir: output directory
    #' @param SeuratObject: seurat object containing different samples
    #' @param sampleName: sample name
    #' @param ct_column: cell type column
    #' @param aberrant_ct: tumor cells
    #' @param ref_ct: reference cell type
    #' @param HMM: TRUE or FALSE
    library(stringr)
	library(infercnv)
    set.seed(12345)
    dir.create(outdir)
	dir.create(str_glue(outdir, sampleName))
    scObject <- subset(SeuratObject, cells = colnames(SeuratObject)[SeuratObject$orig.ident == sampleName & SeuratObject@meta.data[[ct_column]] %in% c(aberrant_ct, ref_ct)])
    count_matrix = as.matrix(scObject@assays$RNA@counts)
	count_matrix_path = str_glue(outdir, sampleName, "/matfile.txt")
    write.table(round(count_matrix, digits=3), file = count_matrix_path, quote = F, sep = "\t")
    annotations <- data.frame(cell = colnames(count_matrix), status = as.character(scObject@meta.data[[ct_column]]))
    #annotations$status <- ifelse(annotations$status %in% aberrant_ct, "Aberrant", "Reference")
	annotations_path = str_glue(outdir, sampleName, "/anno.txt")
    write.table(annotations, file = annotations_path, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

	infercnv_obj = CreateInfercnvObject(raw_counts_matrix=count_matrix_path,
    	annotations_file=annotations_path,
    	delim="\t",
        gene_order_file="~/ref/GRCh38_position.txt",
        ref_group_names=ref_ct)
    infercnv_obj = infercnv::run(infercnv_obj,
        cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
        out_dir=str_glue(outdir, sampleName), 
        HMM = HMM, HMM_report_by = "subcluster", analysis_mode = "subclusters", denoise = TRUE, num_threads = 8)
    file.remove(
        str_glue(outdir, sampleName, "/matfile.txt"), 
        str_glue(outdir, sampleName, "/anno.txt")
    )
    # seurat_obj = infercnv::add_to_seurat(
    #     infercnv_output_path=outdir,
    #     seurat_obj=scObject, # optional
    #     top_n=10
    # )
    # saveRDS(seurat_obj, str_glue(outdir, sampleName, "/seuratObject.RDS"))
}



################ for cNMF ####################
getHVG <- function(seuratObject, name, outdir, sample_numbers, sample_ratio) {
    #' @param seuratObject: seurat object containing different samples
    #' @param name: sample name
    #' @param outdir: output directory
    #' @param sample_numbers: 抽样次数
    #' @param sample_ratio: sample_ratio：抽样比例
    library(stringr)
    library(Seurat)
    #dir.create(str_glue(outdir, "/", name))
    seuratObject1 = seuratObject
    seuratObject1 <- CreateSeuratObject(as.matrix(seuratObject1@assays$RNA@counts), min.cells = 10)
    HVG_list = lapply(1:sample_numbers, function(x) {
        set.seed(x)
        cells_subset = sample(colnames(seuratObject1), size = ncol(seuratObject1) * sample_ratio, replace=F)
        seuratObject2 = subset(seuratObject1, cells = cells_subset)
        seuratObject2 <- Seurat::FindVariableFeatures(seuratObject2, selection.method = "vst", verbose = F)
        return(seuratObject2@assays$RNA@var.features)
    })
    tmp = table(unlist(HVG_list))
    genes_selected = names(tmp)[tmp > 0.7 * sample_numbers]
    writeLines(str_c(genes_selected, collapse ="\n"), str_glue(outdir, "/", name,"_HVG.txt"), sep ="")
}
getExpData <- function(seuratObject, outdir, name) {
    library(Seurat)
    #dir.create(str_glue(outdir, "/", name))
    #seuratObject1 = subset(seuratObject, cells = colnames(seuratObject)[seuratObject$orig.ident == name])
    counts1 <- as.matrix(seuratObject@assays$RNA@counts)
    selected_genes <- rowSums(counts1>0) >= 10
    counts2 <- counts1[selected_genes, ]
    write.table(x = t(counts2), file = str_glue(outdir, "/", name, "_exp.tsv"), sep = "\t", quote = F)
}


#### 主要记录单细胞分析的一些函数，方便重复调用
run_RCAv2 <- function(seuratObject, panel) {
    set.seed(12345)
    library(RCAv2)
    RCA_from_Seurat <- RCAv2::createRCAObject(rawData = seuratObject@assays$RNA@counts, normData = seuratObject@assays$RNA@data)
    RCA_from_Seurat <- RCAv2::dataProject(rca.obj = RCA_from_Seurat, method = panel)
    RCA_from_Seurat <- RCAv2::dataSClust(RCA_from_Seurat, res = 1)
    RCA_from_Seurat <- estimateCellTypeFromProjection(RCA_from_Seurat,confidence = NULL)
    RCA_from_Seurat<-computeUMAP(RCA_from_Seurat)
    return(RCA_from_Seurat)
}
run_RCA <- function(seuratObject, panel) {
    library(stringr)
    library(dplyr)
    library(Seurat)
    library(WGCNA)
    library(flashClust)
    library(gplots)
    library(preprocessCore)
    library(RCA)
    options(stringsAsFactors = FALSE)
    set.seed(12345)
    DefaultAssay(seuratObject) <- "RNA"
    data_obj <- dataConstruct(seuratObject@assays$RNA@counts)
    data_obj1 <- geneFilt(data_obj)
    data_obj1 <- cellNormalize(data_obj, method = "scQ")
    data_obj1 <- dataTransform(data_obj1, "log10")
    data_obj1 <- featureConstruct(data_obj1, method = panel)
    data_obj1 <- cellClust(data_obj1)
    return(data_obj1)
}

importCDS <- function (otherCDS, seurat_scale=F, import_all = FALSE) 
{
    #' @param otherCDS: 
    #' @param seurat_scale: 
    #' @param import_all: 
    if (class(otherCDS)[1] == "Seurat") {
        requireNamespace("Seurat")
        if (!seurat_scale) {
          data <- otherCDS@assays$RNA@counts
        } else {
          data <- otherCDS@assays$RNA@scale.data
        }
        if (class(data) == "data.frame") {
            data <- as(as.matrix(data), "sparseMatrix")
        }
        pd <- tryCatch({
            pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
            pd
        }, error = function(e) {
            pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
            pd <- new("AnnotatedDataFrame", data = pData)
            message("This Seurat object doesn't provide any meta data")
            pd
        })
        if (length(setdiff(colnames(data), rownames(pd))) > 0) {
            data <- data[, rownames(pd)]
        }
        fData <- data.frame(gene_short_name = row.names(data), 
            row.names = row.names(data))
        fd <- new("AnnotatedDataFrame", data = fData)
        #lowerDetectionLimit <- otherCDS@is.expr
        if (all(data == floor(data))) {
            expressionFamily <- negbinomial.size()
            expr <- "negbinomial.size"
        }
        else if (any(data < 0)) {
            expressionFamily <- uninormal()
            expr <- "unimormal"
        }
        else {
            expressionFamily <- tobit()
            expr <- "tobit"
        }
        print(paste0("expressionFamily ",expr))
        # valid_data <- data[, row.names(pd)]
        monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
            #lowerDetectionLimit = lowerDetectionLimit,
            expressionFamily = expressionFamily)
        if (import_all) {
            if ("Monocle" %in% names(otherCDS@misc)) {
                otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
                otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
                monocle_cds <- otherCDS@misc$Monocle
                mist_list <- otherCDS
            }
            else {
                mist_list <- otherCDS
            }
        }
        else {
            mist_list <- list()
        }
        if ("var.genes" %in% slotNames(otherCDS)) {
            var.genes <- setOrderingFilter(monocle_cds, otherCDS@var.genes)
        }
        monocle_cds@auxClusteringData$seurat <- mist_list
    }
    else if (class(otherCDS)[1] == "SCESet") {
        requireNamespace("scater")
        message("Converting the exprs data in log scale back to original scale ...")
        data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
        fd <- otherCDS@featureData
        pd <- otherCDS@phenoData
        experimentData = otherCDS@experimentData
        if ("is.expr" %in% slotNames(otherCDS)) 
            lowerDetectionLimit <- otherCDS@is.expr
        else lowerDetectionLimit <- 1
        if (all(data == floor(data))) {
            expressionFamily <- negbinomial.size()
        }
        else if (any(data < 0)) {
            expressionFamily <- uninormal()
        }
        else {
            expressionFamily <- tobit()
        }
        if (import_all) {
            mist_list <- otherCDS
        }
        else {
            mist_list <- list()
        }
        monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
            lowerDetectionLimit = lowerDetectionLimit, expressionFamily = expressionFamily)
        monocle_cds@auxOrderingData$scran <- mist_list
    }
    else {
        stop("the object type you want to export to is not supported yet")
    }
    return(monocle_cds)
}

## monocle2进行轨迹分析
monocle2_run <- function(seuratObject, cluster, deg_number, cds = NULL, reset_root = NULL) {
    #' @param seuratObject: seurat object
    #' @param cluster: pre-defined clusters, used to find ordering genes
    #' @param deg_number: topN genes
    #' @param cds: if you want to reset start point, you should provide monocle result that you got before
    #' @param reset_root: reset root, state number
    library(monocle)
    library(stringr)
    if (is.null(reset_root) & is.null(cds)) {
        DefaultAssay(seuratObject) = "RNA"
        dat = Seurat::as.CellDataSet(seuratObject)
        dat = estimateSizeFactors(dat)
        dat = estimateDispersions(dat)
        ##过滤低质量的细胞 
        dat = detectGenes(dat, min_expr = 0.1)
        fData(dat)$use_for_ordering = fData(dat)$num_cells_expressed > 0.05 * ncol(dat)
        expressed_genes = row.names(subset(fData(dat),num_cells_expressed >= 10))
        # dat <- reduceDimension(dat,max_components = 2,
        #     norm_method = 'log',num_dim = 20,reduction_method = 'tSNE',
        #     verbose = T,check_duplicates=F
        # )
        # dat <- clusterCells(dat,verbose = F)
        clustering_DEG_genes = differentialGeneTest(dat[expressed_genes,],
                fullModelFormulaStr = str_c("~", cluster), cores = 10)
        ordering_genes = row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:deg_number]
        dat = setOrderingFilter(dat, ordering_genes = ordering_genes)
        dat = reduceDimension(dat, method = 'DDRTree')
        dat = orderCells(dat)
    }
    if (!is.null(reset_root) & !is.null(cds)) {
        cds = orderCells(cds, root_state = reset_root)
        dat = cds
    }
    return(dat)
}



################ Seurat标准流程 ####################
{
    run_Seurat <- function(dat, resolution = 1, vars = NULL) {
        #' @param dat: seurat object
        #' @param resolution: resolution for clustering
        #' @param vars: Variables to regress out
        #' @return seurat object
        
        dat <- NormalizeData(object = dat, normalization.method = "LogNormalize", 
            scale.factor = 10000)
        dat <- FindVariableFeatures(object = dat, selection.method = "vst", nfeatures = 2000)
        dat <- ScaleData(object = dat, vars.to.regress = vars)
        dat <- RunPCA(object = dat, features = VariableFeatures(object = dat))
        dat <- FindNeighbors(object = dat, dims = 1:20, force.recalc = TRUE)
        dat <- FindClusters(object = dat, resolution = resolution)
        dat <- RunUMAP(object = dat, dims = 1:20)
        return(dat)
    }    
}

{
    add_percentage <- function(dat) {
        ribo.genes <- read.table("~/ref/rRNA_genes_RPG.txt")
        dat[["percent.ribo"]] <- PercentageFeatureSet(object = dat, 
            features = rownames(dat) %in% as.character(ribo.genes$V1))
        dat[["percent.mt"]] <- PercentageFeatureSet(object = dat, pattern = "^MT-")
        return(dat)
    }
}


{
    seurat2loom <- function(inRDS, outloom){
        dat.loom <- SeuratDisk::as.loom(inRDS, filename = outloom, assay = "RNA", verbose = FALSE, overwrite = TRUE)
        dat.loom$close_all()
    }
}






convertMouseGeneList <- function(x){
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = getLDS(attributes = c("mgi_id"), filters = "mgi_symbol", 
        values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T
    )
    row.names(genesV2) <- genesV2$HGNC.symbol
    return(genesV2)
}
# dat <- read_delim(gzfile("./ref/Extra_data/GEO/GSE155364/GSE155364_20200715_Sfrp1_tumor_cells_invitro_invivo.txt.gz"), delim = "\t")
# dat1 <- dat[,c(1, 3:14)]
# convertMouseGeneList(dat1$ensembl_gene_id)










####################### draw boxplot of cibersortX result ####################### 
draw_cs_box <- function(mat, status){
    #' @param mat: cibersort结果
    #' @param status: 分组，两组
    dat <- cbind(mat, status)
    dat <- dat[,-c(1,24:26)]
    pvalues_ <- apply(dat[, -c(23)], 2, function(x){
        tmp = data.frame(
        disease = status,
        value = x
        )
        res <- wilcox.test(value ~ disease, data = tmp)
        return(res$p.value)
    })
        dat1 <- tidyr::pivot_longer(data = dat, cols= -status, names_to = "Celltype",
            values_to = "Fraction")
        dat1$status <- factor(dat1$status, levels=c(1, 0))
    pv <- data.frame(celltype = unique(dat1$Celltype), pvalue = pvalues_)
    pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                        labels=c('***', '**', '*', ' ', ' '))
    p <- ggplot(data = dat1, aes(x = Celltype, y = Fraction)) + 
        geom_boxplot(aes(fill = status), outlier.shape = NA) + scale_fill_manual(values = c("#DE5F5A", "#4C8BAD")) +
        geom_text(aes(x=pv$celltype, y=max(dat1$Fraction) * 1.1, 
                    label=pv$sigcode), data = pv) + 
        theme(axis.text.x = element_text(angle=45,size=8, hjust = 1.0))
    return(p)
}

draw_cs_box_single <- function(mat, status) {
    library(ggpubr)
    library(ggplot2)
    dat <- mat[,-c(1,24:26)]
    plots = list()
    for (i in 1:ncol(dat)) {
        tmp = data.frame(
            value = as.numeric(dat[,i]),
            group = status
        )
        tmp$group = factor(tmp$group, levels=c(1,0))
        plots[[i]] <- ggplot(tmp, aes(x = group, y = value)) + geom_boxplot(aes(fill = group)) +
            ylab("Fraction") + xlab('group') + theme_classic() + ggtitle(colnames(dat)[i])
    }
    return(plots)
}

BaCo <- function(X){
  
  alpha0 <- rep(1/nrow(X),ncol(X))
  beta0=1-alpha0
  nrowsX <- nrow(X)
  k <- ncol(X)
  cs <- colSums(X)
  alphas <- alpha0 + X
  betas  <- matrix(rep(beta0,nrowsX), nrow=nrowsX, byrow=TRUE) + matrix(rep(cs,nrowsX), nrow=nrowsX, byrow=TRUE) - X
  alphasPLUSbetas <- alphas + betas
  Psi <- alphas/alphasPLUSbetas - matrix(rep(rowSums(alphas/alphasPLUSbetas)/k, k), ncol=k, byrow=FALSE) 
  var_vec <- as.matrix( ( rowSums( (alphas*betas)/( (alphasPLUSbetas^2)*(alphasPLUSbetas+1) ) ) + rowSums(Psi^2) )/k )
  cov_mtrx <- (Psi %*% t(Psi))/k
  Bcorrvals <- cov_mtrx / sqrt( var_vec %*% t(var_vec) )
  diag(Bcorrvals) <- 1
  Bcorrvals
}



