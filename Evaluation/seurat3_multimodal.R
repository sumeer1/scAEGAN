#!/usr/bin/env Rscript
# atac / rnaseq dataset
# use seurat3 

# conda activate seurat4
# export R_LIBS_USER="/home/lehmanr/R4_libs"

########################################
# input

library(dplyr,     quietly = T)
library(optparse,  quietly = T)
library(Seurat,    quietly = T)
library(tidyverse, quietly = T)
library(cowplot,   quietly = T)
library(pdfCluster, quietly = T)
library(pheatmap,  quietly = T)
library(clusteval, quietly = T)

########################################
# example
# Rscript 
 
option_list = list(
  make_option(c("-r", "--rnaseq"), type="character", default=NULL, 
              help="rnaseq file name", metavar="character"),
  make_option(c("-a", "--atac"), type="character", default=NULL, 
              help="atac file name", metavar="character"),
    make_option(c("-p", "--prefix"), type="character", default="", 
              help="file prefix", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", 
              help="output dir name", metavar="character")
); 
 
opt_parser  <- OptionParser(option_list=option_list);
opt     <- parse_args(opt_parser);

if (is.null(opt$rnaseq)) {
  print_help(opt_parser)
  stop("rnaseq (-r) dataset file required", call.=FALSE)
}
if (is.null(opt$atac)) {
  print_help(opt_parser)
  stop("atac (-a) dataset file required", call.=FALSE)
}
if (is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("output directory (-o) required", call.=FALSE)
}


rnaseq.file     <- opt$rnaseq
atac.file       <- opt$atac
outDir          <- opt$outdir
prefix          <- opt$prefix


# permutation test results
# # the input rnaseq and atac data
# rnaseq.file    <- 'RNA.csv.gz'
# atac.file    <- 'ATAC.csv.gz'
# outDir        <- './RNA_ATAC_int_seurat3'
# 
########################################


########################################
# load data
########################################
cat('Loading data...\n')
print(Sys.time())

rnaseq <- read.table(rnaseq.file, sep = ',', header = T, 
            row.names = 1, stringsAsFactors = F)

atac   <- read.table(atac.file, sep = ',', header = T, 
          row.names = 1, stringsAsFactors = F)

# make sure the output directory exists
if(!dir.exists(outDir)) {
	cat('creating output dir: ',outDir,'\n')
	dir.create(outDir, recursive = T)
}

#######################################
# do seurat processing for plotting etc
cat('Making Seurat Objects...\n')
print(Sys.time())


# RNA-seq dataset
rnaseq.seurat <- CreateSeuratObject(
  counts = rnaseq, project = "rnaseq",
  assay = "RNA", min.cells = 0, min.features = 0, names.field = 1)

rnaseq.seurat <- rnaseq.seurat %>% 
  NormalizeData(verbose = T)  %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:10) %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dims = 1:10, reduction = 'pca')

rnaseq.idents <- Idents(rnaseq.seurat)

# add ATAC dataset
atac.seurat <- CreateSeuratObject(
  counts = atac, project = "atac",
  assay = "RNA", min.cells = 0, min.features = 0, names.field = 1)

atac.seurat <- atac.seurat %>% 
  NormalizeData(verbose = T)  %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:10) %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dims = 1:10, reduction = 'pca')

atac.idents <- Idents(atac.seurat)

#################
# integration
anchors <- FindIntegrationAnchors(object.list = 
  list(A = rnaseq.seurat, B = atac.seurat), dims = 1:30)
A.B.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(A.B.integrated) <- "integrated"

A.B.integrated <- A.B.integrated %>%
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:10)

# do label transfer from integrated to input domains
# integrated to RNAseq
trans.anchors <- FindTransferAnchors(reference = A.B.integrated, query = rnaseq.seurat, 
    dims = 1:30)
predictions <- TransferData(anchorset = trans.anchors, refdata = A.B.integrated$seurat_clusters, 
    dims = 1:30)
rnaseq.seurat <- AddMetaData(rnaseq.seurat, metadata = predictions)

# integrated to ATAC
trans.anchors <- FindTransferAnchors(reference = A.B.integrated, query = atac.seurat, 
    dims = 1:30)
predictions <- TransferData(anchorset = trans.anchors, refdata = A.B.integrated$seurat_clusters, 
    dims = 1:30)
atac.seurat <- AddMetaData(atac.seurat, metadata = predictions)


#######################################
# calc clustering similarities
ari <- list()
ji  <- list()

ari$rna_atac 		<- adj.rand.index(as.character(rnaseq.seurat$seurat_clusters), 
  as.character(atac.seurat$seurat_clusters))
ari$rna_integrated 	<- adj.rand.index(as.character(rnaseq.seurat$seurat_clusters), 
  as.character(rnaseq.seurat$predicted.id))
ari$atac_integrated <- adj.rand.index(as.character(atac.seurat$seurat_clusters), 
  as.character(atac.seurat$predicted.id))

ji$rna_atac 		<- cluster_similarity(rnaseq.seurat$seurat_clusters, atac.seurat$seurat_clusters, 
  similarity = "jaccard")

ji$rna_integrated 	<- cluster_similarity(as.factor(rnaseq.seurat$seurat_clusters),
  as.factor(rnaseq.seurat$predicted.id), 
  similarity = "jaccard")

ji$atac_integrated 	<- cluster_similarity(as.factor(atac.seurat$seurat_clusters), 
  as.factor(atac.seurat$predicted.id), 
  similarity = "jaccard")


#######################################
# plotting low dim reps
# first show independent clusterings of each domain
plts <- list()
plts[["RNAseq"]]    <- DimPlot(rnaseq.seurat, reduction = "umap", label = TRUE,
						group.by = 'RNA_snn_res.0.6')
plts[["Atac"]]      <- DimPlot(atac.seurat, reduction = "umap", label = TRUE,
						group.by = 'RNA_snn_res.0.6')
plts[["Integrated"]]      <- DimPlot(A.B.integrated, reduction = "umap", label = TRUE,
						group.by = 'RNA_snn_res.0.6')

p.vae.grid  <- plot_grid(plotlist = plts, nrow = 1, ncol = 3, labels = names(plts))
ggsave(file.path(outDir, paste0(prefix, '_RNASEQ_ATAC_IntegratedSeurat3_independentClustering_umap.pdf')), 
  p.vae.grid, width = 40, height = 12, units = "cm")

# now show how the clustering in the integrated domain looks like in the input domain
plts <- list()
plts[["RNAseq"]]    <- DimPlot(rnaseq.seurat, reduction = "umap", label = TRUE,
            group.by = 'predicted.id')
plts[["Atac"]]      <- DimPlot(atac.seurat, reduction = "umap", label = TRUE,
            group.by = 'predicted.id')

p.vae.grid  <- plot_grid(plotlist = plts, nrow = 1, ncol = 2, labels = names(plts))
ggsave(file.path(outDir, paste0(prefix, '_RNASEQ_ATAC_IntegratedLabels_Seurat3_umap.pdf')), p.vae.grid,
  width = 30, height = 12, units = "cm")


#######################################
# plot heatmaps of cluster similarity
scale.type <- 'row'

tab 	 <- table(rnaseq.seurat$seurat_clusters, atac.seurat$seurat_clusters)
RNAseq.atac.heatmap <- pheatmap(tab, cluster_rows = F, cluster_cols = F, 
    						display_numbers = T, scale = scale.type, fontsize_number = 7,
                main = 'RNAseq Vs ATAC')

tab 	 	 <- table(rnaseq.seurat$seurat_clusters, rnaseq.seurat$predicted.id)
rnaseq.integrated.heatmap 	 <- pheatmap(tab, cluster_rows = F, cluster_cols = F, 
    						display_numbers = T, scale = scale.type, fontsize_number = 7,
                main = 'RNAseq Vs Integrated')

tab 	 	 <- table(atac.seurat$seurat_clusters, atac.seurat$predicted.id)
atac.integrated.heatmap 	 <- pheatmap(tab, cluster_rows = F, cluster_cols = F, 
    						display_numbers = T, scale = scale.type, fontsize_number = 7,
                main = 'ATAC Vs Integrated')


p.vae.grid  <- plot_grid(plotlist = list(RNAseq.atac.heatmap[[4]], rnaseq.integrated.heatmap[[4]],
	atac.integrated.heatmap[[4]]), nrow = 1, ncol = 3, label_size = 10)
ggsave(file.path(outDir, paste0(prefix, '_ingrated_vs_inputClass_heatmap.pdf')), p.vae.grid,
	width = 30, height = 10, units = "cm")

#######################################
# similarity before integration
cat("########## RESULT ##########\n")
cat("ARI: \n")
print(ari)
cat("JI: \n")
print(ji)
cat("\n############################\n")

# write out clusterings for further analyses
write.table(rnaseq.seurat$seurat_clusters,       file = file.path(outDir, 'rnaseq_seurat3_clustering.csv'),     quote = F, sep = ',')
write.table(atac.seurat$seurat_clusters,         file = file.path(outDir, 'atac_seurat3_clustering.csv'),       quote = F, sep = ',')
write.table(rnaseq.seurat$predicted.id,   file = file.path(outDir, 'rnaseq_transfered_seurat3_clustering.csv'), quote = F, sep = ',')
write.table(atac.seurat$predicted.id,     file = file.path(outDir, 'atac_transfered_seurat3_clustering.csv'), quote = F, sep = ',')
  