#!/usr/bin/env Rscript
# atac / rnaseq dataset
# seurat 4

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
# outDir        <- './RNA_ATAC_int_seurat4'

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
  NormalizeData(verbose = T) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(reduction.name = 'rpca') %>% 
  FindNeighbors(dims = 1:10, reduction = "rpca") %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dims = 1:10, reduction.name = 'rumap', reduction = 'rpca')

rnaseq.idents <- Idents(rnaseq.seurat)

# add ATAC dataset
rnaseq.seurat[["ATAC"]] 	<- CreateAssayObject(counts = atac)
DefaultAssay(rnaseq.seurat) <- 'ATAC'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting
VariableFeatures(rnaseq.seurat) <- rownames(rnaseq.seurat[["ATAC"]])
rnaseq.seurat <- rnaseq.seurat %>% #, normalization.method = 'CLR', margin = 2) %>% 
  NormalizeData() %>% 
  ScaleData() %>% 
  RunPCA(reduction.name = 'apca') %>%
  FindNeighbors(dims = 1:10, reduction = "apca") %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dims = 1:10, reduction.name = 'aumap', reduction = 'apca')

atac.idents <- Idents(rnaseq.seurat)

# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using bm[['weighted.nn']]
# The WNN graph can be accessed at bm[["wknn"]], 
# and the SNN graph used for clustering at bm[["wsnn"]]
# Cell-specific modality weights can be accessed at bm$RNA.weight
rnaseq.seurat <- FindMultiModalNeighbors(
  rnaseq.seurat, reduction.list = list("rpca", "apca"), 
  dims.list = list(1:30, 1:30), modality.weight.name = c("RNA.weight", 'ATAC.weight'))

rnaseq.seurat <- FindClusters(rnaseq.seurat, graph.name = "wsnn", algorithm = 1, 
	resolution = 0.6, verbose = FALSE)

#######################################
# plot integration result
rnaseq.seurat <- RunUMAP(rnaseq.seurat, nn.name = "weighted.nn", 
	reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

integrated.idents <- Idents(rnaseq.seurat)


ari <- list()
ji  <- list()

ari$rna_atac 		<- adj.rand.index(as.character(rnaseq.idents), as.character(atac.idents))
ari$rna_integrated 	<- adj.rand.index(as.character(rnaseq.idents), as.character(integrated.idents))
ari$atac_integrated <- adj.rand.index(as.character(atac.idents), as.character(integrated.idents))

ji$rna_atac 		<- cluster_similarity(rnaseq.idents, atac.idents, similarity = "jaccard")
ji$rna_integrated 	<- cluster_similarity(rnaseq.idents, integrated.idents, similarity = "jaccard")
ji$atac_integrated 	<- cluster_similarity(atac.idents, integrated.idents, similarity = "jaccard")


#######################################
# plotting low dim reps
plts <- list()
plts[["RNAseq"]]    <- DimPlot(rnaseq.seurat, reduction = "rumap", label = TRUE,
						group.by = 'RNA_snn_res.0.6')
plts[["Atac"]]      <- DimPlot(rnaseq.seurat, reduction = "aumap", label = TRUE,
						group.by = 'ATAC_snn_res.0.6')
plts[["Integrated"]]      <- DimPlot(rnaseq.seurat, reduction = "wnn.umap", label = TRUE,
						group.by = 'wsnn_res.0.6')

p.vae.grid  <- plot_grid(plotlist = plts, nrow = 1, ncol = 3, labels = names(plts))
ggsave(file.path(outDir, paste0(prefix, '_RNASEQ_ATAC_umap.pdf')), p.vae.grid,
	width = 40, height = 12, units = "cm")


# now show how the clustering in the integrated domain looks like in the input domain
plts <- list()
plts[["RNAseq"]]    <- DimPlot(rnaseq.seurat, reduction = "rumap", label = TRUE,
            group.by = 'wsnn_res.0.6')
plts[["Atac"]]      <- DimPlot(rnaseq.seurat, reduction = "aumap", label = TRUE,
            group.by = 'wsnn_res.0.6')

p.vae.grid  <- plot_grid(plotlist = plts, nrow = 1, ncol = 2, labels = names(plts))
ggsave(file.path(outDir, paste0(prefix, '_RNASEQ_ATAC_IntegratedLabels_Seurat4_umap.pdf')), p.vae.grid,
  width = 30, height = 12, units = "cm")

#######################################
# plot heatmaps of cluster similarity
scale.type <- 'row'

tab 	 <- table(rnaseq.idents, atac.idents)
RNAseq.atac.heatmap <- pheatmap(tab, cluster_rows = F, cluster_cols = F, 
    						display_numbers = T, scale = scale.type, fontsize_number = 7,
                main = 'RNAseq Vs ATAC')

tab 	 	 <- table(rnaseq.idents, integrated.idents)
rnaseq.integrated.heatmap 	 <- pheatmap(tab, cluster_rows = F, cluster_cols = F, 
    						display_numbers = T, scale = scale.type, fontsize_number = 7,
                main = 'RNAseq Vs Integrated')

tab 	 	 <- table(atac.idents, integrated.idents)
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
write.table(rnaseq.idents,       file = file.path(outDir, 'rnaseq_seurat4_clustering.csv'),     quote = F, sep = ',')
write.table(atac.idents,         file = file.path(outDir, 'atac_seurat4_clustering.csv'),       quote = F, sep = ',')
write.table(integrated.idents,   file = file.path(outDir, 'integrated_seurat4_clustering.csv'), quote = F, sep = ',')

