#!/usr/bin/env Rscript
# atac / rnaseq dataset

########################################
# input

library(optparse,  quietly = T)
library(Seurat,    quietly = T)
library(tidyverse, quietly = T)
library(cowplot,   quietly = T)
library(pdfCluster, quietly = T)
library(pheatmap,  quietly = T)
library(clusteval, quietly = T)

########################################
# example
# Rscript evalCycleGAN_Unpaired.R --dA data/Autoencoder/RNA.csv.gz --dB data/Autoencoder/ATAC.csv.gz --dALatent data/CycleGan/AE_RNA.csv.gz --dBLatent data/CycleGan/AE_ATAC.csv.gz --A2B data/CycleGan/AE_AtoB.csv.gz --B2A data/CycleGan/AE_BtoA.csv.gz -o tmp
 
option_list = list(
  make_option(c("--dA"), type="character", default=NULL, 
              help="domain A data file name", metavar="character"),
  make_option(c("--dB"), type="character", default=NULL, 
              help="domain B data file name", metavar="character"),
  make_option(c("--dALatent"), type="character", default=NULL, 
              help="A latent space file name", metavar="character"),
  make_option(c("--dBLatent"), type="character", default=NULL, 
              help="B latent space file name", metavar="character"),
  make_option(c("--A2B"), type="character", default=NULL, 
              help="A2B file name", metavar="character"),
  make_option(c("--B2A"), type="character", default=NULL, 
              help="B2A file name", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default="", 
              help="file prefix", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", 
              help="output dir name", metavar="character")
); 
 
opt_parser  <- OptionParser(option_list=option_list);
opt     <- parse_args(opt_parser);

print(opt)

if (is.null(opt$dA)) {
  print_help(opt_parser)
  stop("domain A (-dA) dataset file required", call.=FALSE)
}
if (is.null(opt$dB)) {
  print_help(opt_parser)
  stop("domain B (-dB) dataset file required", call.=FALSE)
}
if (is.null(opt$dALatent)) {
  print_help(opt_parser)
  stop("A Latent (-ALatent) dataset file required", call.=FALSE)
}
if (is.null(opt$dBLatent)) {
  print_help(opt_parser)
  stop("B Latent (-BLatent) dataset file required", call.=FALSE)
}
if (is.null(opt$A2B)) {
  print_help(opt_parser)
  stop("A2B (--A2B) dataset file required", call.=FALSE)
}
if (is.null(opt$B2A)) {
  print_help(opt_parser)
  stop("B2A (--B2A) dataset file required", call.=FALSE)
}
if (is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("output directory (-o) required", call.=FALSE)
}


a.file     <- opt$dA
b.file     <- opt$dB
ae.a.file  <- opt$dALatent
ae.b.file  <- opt$dBLatent
A2B.file   <- opt$A2B
B2A.file   <- opt$B2A
outDir     <- opt$outdir
prefix     <- opt$prefix








# permutation test results
# ae.a.file    <- ''
# ae.b.file    <- ''
# A2B.file  <- ''
# B2A.file  <- ''
# # the input rnaseq and atac data to the AE
# a.file    <- ''
# b.file    <- ''
# outDir        <- ''

########################################


########################################
# load data
########################################
cat('Loading data...\n')
print(Sys.time())

ae.A	<- (read.table(ae.a.file, sep = ',', header = T, 
            row.names = 1, stringsAsFactors = F))

ae.B 	<- (read.table(ae.b.file, sep = ',', header = T, 
          row.names = 1, stringsAsFactors = F))

A2B 	<-  read.table(A2B.file, sep = ',', header = T, 
          row.names = 1, stringsAsFactors = F)

B2A   <-  read.table(B2A.file, sep = ',', header = T, 
          row.names = 1, stringsAsFactors = F)

a <- (read.table(a.file, sep = ',', header = T, 
            row.names = 1, stringsAsFactors = F))

b   <- (read.table(b.file, sep = ',', header = T, 
          row.names = 1, stringsAsFactors = F))

# make sure the output directory exists
if(!dir.exists(outDir)) {
	cat('creating output dir: ',outDir,'\n')
	dir.create(outDir, recursive = T)
}

#######################################
# do seurat processing for plotting etc
cat('Making Seurat Objects...\n')
print(Sys.time())

# AE compressed RNA-seq dataset
ae.A.seurat <- CreateSeuratObject(
  counts = ae.A, project = "ae.A",
  assay = "RNA", min.cells = 0, min.features = 0, names.field = 1)

ae.A.seurat <- ae.A.seurat %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(npcs = 20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.6, algorithm = 1) %>% 
  RunUMAP(dims = 1:20)


# AE compressed ATAC dataset
ae.B.seurat <- CreateSeuratObject(
  counts = ae.B, project = "ae.B",
  assay = "RNA", min.cells = 0, min.features = 0, names.field = 1)
ae.B.seurat <- ae.B.seurat %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(npcs = 20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.6, algorithm = 1) %>%
  RunUMAP(dims = 1:20)


# A2B
A2B.seurat <- CreateSeuratObject(
  counts = A2B, project = "A2B",
  assay = "RNA", min.cells = 0, min.features = 0, names.field = 1)
A2B.seurat <- A2B.seurat %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.6) %>%
  RunUMAP(dims = 1:20)

# add rnaseq clustering to translated dataset
A2B.seurat <- A2B.seurat %>% 
	AddMetaData(Idents(ae.A.seurat), col.name = 'RNAseqID')

# B2A
B2A.seurat <- CreateSeuratObject(
  counts = B2A, project = "B2A",
  assay = "RNA", min.cells = 0, min.features = 0, names.field = 1)
B2A.seurat <- B2A.seurat %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dims = 1:20)

# add rnaseq clustering to translated dataset
B2A.seurat <- B2A.seurat %>% 
  AddMetaData(Idents(ae.B.seurat), col.name = 'AtacID')

# RNA-seq dataset
a.seurat <- CreateSeuratObject(
  counts = a, project = "A",
  assay = "RNA", min.cells = 0, min.features = 0, names.field = 1)

a.seurat <- a.seurat %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dims = 1:20)

# ATAC dataset
b.seurat <- CreateSeuratObject(
  counts = b, project = "B",
  assay = "RNA", min.cells = 0, min.features = 0, names.field = 1)
b.seurat <- b.seurat %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dims = 1:20)

#######################################
# plotting all
plts <- list()
plts[["RNAseq"]]    <- DimPlot(a.seurat,    reduction = "umap", label = TRUE)
plts[["Atac"]]      <- DimPlot(b.seurat,    reduction = "umap", label = TRUE)
plts[["AE_RNAseq"]] <- DimPlot(ae.A.seurat, reduction = "umap", label = TRUE)
plts[["AE_Atac"]]   <- DimPlot(ae.B.seurat, reduction = "umap", label = TRUE)
plts[["A2B"]]       <- DimPlot(A2B.seurat,  reduction = "umap", label = TRUE)
plts[["B2A"]]       <- DimPlot(B2A.seurat,  reduction = "umap", label = TRUE)

idx <- c(1,3,5,2,4,6)
p.vae.grid  <- plot_grid(plotlist = plts[idx], nrow = 2, ncol = 3, labels = names(plts)[idx])
ggsave(file.path(outDir, paste0( 'CycleGan_in_vs_reconst1.pdf')), p.vae.grid,
	width = 30, height = 16, units = "cm")

# plot heatmaps of cluster similarity
scale.type <- 'row'
ae.B2A.vs.RNAseq.tab 	 <- table(Idents(ae.A.seurat), Idents(B2A.seurat))
ae.B2A.vs.RNAseq.heatmap <- pheatmap(ae.B2A.vs.RNAseq.tab, cluster_rows = F, cluster_cols = F, 
    						display_numbers = T, scale = scale.type, fontsize_number = 7,
                main = 'AE RNAseq Vs B2A')

ae.A2B.vs.atac.tab 	 	 <- table(Idents(ae.B.seurat), Idents(A2B.seurat))
ae.A2B.vs.atac.heatmap 	 <- pheatmap(ae.A2B.vs.atac.tab, cluster_rows = F, cluster_cols = F, 
    						display_numbers = T, scale = scale.type, fontsize_number = 7,
                main = 'AE ATAC Vs A2B')

p.vae.grid  <- plot_grid(plotlist = list(ae.B2A.vs.RNAseq.heatmap[[4]], ae.A2B.vs.atac.heatmap[[4]]), 
					nrow = 1, ncol = 2, label_size = 10)
ggsave(file.path(outDir, paste0( 'cGANtranslated_vs_inputClass_heatmap1.pdf')), p.vae.grid,
	width = 22, height = 10, units = "cm")

########################################
# cross classification precision and correlation
########################################
cat('SVM Label Transfer...\n')
print(Sys.time())

# classify generated cells into populations using a SVM
classifyNew.SVM <- function(trainingData, cellIds, newData, scale = T, kernel = 'linear') {
    require(e1071)

    cat('training data: ' , nrow(trainingData), 'genes / ', ncol(trainingData), ' cells\n')
    cat('new data: ' , nrow(newData), 'genes / ', ncol(newData), ' cells\n')
    svm.data    <- as.data.frame(t(trainingData))
    newData         <- t(newData)

    # depending on whether the cell ids dataframe or the seurat object Idents are provided
    if(class(cellIds) == 'factor') {
      svm.data$pop    <- cellIds[rownames(svm.data)]
    } else if(class(cellIds) == 'data.frame') {
      svm.data$pop    <- cellIds[rownames(svm.data), 'pop']
    } else {
      stop('unclear cellIds format...')
    }

    # svm.data$pop    <- relevel(svm.data$pop, ref = "6_1")
    cat('training SVM...\n')
    ml.model        <- svm(formula = pop ~ .,
                 data = svm.data, 
                 type = 'C-classification', 
                 kernel = kernel,
                 scale = scale)
    cat('predicting...\n')
    svm.prediction  <- predict(ml.model, newdata = newData, type = 'class')
    res             <- data.frame(pop = svm.prediction, row.names = rownames(newData))
    return(res)
}

############################
cat('svm cross class ...\n')

joint.var.feats <- intersect(intersect(intersect(VariableFeatures(ae.A.seurat), VariableFeatures(ae.B.seurat)),
  VariableFeatures(A2B.seurat)), VariableFeatures(B2A.seurat))
cat('found', length(joint.var.feats),' joint variable features...\n')
feature.idx.A <- match(joint.var.feats, rownames(ae.A.seurat@assays$RNA@scale.data))
feature.idx.B <- match(joint.var.feats, rownames(ae.B.seurat@assays$RNA@scale.data))
feature.idx.A2B <- match(joint.var.feats, rownames(A2B.seurat@assays$RNA@scale.data))
feature.idx.B2A <- match(joint.var.feats, rownames(B2A.seurat@assays$RNA@scale.data))

# label transfer BEFORE integration
cellIds.Pred <- classifyNew.SVM(ae.A.seurat@assays$RNA@scale.data[feature.idx.A,], Idents(ae.A.seurat), 
  as.data.frame(ae.B.seurat@assays$RNA@scale.data[feature.idx.B,]), 
  kernel = 'linear', scale = F)
ae.A.vs.B.adj.rand  <- adj.rand.index(Idents(ae.B.seurat),     cellIds.Pred$pop)
ae.A.vs.B.ji        <- cluster_similarity(Idents(ae.B.seurat), cellIds.Pred$pop)

cellIds.Pred <- classifyNew.SVM(ae.B.seurat@assays$RNA@scale.data[feature.idx.B,], Idents(ae.B.seurat), 
  as.data.frame(ae.A.seurat@assays$RNA@scale.data[feature.idx.A,]), 
  kernel = 'linear', scale = F)
ae.B.vs.A.adj.rand  <- adj.rand.index(Idents(ae.A.seurat), cellIds.Pred$pop)
ae.B.vs.A.ji        <- cluster_similarity(Idents(ae.A.seurat), cellIds.Pred$pop)

# label transfer AFTER integration
cat('check integration result A to B2A...\n')
cellIds.Pred <- classifyNew.SVM(ae.A.seurat@assays$RNA@scale.data[feature.idx.A,], Idents(ae.A.seurat), 
  as.data.frame(B2A.seurat@assays$RNA@scale.data[feature.idx.B2A,]), 
  kernel = 'linear', scale = F)
ae.A.vs.B2A.adj.rand  <- adj.rand.index(Idents(ae.A.seurat), cellIds.Pred$pop)
ae.A.vs.B2A.ji        <- cluster_similarity(Idents(ae.A.seurat), cellIds.Pred$pop)

# inA.vs.B2.class <- table(Idents(ae.A.seurat), cellIds.Pred$pop)
# cross.class.plt  <- pheatmap(inA.vs.generatedA.class, cluster_rows = F, cluster_cols = F, 
#     scale = 'none', display_numbers = T)

# ggsave(file.path(outDir, paste0(prefix, 'domain1_to_2_crossClassSVM.pdf')), cross.class.plt,
#   width = 10, height = 8, units = "cm")

cat('check integration result B to A2B...\n')
cellIds.Pred <- classifyNew.SVM(ae.B.seurat@assays$RNA@scale.data[feature.idx.B,], Idents(ae.B.seurat), 
  as.data.frame(A2B.seurat@assays$RNA@scale.data[feature.idx.A2B,]), 
  kernel = 'linear', scale = F)
ae.B.vs.A2B.adj.rand  <- adj.rand.index(Idents(ae.B.seurat), cellIds.Pred$pop)
ae.B.vs.A2B.ji        <- cluster_similarity(Idents(ae.B.seurat), cellIds.Pred$pop)




#######################################
# check ajusted rand index of integrated data clustering vs. input data clusterings

cat("########## RESULT ##########\n")
cat("# Label transfer across domain with SVM\n")
cat("\nBEFORE integration: \nARI ae.A vs. ae.B: ", ae.A.vs.B.adj.rand, 
  "\nARI ae.B vs. ae.A: ",  ae.B.vs.A.adj.rand, 
  "\nJI ae.A vs. ae.B: ",   ae.A.vs.B.ji, 
  "\nJI ae.B vs. ae.A: ",   ae.B.vs.A.ji, "\n")
cat("\nAFTER integration: \nARI ae.A vs. ae.B2A: ", ae.A.vs.B2A.adj.rand, 
  "\nARI ae.B vs. ae.A2B: ",  ae.B.vs.A2B.adj.rand, 
  "\nJI ae.A vs. ae.B2A: ",   ae.A.vs.B2A.ji, 
  "\nJI ae.B vs. ae.A2B: ",   ae.B.vs.A2B.ji, "\n")

cat('Done...\n')
print(Sys.time())
