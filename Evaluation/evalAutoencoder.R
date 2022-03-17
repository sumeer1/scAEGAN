#!/usr/bin/env Rscript
# autoencoder integration results

########################################
# example
# Rscript evalAutoencoder.R -r data/Autoencoder/RNA.csv.gz -a data/Autoencoder/ATAC.csv.gz -i data/Autoencoder/RNA_ATAC_Integrated.csv.gz -o tmp





ae.rnaseq.file      <- ''
ae.atac.file        <- ''
ae.integrated.file  <- ''
outDir              <- ''





cat('Loading data...\n')
print(Sys.time())

ae.rnaseq	<- t(read.table(ae.rnaseq.file, sep = ',', header = T, 
                          row.names = 1, stringsAsFactors = F))

ae.atac 	<- (read.table(ae.atac.file, sep = ',', header = T, 
                         row.names = 1, stringsAsFactors = F))

ae.integrated 	<-  (read.table(ae.integrated.file, sep = ',', header = T, 
                               row.names = 1, stringsAsFactors = F))

#######################################
# do seurat processing for plotting etc
cat('Making Seurat Objects...\n')
print(Sys.time())

# RNA-seq dataset
ae.rnaseq.seurat <- CreateSeuratObject(
  counts = ae.rnaseq, project = "ae.rnaseq",
  assay = "RNA", min.cells = 0, min.features = 0, names.field = 1)

ae.rnaseq.seurat <- ae.rnaseq.seurat %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dims = 1:20)

# ATAC dataset
ae.atac.seurat <- CreateSeuratObject(
  counts = ae.atac, project = "ae.atac",
  assay = "RNA", min.cells = 0, min.features = 0, names.field = 1)
ae.atac.seurat <- ae.atac.seurat %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:20) %>% #50 for CP
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dims = 1:20)

# integrated
ae.integrated.seurat <- CreateSeuratObject(
  counts = ae.integrated, project = "ae.integrated",
  assay = "RNA", min.cells = 0, min.features = 0, names.field = 1)
ae.integrated.seurat <- ae.integrated.seurat %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(npcs=20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.6, algorithm = 1) %>% 
  RunUMAP(dims = 1:20)

# add rnaseq and atac-derived clusterings to integrated dataset
ae.integrated.seurat <- ae.integrated.seurat %>% 
  AddMetaData(Idents(ae.rnaseq.seurat), col.name = 'RNAseqID') %>%
  AddMetaData(Idents(ae.atac.seurat), col.name = 'AtacID')


#######################################
# plotting all
plts <- list()
plts[["Old"]] 		<- DimPlot(ae.rnaseq.seurat, reduction = "umap")
plts[["Young"]] 			<- DimPlot(ae.atac.seurat, reduction = "umap")
plts[["Integrated"]] 	<- DimPlot(ae.integrated.seurat, reduction = "umap")
plts[["IntOldIDs"]] 	<- DimPlot(ae.integrated.seurat, reduction = "umap", group.by = 'RNAseqID')
plts[["IntYoungIDs"]] 	<- DimPlot(ae.integrated.seurat, reduction = "umap", group.by = 'AtacID')

p.vae.grid  <- plot_grid(plotlist = plts, nrow = 2, ncol = 3, labels = names(plts))
ggsave(file.path(outDir, paste0('AE_input_vs_reconst.pdf')), p.vae.grid,
       width = 30, height = 16, units = "cm")


# plot heatmaps of cluster similarity
scale.type <- 'row'
ae.int.vs.RNAseq.tab 	 <- table(Idents(ae.integrated.seurat), ae.integrated.seurat@meta.data$RNAseqID)
ae.int.vs.RNAseq.heatmap <- pheatmap(ae.int.vs.RNAseq.tab, cluster_rows = F, cluster_cols = F, 
                                     scale = scale.type, display_numbers = T, fontsize_number = 7, 
                                     main = 'Old Vs AE integrated')

ae.int.vs.atac.tab 	 	 <- table(Idents(ae.integrated.seurat),
                                ae.integrated.seurat@meta.data$AtacID)
ae.int.vs.atac.heatmap 	 <- pheatmap(ae.int.vs.atac.tab, cluster_rows = F, cluster_cols = F, 
                                     scale = scale.type, display_numbers = T, fontsize_number = 7, 
                                     main = 'Young Vs AE integrated')

p.vae.grid  <- plot_grid(plotlist = list(ae.int.vs.RNAseq.heatmap[[4]], ae.int.vs.atac.heatmap[[4]]), 
                         nrow = 1, ncol = 2, label_size = 10)

#######################################
# check numeric quality measures
# check ajusted rand index of integrated data clustering vs. input data clusterings
ae.int.vs.RNAseq.ari  <- adj.rand.index(Idents(ae.integrated.seurat), 
                                        ae.integrated.seurat@meta.data$RNAseqID)

ae.int.vs.atac.ari    <- adj.rand.index(Idents(ae.integrated.seurat), 
                                        ae.integrated.seurat@meta.data$AtacID)

ae.int.vs.RNAseq.ji <- cluster_similarity(Idents(ae.integrated.seurat), 
                                          ae.integrated.seurat@meta.data$RNAseqID, 
                                          similarity = "jaccard")

ae.int.vs.atac.ji <-cluster_similarity(Idents(ae.integrated.seurat), 
                                       ae.integrated.seurat@meta.data$AtacID, 
                                       similarity = "jaccard")

cat("########## RESULT ##########\n")
cat("ARI RNA-seq to integrated: ",  ae.int.vs.RNAseq.ari, 
    "\nARI ATAC-seq to integrated: ", ae.int.vs.atac.ari, 
    "\nJI RNA-seq to integrated: ",   ae.int.vs.RNAseq.ji, 
    "\nJI ATAC-seq to integrated: ",  ae.int.vs.atac.ji)
cat("\n############################\n")


cat('Done...\n')
print(Sys.time())
