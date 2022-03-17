
library(optparse,  quietly = T)
library(Seurat,    quietly = T)
library(tidyverse, quietly = T)
library(cowplot,   quietly = T)
library(pdfCluster, quietly = T)
library(pheatmap,  quietly = T)
library(clusteval, quietly = T)
library(scclusteval,quietly = T)


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

A2B 	<-  (read.table(A2B.file, sep = ',', header = T,
                    row.names = 1, stringsAsFactors = F))

B2A   <-  read.table(B2A.file, sep = ',', header = T, 
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




jacard_vector_rna = c()
jacard_vector_rna_test = c()


jacard=PairWiseJaccardSetsHeatmap(ae.A.seurat@active.ident, A2B.seurat@active.ident, best_match = FALSE,
                                                                     title = NULL,
                                                                     col_low = "white",
                                                                     col_high = "red",
                                                                     cluster_rows = F,
                                                                     cluster_columns = F,
                                                                     show_row_dend = F,
                                                                    show_column_dend = F,) 
 jacard_distance = jacard@matrix
 is.na(jacard_distance) <- jacard_distance==0
 jacard_distance = jacard_distance[,colSums(is.na(jacard_distance))<nrow(jacard_distance)]
 jacard_vector_rna = append(jacard_vector_rna, mean(apply(jacard_distance,1,sum, na.rm=TRUE)))
 validation_split_cells = tail(names(ae.A.seurat@active.ident),floor((length(names(B2A.seurat@active.ident))*20)/100))
 A = ae.A.seurat@active.ident[names(ae.A.seurat@active.ident) %in% validation_split_cells]
 B = B2A.seurat@active.ident[names(A2B.seurat@active.ident) %in% validation_split_cells]
 jacard=PairWiseJaccardSetsHeatmap(A, B,show_row_dend = F, show_column_dend = F,cluster_row = F, cluster_column =F)
 jacard_distance = jacard@matrix
 is.na(jacard_distance) <- jacard_distance==0
 jacard_distance = jacard_distance[,colSums(is.na(jacard_distance))<nrow(jacard_distance)]                                                                    
 jacard_vector_rna_test = append(jacard_vector_rna_test, mean(apply(jacard_distance,1,sum, na.rm=TRUE)))                                                                   
 print(jacard_vector_rna_test)
