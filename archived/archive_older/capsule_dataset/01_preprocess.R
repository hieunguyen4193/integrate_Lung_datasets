gc()
rm(list = ls())

set.seed(411)

library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)

# BiocManager::install("glmGamPoi", update = FALSE)

nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 30
pHB_lower <- 0
pHB_upper <- 5

theme_set(theme_cowplot())

#color scheme
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E",
  p028 = "#E1BD6D",
  p029 = "#35274A",
  p030 = "#F2300F",
  p031 = "#7294D4",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467")


# Data loading and QC


### sample list
maindir <- "/media/hieunguyen/HNSD_mini/data/UKK_Lung_integrated_datasets/capsule-8321305"
path.to.main.output <- file.path(maindir, "output")
path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

if (file.exists(file.path(path.to.01.output, "seurat_object.preprocessed.rds")) == FALSE){
  
  samples <- read_excel(file.path(maindir, "data/metadata/patients_metadata.xlsx"), range = cell_cols("A:A")) %>% .$sample_id
  
  ### import cellranger files from different data sets
  for (i in seq_along(samples)){
    print(sprintf("working on sample %s", i))
    assign(paste0("scs_data", i), Read10X(data.dir = paste0(file.path(maindir, "/data/cellranger/"), samples[i], "/filtered_feature_bc_matrix")))
  }
  
  ### create seurat objects from cellranger files
  all.s.obj <- list()
  for (i in seq_along(samples)){
    print(sprintf("working on sample %s", i))
    obj.name <- sprintf("seu_obj%s", i)
    all.s.obj[[obj.name]] <- CreateSeuratObject(counts = eval(parse(text = paste0("scs_data", i))), project = samples[i], min.cells = 3)
  }
  
  ### merge data sets
  seu_obj <- merge(x = all.s.obj[[names(all.s.obj)[[1]]]], 
                   y = all.s.obj[names(all.s.obj)[2:length(all.s.obj)]], add.cell.ids = samples, project = "lung")
  
  ### calculate mitochondrial, hemoglobin and ribosomal gene counts
  seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^MT-", col.name = "pMT")
  seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^HBA|^HBB", col.name = "pHB")
  seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^RPS|^RPL", col.name = "pRP")
  
  # Data Filtering 
  
  seu_obj_unfiltered <- seu_obj
  
  #saveRDS(seu_obj_unfiltered, file = "seurat_objects/all_unfiltered.RDS")
  
  ## After filtering
  
  seu_obj <- subset(seu_obj_unfiltered, subset = 
                      nFeature_RNA > nFeature_lower & 
                      nFeature_RNA < nFeature_upper & 
                      nCount_RNA > nCount_lower & 
                      nCount_RNA < nCount_upper & 
                      pMT < pMT_upper & 
                      pHB < pHB_upper)
  
  seu_obj <- SCTransform(seu_obj, verbose = T, vars.to.regress = c("nCount_RNA", "pMT"), conserve.memory = T)
  
  seu_obj <- RunPCA(seu_obj)
  
  seu_obj <- RunUMAP(seu_obj, dims = 1:15, verbose = T)
  seu_obj <- FindNeighbors(seu_obj, dims = 1:15)
  
  for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
    seu_obj <- FindClusters(seu_obj, resolution = i)
    print(DimPlot(seu_obj, reduction = "umap") + labs(title = paste0("resolution: ", i)))
  }
  
  # Cell cycle scoring
  
  ### add cell cycle, cc.genes loaded with Seurat
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  score_cc <- function(seu_obj) {
    seu_obj <- CellCycleScoring(seu_obj, s.genes, g2m.genes)
    seu_obj@meta.data$CC.Diff <- seu_obj@meta.data$S.Score - seu_obj@meta.data$G2M.Score
    return(seu_obj)
  }
  
  seu_obj <- score_cc(seu_obj)
  
  saveRDS(seu_obj, file.path(path.to.01.output, "seurat_object.preprocessed.rds"))
} else {
  seu_obj <- readRDS(file.path(path.to.01.output, "seurat_object.preprocessed.rds"))
}

seu_obj <- PrepSCTFindMarkers(object = seu_obj)

##### Find all cluster markers
Idents(seu_obj) <- "SCT_snn_res.0.2"
cluster.markers <- FindAllMarkers(object = seu_obj, assay = "SCT", test.use = "wilcox")

saveRDS(cluster.markers, file.path(path.to.01.output, "cluster_markers.rds"))
