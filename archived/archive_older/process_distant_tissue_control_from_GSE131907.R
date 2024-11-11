gc()
rm(list = ls())
path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/scrna_gex/scrna_gex_pipeline"

maindir <- "/media/hieunguyen/HNSD_mini/data/UKK_Lung_integrated_datasets"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

MINCELLS <- 5
MINGENES <- 50
cluster.resolution <- 0.5
num.PC.used.in.UMAP <- 30
num.PC.used.in.Clustering <- 30
num.PCA <- 30
my_random_seed <- 42

PROJECT <- "GSE131907_normal"

if (file.exists(file.path(maindir, "integrated_dataset", sprintf("%s_distant_tissue_control.added_cell_type.rds", PROJECT))) == FALSE){
  if (file.exists(file.path(maindir, "integrated_dataset", sprintf("%s_distant_tissue_control.rds", PROJECT))) == FALSE){
    path.to.input <- file.path("/media/hieunguyen/HNSD_mini/data/UKK_Lung_integrated_datasets/integrated_dataset/GSE131907_Lung_Cancer_raw_UMI_matrix.controlOnly.rds")
    input.data <- readRDS(path.to.input)
    s.obj <- CreateSeuratObject(counts = input.data , 
                                min.cells = MINCELLS, 
                                min.features = MINGENES, 
                                project = PROJECT)
    
    # estimate the percentage of mapped reads to Mitochondrial and Ribosome genes
    s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, 
                                                  pattern = "^mt-|^MT-")
    s.obj[["percent.ribo"]] <- PercentageFeatureSet(s.obj, 
                                                    pattern = "^Rpl|^Rps|^RPL|^RPS")
    
    s.obj <- NormalizeData(s.obj) # ---> use Log Normalized
    s.obj <- FindVariableFeatures(s.obj, selection.method = "vst")
    s.obj <- ScaleData(s.obj, features = rownames(s.obj))
    
    pca_reduction_name <- "RNA_PCA"
    umap_reduction_name <- "RNA_UMAP"
    
    s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = FALSE, reduction.name=pca_reduction_name)
    s.obj <- RunUMAP(s.obj, reduction = pca_reduction_name, 
                     dims = 1:num.PC.used.in.UMAP, reduction.name=umap_reduction_name,
                     seed.use = my_random_seed, umap.method = "uwot")
    # clustering 
    s.obj <- FindNeighbors(s.obj, reduction = pca_reduction_name, dims = 1:num.PC.used.in.Clustering)
    s.obj <- FindClusters(s.obj, resolution = cluster.resolution, random.seed = 0)
    
    meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
      rowwise() %>%
      mutate(name = str_split(barcode, "_")[[1]][[3]]) %>%
      column_to_rownames("barcode")
    
    meta.data <- meta.data[row.names(s.obj@meta.data),]
    s.obj <- AddMetaData(object = s.obj, metadata = meta.data$name, col.name = "name")
    
    saveRDS(s.obj, file.path(maindir, "integrated_dataset", sprintf("%s_distant_tissue_control.rds", PROJECT)))
  } else {
    s.obj <- readRDS(file.path(maindir, "integrated_dataset", sprintf("%s_distant_tissue_control.rds", PROJECT)))
  }
  
  if (file.exists(file.path(maindir, "integrated_dataset", "integrated_datasets_Myeloid_only.rds")) == FALSE){
    path.to.input.s.obj.integrated <- file.path(maindir, "integrated_dataset", "refquery_final.rds")
    s.obj.integrated <- readRDS(path.to.input.s.obj.integrated)
    s.obj.myeloid <- subset(s.obj.integrated, Cell_Cluster_level1 == "Myeloid")
    saveRDS(s.obj.myeloid, file.path(maindir, "integrated_dataset", "integrated_datasets_Myeloid_only.rds"))  
  } else {
    s.obj.myeloid <- readRDS(file.path(maindir, "integrated_dataset", "integrated_datasets_Myeloid_only.rds"))
  }
  
  #####----------------------------------------------------------------------#####
  ##### Transfer cell annotation to the control dataset
  #####----------------------------------------------------------------------#####
  
  path.to.input.s.obj.integrated <- file.path(maindir, "integrated_dataset", "refquery_final.rds")
  s.obj.integrated <- readRDS(path.to.input.s.obj.integrated)
  DefaultAssay(s.obj.integrated) <- "RNA"
  
  s.obj.integrated <- FindVariableFeatures(s.obj.integrated)
  
  anchors <- FindTransferAnchors(reference = s.obj.integrated, query = s.obj, dims = 1:30)
  
  predictions1 <- TransferData(anchorset = anchors, refdata = s.obj.integrated$Cell_Cluster_level2, dims = 1:30)
  s.obj <- AddMetaData(s.obj, metadata = predictions1)
  
  DimPlot(object = s.obj, reduction = "RNA_UMAP", label = TRUE, label.box = TRUE, group.by = "predicted.id")
  
  ##### add true cell type from dataset annotation (GSE deposited data)
  lung.cell.annotation <- read.csv("/media/hieunguyen/HNSD_mini/data/UKK_Lung_integrated_datasets/integrated_dataset/GSE131907_Lung_Cancer_cell_annotation.txt", sep = "\t")
  
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
  
  meta.data <- merge(meta.data, subset(lung.cell.annotation, select = c(Index, Cell_type)), by.x = "barcode", by.y = "Index") %>%
    column_to_rownames("barcode")
  
  meta.data <- meta.data[row.names(s.obj@meta.data),]
  
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$Cell_type, col.name = "true_Cell_type")
  
  # DimPlot(object = s.obj, reduction = "RNA_UMAP", label = TRUE, label.box = TRUE, group.by = "true_Cell_type")
  
  saveRDS(s.obj, file.path(maindir, "integrated_dataset", sprintf("%s_distant_tissue_control.added_cell_type.rds", PROJECT)))
} else {
  print("Object existed! Reading in....")
  s.obj <- readRDS(file.path(maindir, "integrated_dataset", sprintf("%s_distant_tissue_control.added_cell_type.rds", PROJECT)))
}

##### merge
if (file.exists(file.path(maindir, "integrated_dataset", "merge_myloid_from_integrated_dataset_and_control_dataset.rds")) == FALSE){
  s.obj.myeloid <- readRDS(file.path(maindir, "integrated_dataset", "integrated_datasets_Myeloid_only.rds"))
  
  s.obj.merge <- merge(s.obj, s.obj.myeloid)
  s.obj.merge <- DietSeurat(s.obj.merge)
  
  remove(s.obj)
  remove(s.obj.myeloid)
  
  s.obj.merge <- NormalizeData(s.obj.merge) # ---> use Log Normalized
  s.obj.merge <- FindVariableFeatures(s.obj.merge, selection.method = "vst")
  s.obj.merge <- ScaleData(s.obj.merge)
  
  pca_reduction_name <- "RNA_PCA"
  umap_reduction_name <- "RNA_UMAP"
  
  s.obj.merge <- RunPCA(s.obj.merge, npcs = num.PCA, verbose = FALSE, reduction.name=pca_reduction_name)
  s.obj.merge <- RunUMAP(s.obj.merge, reduction = pca_reduction_name, 
                         dims = 1:num.PC.used.in.UMAP, reduction.name=umap_reduction_name,
                         seed.use = my_random_seed, umap.method = "uwot")
  # clustering 
  s.obj.merge <- FindNeighbors(s.obj.merge, reduction = pca_reduction_name, dims = 1:num.PC.used.in.Clustering)
  s.obj.merge <- FindClusters(s.obj.merge, resolution = cluster.resolution, random.seed = 0)
  
  meta.data <- s.obj.merge@meta.data %>% rownames_to_column("barcode") %>%
    mutate(sampleid = ifelse(is.na(name), "integrated", "control")) %>%
    mutate(final.celltype = ifelse(sampleid == "control", predicted.id, Cell_Cluster_level2)) %>%
    column_to_rownames("barcode")
  
  meta.data <- meta.data[row.names(s.obj.merge@meta.data),]
  s.obj.merge <- AddMetaData(object = s.obj.merge, metadata = meta.data$sampleid, col.name = "sampleid")
  s.obj.merge <- AddMetaData(object = s.obj.merge, metadata = meta.data$final.celltype, col.name = "final.celltype")
  
  saveRDS(s.obj.merge, file.path(maindir, "integrated_dataset", "merge_myloid_from_integrated_dataset_and_control_dataset.rds"))
  write.csv(data.frame(status = c("finished!")), file.path(maindir, "integrated_dataset", "finished.csv"))
} else {
  s.obj.merge <- readRDS(file.path(maindir, "integrated_dataset", "merge_myloid_from_integrated_dataset_and_control_dataset.rds"))
}

remove(s.obj)

##### preparation
path.to.save.output <- file.path(maindir, "integrated_dataset", "output")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

dir.create(file.path(path.to.save.output, "DE_results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.save.output, "violin_plot"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.save.output, "dot_plot"), showWarnings = FALSE, recursive = TRUE)

kinase.genelist <- readxl::read_excel("/media/hieunguyen/HNSD01/src/UKK_Lung_cancer_datasets/gene_list.xlsx", col_names = FALSE)$`...1`
check.celltypes <- c("Proliferating Mac", 
                     "Alveolar Mac", 
                     "cDC2/moDCs",
                     "Monocytes",
                     "Low quality Mac", 
                     "Neutrophils",
                     "Lipid-associated Mac")

s.obj.merge <- JoinLayers(s.obj.merge)

# install presto for faster test
# devtools::install_github('immunogenomics/presto')

##### Differential gene expression for each cell type between "integrated" vs "control"
# DimPlot(object = s.obj.merge, reduction = "RNA_UMAP", group.by = "final.celltype", label = TRUE, label.box = TRUE)
# DimPlot(object = s.obj.merge, reduction = "RNA_UMAP", group.by = "sampleid", label = TRUE, label.box = TRUE)

for (input.celltype in check.celltypes){
  print(sprintf("Working on %s", input.celltype))
  tmp.diff.markers <- FindMarkers(object = subset(s.obj.merge, final.celltype == input.celltype), 
                                  ident.1 = "integrated", ident.2 = "control", group.by = "sampleid", assay = "RNA", slot = "data", min.pct = 0.1)
  tmp.diff.markers <- tmp.diff.markers %>% 
    rownames_to_column("Gene") %>% 
    rowwise() %>%
    mutate(abs.log2FC = abs(avg_log2FC)) %>%
    subset(p_val_adj <= 0.05)  %>%
    arrange(desc(abs.log2FC))
    
  tmp.diff.markers <- subset(tmp.diff.markers, substr(tmp.diff.markers$Gene, start = 1, stop = 2) != "MT")
  tmp.diff.markers <- subset(tmp.diff.markers, substr(tmp.diff.markers$Gene, start = 1, stop = 3) != "HLA")
  
  writexl::write_xlsx(tmp.diff.markers, file.path(path.to.save.output, "DE_results", sprintf("DEA_in_%s.xlsx", str_replace(input.celltype, "/", "-"))))
}

all.s.obj <- hash()
for (input.celltype in check.celltypes){
  all.s.obj[[input.celltype]] <- subset(s.obj.merge, final.celltype == input.celltype)
}

##### plot
library("ggpubr")
kinase.genelist <- intersect(kinase.genelist, row.names(s.obj.merge))

DefaultAssay(s.obj.merge) <- "RNA"

for (input.celltype in names(all.s.obj)){
  print(sprintf("%s: %s", input.celltype, ncol(all.s.obj[[input.celltype]])))
  print(table(all.s.obj[[input.celltype]]$sampleid))
}

# for (input.celltype in check.celltypes){
#   print(sprintf("working on %s", input.celltype))
#   dir.create(file.path(path.to.save.output, "violin_plot", input.celltype), showWarnings = FALSE, recursive = TRUE)
#   for (gene.id in kinase.genelist){
#     print(sprintf("working on gene %s", gene.id))
#     p <- VlnPlot(object = all.s.obj[[input.celltype]], features = c(gene.id), group.by = "sampleid", pt.size = 1) + 
#       stat_compare_means(method = "wilcox")
#     ggsave(plot = p, filename = sprintf("violinplot_%s_%s.png", gene.id, str_replace(input.celltype, "/", "-")), path = file.path(path.to.save.output, "violin_plot", input.celltype), device = "png", width = 14, height = 10, dpi = 300)
#   }
# }

new.gene.list <- read.csv(file.path("/media/hieunguyen/HNSD_mini/data/UKK_Lung_integrated_datasets/new_gene_set_20240409.csv"))

for (input.celltype in check.celltypes){
  print(sprintf("working on %s", input.celltype))
  dir.create(file.path(path.to.save.output, "violin_plot_20240409", input.celltype), showWarnings = FALSE, recursive = TRUE)
  for (gene.id in new.gene.list$Gene){
    if (gene.id %in% row.names(all.s.obj[[input.celltype]])){
      print(sprintf("working on gene %s", gene.id))
      p <- VlnPlot(object = all.s.obj[[input.celltype]], features = c(gene.id), group.by = "sampleid", pt.size = 1) +
        stat_compare_means(method = "wilcox")
      ggsave(plot = p, filename = sprintf("violinplot_%s_%s.png", gene.id, str_replace(input.celltype, "/", "-")), path = file.path(path.to.save.output, "violin_plot_20240409", str_replace(input.celltype, "/", "-")), device = "png", width = 14, height = 10, dpi = 300)      
    }

  }
}
