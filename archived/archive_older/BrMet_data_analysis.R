gc()
rm(list = ls())


control.groups <- list(
  `Klemm 2020` = c("non tumour MG"),
  `Wischnewski 2023` = c("non tumour MG", "non tumour MDM")
)

cancer.groups <- list(
  `Klemm 2020` = c("Glioma MG", 
                   "Glioma MDM",
                   "BrM MG",
                   "BrM MDM"),
  `Wischnewski 2023` = c("Glioma MG",
                         "Glioma MDM",
                         "BrM MG",
                         "BrM MDM")
)

data.version <- "20240406" 

maindir <- file.path("/media/hieunguyen/HNSD_mini/data/UKK_Lung_integrated_datasets/BrainMet", data.version)
path.to.metadata <- file.path(maindir, "Brain_Met_dataset_metadata.xlsx")

meta.data <- readxl::read_excel(path.to.metadata, sheet = "metadata")
org.metadata <- meta.data

tk.genelist <- readxl::read_excel(path.to.metadata, sheet = "TK_gene_list")
colnames(tk.genelist) <- c("Gene")
new.gene.list <- read.csv(file.path("/media/hieunguyen/HNSD_mini/data/UKK_Lung_integrated_datasets/new_gene_set_20240409.csv"))
colnames(new.gene.list) <- c("Gene")

all.gene.lists <- list(TK_genes = tk.genelist,
                       new_gene_list = new.gene.list)

path.to.input.data <- file.path(maindir, "BrM dataset.xlsx")
path.to.main.src <- "/media/hieunguyen/HNSD01/src/UKK_Lung_cancer_datasets"

source(file.path(path.to.main.src, "import_libraries.R"))
source(file.path(path.to.main.src, "helper_functions.R"))
library(DESeq2)
library(ggpubr)

count.thres <- 10

#####----------------------------------------------------------------------#####
##### HELPER FUNCTIONS
#####----------------------------------------------------------------------#####
run_deseq_and_plot_selected_genes <- function(data.source, selected.control, selected.cancer, genelist){
  path.to.save.output <- file.path(maindir, "output_20240411", data.source, genelist, sprintf("%s_vs_%s", selected.cancer, selected.control))
  dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path.to.save.output, "figures"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path.to.save.output, "summary_table"), showWarnings = FALSE, recursive = TRUE)
  
  maindf <- readxl::read_excel(file.path(path.to.input.data), sheet = data.source)
  meta.data <- subset(org.metadata, org.metadata$Source == data.source)
  meta.data <- subset(meta.data, meta.data$Label %in% c(selected.control, selected.cancer))
  
  maindf <- maindf %>% column_to_rownames("...1")
  
  meta.data <- meta.data %>% rowwise() %>%
    mutate(condition = ifelse(Label == selected.cancer, "positive", "control"))
  meta.data$condition <- factor(meta.data$condition, levels = c("control", "positive"))
  
  maindf <- maindf[, meta.data$SampleID]
  if (data.source == "Klemm 2020"){
    dds <- DESeqDataSetFromMatrix(countData = maindf %>% as.matrix(),
                                  colData = meta.data,
                                  design = ~ condition)
  } else if (data.source == "Wischnewski 2023"){
    dds <- DESeqDataSetFromMatrix(countData = floor(maindf) %>% as.matrix(),
                                  colData = meta.data,
                                  design = ~ condition)
  }
  
  dds <- dds[rowSums(counts(dds)) > count.thres, ]
  
  dds <- estimateSizeFactors(dds)
  
  vst_mat <- vst(dds, blind = TRUE)
  count.matrix <- assay(vst_mat)
  
  resdf <- data.frame()
  
  for (gene.id in all.gene.lists[[genelist]]$Gene){
    if (gene.id %in% row.names(dds)){
      # tmp.countdf <- norm.counts[gene.id, ] %>% as.data.frame() %>% rownames_to_column("SampleID")
      tmp.countdf <- count.matrix[gene.id, ] %>% as.data.frame() %>% rownames_to_column("SampleID")
      
      colnames(tmp.countdf) <- c("SampleID", gene.id)
      tmp.countdf <- merge(tmp.countdf, meta.data, by.x = "SampleID", by.y = "SampleID")
      max.val <- tmp.countdf[[gene.id]] %>% max()
      
      count.sample <- table(tmp.countdf$Label)
      
      test.res <- wilcox.test(subset(tmp.countdf, tmp.countdf$Label == selected.cancer)[[gene.id]], 
                              subset(tmp.countdf, tmp.countdf$Label == selected.control)[[gene.id]], 
                              paired = FALSE, alternative = "two.sided")
      
      logFC <- log2(mean(subset(tmp.countdf, tmp.countdf$Label == selected.cancer)[[gene.id]])/mean(subset(tmp.countdf, tmp.countdf$Label == selected.control)[[gene.id]]))
      
      title <- sprintf(" Gene: %s, wilcox.p-value = %s, log2FC = %s, \n %s: %s, %s: %s", gene.id, test.res$p.value %>% round(6), round(logFC, 6),
                       names(count.sample)[[1]],
                       count.sample[[1]],
                       names(count.sample)[[2]], count.sample[[2]])
      violin.plot <- tmp.countdf %>% ggplot(aes_string( x = "Label", y = gene.id, fill = "Label")) + geom_boxplot() + 
        theme_pubr(legend = "none") + 
        geom_jitter(color="black", size=2, alpha=0.9, height = 0, width = 0.1) + 
        # stat_summary(fun = mean, geom='crossbar', colour = "black", position = position_dodge()) + 
        ggtitle(title)
      
      ggsave(plot = violin.plot, path = file.path(path.to.save.output, "figures"), filename = sprintf("Gene_%s_in_%s_vs_%s.png", gene.id, selected.cancer, selected.control), device = "png", dpi = 300, width = 14, height = 10)
      
      tmp.resdf <- data.frame(gene = c(gene.id), pval = c(test.res$p.value), log2FC = c(logFC))
      resdf <- rbind(resdf, tmp.resdf)
    }
  }
  
  resdf <- resdf %>% rowwise() %>%
    mutate(significant_diff = ifelse(pval <= 0.05, "yes", "no"))
  
  writexl::write_xlsx(resdf, file.path(path.to.save.output, "summary_table", sprintf("All_selected_genes_in_%s_vs_%s.table.xlsx", selected.cancer, selected.control)))
}

#####----------------------------------------------------------------------#####
##### MAIN RUN
#####----------------------------------------------------------------------#####

# for (data.source in c("Klemm 2020", "Wischnewski 2023")){
#   for (selected.control in control.groups[[data.source]]) {
#     for (selected.cancer in cancer.groups[[data.source]]){
#       for (genelist in names(all.gene.lists)){
#         print(sprintf("Working on %s, %s, %s, %s", data.source, selected.control, selected.cancer, genelist))
#         run_deseq_and_plot_selected_genes(data.source, selected.control, selected.cancer, genelist)
#       }
#     }
#   }  
# }  

for (data.source in c("Klemm 2020", "Wischnewski 2023")){
  for (genelist in names(all.gene.lists)){
    run_deseq_and_plot_selected_genes(data.source, "BrM MDM", "BrM MG", genelist)
  }
}

for (data.source in c("Klemm 2020", "Wischnewski 2023")){
  for (genelist in names(all.gene.lists)){
    run_deseq_and_plot_selected_genes(data.source, "Glioma MDM", "Glioma MG", genelist)
  }
}
