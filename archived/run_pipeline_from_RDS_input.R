gc()
rm(list = ls())

if ("argparse" %in% installed.packages() == FALSE){
  install.packages("argparse")
}

for (sample.id in c("GSE148071",
                    "GSE153935",
                    "KU_loom")){
  #library(argparse)
  # parser <- ArgumentParser()
  # parser$add_argument("-i", "--sampleid", action="store",
  #                     help="Name of the input sample")
  # args <- parser$parse_args()
  # sample.id <- args$sampleid
  # sample.id <- "GSE119911"
  print("#####-----------------------------------------------------------#####")
  print(sprintf("WORKING ON SAMPLE %s", sample.id))
  print("#####-----------------------------------------------------------#####")
  
  my_random_seed <- 42
  set.seed(my_random_seed)
  
  # __________VDJ DATA ANYLYSIS PIPELINE__________
  PROJECT <- "Lung_integrated_dataset"
  version.name <- sprintf("SeuratV5_%s", sample.id)
  
  PROJECT.with.version <- sprintf("%s_%s", PROJECT, version.name)
  
  outdir <- "/media/hieunguyen/HNSD_mini/data/outdir"
  
  path.to.main.input <- "/media/hieunguyen/HNHD01/data/UKK/Lung/NSCLC_integrated/splitted_samples"
  
  path.to.main.output <- file.path(outdir, PROJECT)
  dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)
  
  path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
  path2src <- file.path(path.to.pipeline.src, "processes_src")
  
  source(file.path(path2src, "import_libraries.R"))
  source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))
  analysis.round <- "1st"
  
  path2input <- file.path(path.to.main.input, sample.id, sample.id)
  
  stage_lst <- list(sample.id)
  names(stage_lst) <- sample.id
  
  MINCELLS  <- 50
  MINGENES  <- 5
  
  save.RDS <- list(s1 = TRUE,
                   s2 = TRUE,
                   s3 = TRUE,
                   s4 = TRUE,
                   s5 = TRUE,
                   s6 = TRUE,
                   s7 = TRUE,
                   s8 = TRUE,
                   s8a = TRUE,
                   s9 = FALSE)
  if (sample.id == "GSE148071"){
    sw <- list(s1 = "on",
               s2 = "on",
               s3 = "on",
               s4 = "on",
               s5 = "on",
               s6 = "off",
               s7 = "on",
               s8 = "off",
               s8a = "on",
               s9 = "off")
  } else {
    sw <- list(s1 = "on",
               s2 = "on",
               s3 = "on",
               s4 = "on",
               s5 = "on",
               s6 = "on",
               s7 = "on",
               s8 = "off",
               s8a = "on",
               s9 = "off")
  }
  
  
  rerun <- list(s1 = FALSE, 
                s2 = FALSE,
                s3 = FALSE,
                s4 = FALSE,
                s5 = FALSE,
                s6 = FALSE,
                s7 = FALSE,
                s8 = FALSE,
                s8a = FALSE,
                s9 = FALSE)
  
  filter.thresholds <- list(nFeatureRNAfloor = NULL,
                            nFeatureRNAceiling = NULL,
                            nCountRNAfloor = NULL, 
                            nCountRNAceiling = NULL,
                            pct_mitofloor = NULL, 
                            pct_mitoceiling = 5,
                            pct_ribofloor = NULL, 
                            pct_riboceiling = NULL,
                            ambientRNA_thres = 0.5,
                            log10GenesPerUMI_thres = NULL)
  
  remove_doublet <- FALSE
  path.to.10X.doublet.estimation <- file.path("/media/hieunguyen/HNSD_mini/data/resources", "DoubletEstimation10X.csv")  
  #####--------------------------------------------------------------------#####
  ##### IMPORTANT INPUT PARAMS
  #####--------------------------------------------------------------------#####
  num.PCA <- 25
  num.PC.used.in.UMAP <- 25
  num.PC.used.in.Clustering <- 25
  regressOut_mode <- NULL
  features_to_regressOut <- NULL
  use.sctransform <- TRUE
  vars.to.regress <- c("percent.mt")
  
  print(sprintf("Working on sample.id %s", sample.id))
  path.to.output <- file.path(path.to.main.output, sprintf("%s_round", analysis.round), sample.id)
  dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)
  
  path.to.anno.contigs <- NULL
  path.to.count.clonaltype <- NULL
  filtered.barcodes <- NULL
  path.to.s3a <- NULL
  input.method <- "from_RDS"
  s.obj <- run_pipeline_GEX(path2src=path2src,
                            path2input=path2input,
                            path.to.logfile.dir=file.path(path.to.output, "logs"),
                            stage_lst=stage_lst,
                            path.to.10X.doublet.estimation=path.to.10X.doublet.estimation,
                            MINCELLS=MINCELLS,
                            MINGENES=MINGENES,
                            PROJECT=PROJECT.with.version,
                            remove_doublet=remove_doublet,
                            save.RDS=save.RDS,
                            path.to.output=path.to.output,
                            rerun=rerun,
                            DE.test="wilcox",
                            num.PCA=num.PCA,
                            num.PC.used.in.UMAP=num.PC.used.in.UMAP,
                            num.PC.used.in.Clustering=num.PC.used.in.Clustering,
                            use.sctransform=use.sctransform,
                            filtered.barcodes=filtered.barcodes,
                            filter.thresholds=filter.thresholds,
                            path.to.anno.contigs=path.to.anno.contigs,
                            path.to.count.clonaltype=path.to.count.clonaltype,
                            input.method = input.method,
                            my_random_seed = my_random_seed,
                            with.VDJ = TRUE, 
                            path.to.s3a.source = path.to.s3a, 
                            regressOut_mode = regressOut_mode,
                            features_to_regressOut = features_to_regressOut,
                            sw = sw,
                            vars.to.regress = vars.to.regress)
  
  #### ALWAYS REMEMBER TO SAVE SESSIONINFO !!!!!!
  writeLines(capture.output(sessionInfo()), file.path(path.to.output, sprintf("%s_sessionInfo.txt", PROJECT)))
  
  save.outdir <- "/media/hieunguyen/HNHD01/outdir"
  dir.create(file.path(save.outdir, PROJECT.with.version), showWarnings = FALSE, recursive = TRUE)
  system(sprintf("mv %s %s", path.to.output, file.path(save.outdir, PROJECT.with.version)))
  
}
