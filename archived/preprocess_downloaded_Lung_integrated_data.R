gc()
rm(list = ls())

library(dplyr)
library(tidyr)
library(tidyverse)

path.to.maindir <- "/media/hieunguyen/HNHD01/data/UKK/Lung/NSCLC_integrated"
dir.create(file.path(path.to.maindir, "splitted_samples"), showWarnings = FALSE, recursive = TRUE)
input.data <- readRDS(file.path(path.to.maindir, "RNA_rawcounts_matrix.rds"))
meta.data <- read.csv(file.path(path.to.maindir, "metadata.csv"))
sumdf <- data.frame()
for (study.id in unique(meta.data$Study)){
    dir.create(file.path(path.to.maindir, "splitted_samples", study.id, study.id), showWarnings = FALSE, recursive = TRUE)
    print(sprintf("Working on %s", study.id))
    selected.cells <- subset(meta.data, meta.data$Study == study.id)$X
    tmpdf <- data.frame(study = c(study.id), num.cell = c(length(selected.cells)))
    sumdf <- rbind(sumdf, tmpdf)
    output.obj <- input.data[, selected.cells]
    saveRDS(output.obj, file.path(path.to.maindir, "splitted_samples", study.id, study.id, sprintf("%s.rds", study.id)))
  }

