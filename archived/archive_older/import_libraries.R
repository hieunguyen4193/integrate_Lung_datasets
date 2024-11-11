################################################################################
# This script is used to clean the envrionment and import all necessary packages
################################################################################
# Specify the list of packages that need to be imported ########################
list.of.packages <- c("Seurat",
                      "SingleCellExperiment",
                      "optparse", 
                      "comprehenr", 
                      "tidyverse", 
                      "ggplot2", 
                      "SoupX",
                      "comprehenr",
                      "DoubletFinder",
                      "vroom",
                      "hash",
                      "DT",
                      "janitor",
                      "knitr",
                      "circlize",
                      "formattable",
                      "htmlwidgets",
                      "plotly",
                      "stringr",
                      "rcompanion",
                      "argparse",
                      "scatterpie", 
                      "scales"
)

## temp removed: rstatix

bioc.packages <- c("celda", 
                   "BiocSingular", 
                   "PCAtools", 
                   "SingleCellExperiment",
                   "sctransform", 
                   "progeny")

# Check if packages are installed ##############################################

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.bioc.packages <- bioc.packages[!(bioc.packages %in% installed.packages()[,"Package"])]

# Install new packages #########################################################

print("Installing new packages by install.packages...")
if(length(new.packages)) install.packages(new.packages)
# if(length(new.packages)) install.packages(new.packages, repos = "http://cran.rstudio.com/")
print("Finish installing new packages by install.packages")
print("Start installing packages from Bioconductor...")

if(length(new.bioc.packages)) BiocManager::install(new.bioc.packages, update = FALSE, ask = TRUE)
print("Finish installing packages from Bioconductor")


# Import all packages ##########################################################
package_loading_Status <- lapply(list.of.packages, 
                                 require, 
                                 character.only = TRUE)

package_loading_Status_bioc <- lapply(bioc.packages, 
                                      require, 
                                      character.only = TRUE)


# EOF ##########################################################################
