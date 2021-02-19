Script for data analysis
================
João C. D. Muzzi, Mauro A. A. Castro <br>

19 February 2021.


Context
----
Despite progress in understanding the biology of adrenocortical carcinoma (ACC), the treatment options have not changed in the last 3 decades. Our goal was to improve the knowledge on immune pathways that could be druggable targets to enhance the immune response of patients. Our strategy was to revisit the differences based on a molecular classification on High and Low Steroid Phenotypes (HSP and LSP respectively) proposed by Zheng et al. (2016), using The Cancer Genome Atlas (TCGA) ACC database and public data sets published by reference studies of Thorsson et al., (2018) and Sturm et al. (2019) for TCGA samples. This script reproduces all results published in Muzzi et al. (2021) and serves as complementary material.
In Preprocessing.Rmd file, we show how to obtain and preprocess the data sets used in the study, at the end of the section the data sets are saved as RData and can be imported for further analysis. Here we use the preprocessed data and show the steps for data analysis and how to obtain the results and plots shown in Muzzi et al. (2021). 

```r
knitr::opts_chunk$set(fig.show="hide")
```

Package Installation and Loading
----
```r
# Packages Installation
## Find out if packages are already installed, if not proceed with installation
### R packages in CRAN repository
Rpackages <- c("tidyverse","scales","igraph", "ggpubr", "ggforce", "ggrepel",
               "cowplot", "survminer", "survival", "gridExtra")
Rpackages.new <- Rpackages[!(Rpackages %in% installed.packages()[,"Package"])]
if(length(Rpackages.new)) install.packages(Rpackages.new)
### R packages in BioConductor
Rpackages <- c("DESeq2","SummarizedExperiment","tmod", "qvalue","ComplexHeatmap",
               "RedeR")
Rpackages.new <- Rpackages[!(Rpackages %in% installed.packages()[,"Package"])]
if(length(Rpackages.new)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install(Rpackages.new)
}
### R packages from other sources
Rpackages <- c("devtools", "ggunchained")
Rpackages.new <- Rpackages[!(Rpackages %in% installed.packages()[,"Package"])]
if(length(Rpackages.new)){
  install.packages("devtools")
  remotes::install_github("JanCoUnchained/ggunchained")
}
```
