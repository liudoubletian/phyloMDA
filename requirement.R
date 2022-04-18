#####R libraries#####
#CRAN library
install.packages(MGLM);library(MGLM)
install.packages(plyr);library(plyr)
install.packages(caper);library(caper)
install.packages(genlasso);library(genlasso)
install.packages(magrittr);library(magrittr)
install.packages(foreach);library(foreach)
install.packages(ape);library(ape)
install.packages(miLineage);library(miLineage)
install.packages(ggplot2);library(ggplot2)
install.packages(dplyr);library(dplyr)
install.packages(readxl);library(readxl)
library(methods)
#Bioconductor library
install.packages("BiocManager")
BiocManager::install("phyloseq");library(phyloseq)
BiocManager::install("ggtree");library(ggtree)
#other libraries
devtools::install_github("ZRChao/adaANCOM");library(adaANCOM)
