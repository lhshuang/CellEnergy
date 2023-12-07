# CellEnergy

## 1. Introduction

![flow chart](https://github.com/lhshuang/CellEnergy/blob/main/flow%20chart.png)

CellEnergy, a novel approach that assigns an energy indicator to individual cells based on Gene Local Network Energy (GLNE). CellEnergy can measure cellular plasticity, identify finer cell subpopulations, automatically detect terminal cell states without prior biological knowledge, and support mapping cell fate and deciphering driver genes. This package has a dependency on R version (R >= 4.3.0). 

## 2. Installation

### 2.1 Dependencies

Most dependencies can be installed  through CRAN. Some may need to be installed manually through BioConductor (https://bioconductor.org/).

Install BiocManager and Bioconductor dependency packages
```r
install.packages("BiocManager")
install.packages("devtools")
BiocManager::install(c('mclust', 'philentropy', 'plot3D', 'reticulate', 'plyr', 'dplyr', 'ggplot2',"limma","statmod"))
```
### 2.2 Install
Install the package from GitHub.
```r
library("devtools")
devtools::install_github("lhshuang/CellEnergy")
```
Load the package
```r
library("CellEnergy")
```

## 3. Tutorials

 [mouse pancreatic endocrine lineage formation.](https://github.com/lhshuang/CellEnergy/tree/main/inst/Tutorials) Data available [here](https://github.com/lhshuang/CellEnergy/tree/main/inst/extdata).