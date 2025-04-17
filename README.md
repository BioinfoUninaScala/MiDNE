<p align="right">
 <img src="https://github.com/BioinfoUninaScala/MiDNE/blob/main/MiDNE_logo.png" width="250" alt="MiDNE Logo">
</p>

# MiDNE: Multi-omics genes and Drugs Network Embedding
#### A novel R package for integrating gene-centered multi-omics data with drug information. 
MiDNE is a computational pipeline to predict condition-specific gene-gene, drug-gene and drug-drug associations by integrating multi-omics data and drug information. MiDNE leverages network-based approach to model the multi-omics relationships between genes and drugs into a heterogeneous network. The neighborhood of each node is explored and numerically encoded through a Random Walk with Restart procedure. MiDNE can also learn a low-dimensional representation of each node in the integrated network, and then facilitates its visualization and interpretation via clustaering and enrichment analyses. 

### Installation 
In R console, run 

```r
library(devtools)
install_github("BioinfoUninaScala/MiDNE", 
               build_vignettes=FALSE, 
               repos=BiocManager::repositories(),
               dependencies=TRUE, type="source")
```
----------
