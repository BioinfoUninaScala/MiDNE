<p align="center">
 <img src="https://github.com/BioinfoUninaScala/MiDNE/blob/main/MiDNE_logo.png" width="350" alt="MiDNE Logo">
</p>

# MiDNE: Multi-omics genes and Drugs Network Embedding
#### A novel R package for integrating gene-centered multi-omics data with drug information. 
MiDNE is a computational pipeline to predict condition-specific gene-gene, drug-gene and drug-drug associations by integrating multi-omics data and drug information. MiDNE leverages network-based approach to model the multi-omics relationships between genes and drugs into a heterogeneous network. The neighborhood of each node is explored and numerically encoded through a Random Walk with Restart procedure. MiDNE can also learn a low-dimensional representation of each node in the integrated network, and then facilitates its visualization and interpretation via clustering and enrichment analyses. 

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

### Code and data
#### `R/` directory
##### Load example files and networks
- `loadAnnotation.R`:
- `loadOmicsMat.R`: 
- `loadDrugs.R`:     
- `loadOmicsNet.R`:    
- `loadDrugNet.R`:       
- `loadBipartiteNet.R`:  

##### Pre-processing omics matrices
- `remove_zero_rows.R`:
- `normalize_omics.R`:
- `get_intersection_matrices.R`:     

##### Network inference
- `gen_coCNVnet.R`:           
- `gen_coExpressionNet.R`:    
- `gen_coAbundanceNet.R`:      
- `gen_coDNAmethNet.R`:       
- `gen_isolatedDrugNet.R`:     
- `fisher_test_post_hoc.R`:   
- `get_filtered_corMat_by_adj_pval.R`:  
- `omics_network_inference.R`:

##### Integration 
- `create_multiplex.R`:    
- `prune_multiplex_network.R`:
- `create_layer_transition_matrix.R`:
- `gen_sim_mat_M.R`:                  
- `gen_sim_mat_MH.R`:             
- `get_embedding.R`:             
- `get_parallel_umap_embedding.R`:      

##### Plotting
- `plot_2D_matrix.R`:

##### Interactive interface
- `MiDNEshiny.R`:

##### Other
- `utils.R`:


#### `inst/extdata/` directory
##### /data
- `/pharmacological/FDAdrugs.RDS`:
- `/pharmacological/ALLdrugs.RDS`:
- `/pharmacological/FDAdrugs.RDS`:
- `/pharmacological/ALLdrugs.RDS`:
- `/biological/BRCA_Methylation_Meth450.RDS`:
- `/biological/BRCA_proteome_CDAP.RDS`:
- `/biological/BRCA_SCNA.RDS`:

##### /annotation
- `Human__TCGA_BRCA__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi`:
- `all_genes_drugs_annotation.RDS`:

##### /networks
- `/biological/BRCA_filt_codel_network.RDS`:
- `/biological/BRCA_filt_coamp_network.RDS`:
- `/biological/BRCA_filt_cometh_network.RDS`:
- `/biological/BRCA_filt_prot_net.RDS`:
- `/biological/BRCA_filt_coexpr_network.RDS`:
- `/pharmacological/FDAdrugs_net.csv`:
- `/bipartite/FDA_active_DRUGBANK_bnet.RDS`:

##### /similarity_matrices
- `emb_FDA_active_drug_5omics_uRWRMHmat.RDS`:
- `umap_emb_FDA_active_drug_5omics_uRWRMHmat.RDS`:

