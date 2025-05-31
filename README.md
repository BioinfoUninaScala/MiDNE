<p align="center">
 <img src="https://github.com/BioinfoUninaScala/MiDNE/blob/main/MiDNE_logo.png" width="350" alt="MiDNE Logo">
</p>


# MiDNE: Multi-omics genes and Drugs Network Embedding
#### A novel R package for integrating gene-centered multi-omics data with drug information. 
MiDNE is a computational pipeline to predict condition-specific gene-gene, drug-gene and drug-drug associations by integrating multi-omics data and drug information. MiDNE leverages network-based approach to model the multi-omics relationships between genes and drugs into a heterogeneous network. The neighborhood of each node is explored and numerically encoded through a Random Walk with Restart procedure. MiDNE can also learn a low-dimensional representation of each node in the integrated network, and then facilitates its visualization and interpretation via clustering and enrichment analyses. 

---

## Installation 
In R console, run 

```r
library(devtools)
install_github("BioinfoUninaScala/MiDNE", 
               build_vignettes=FALSE, 
               dependencies=TRUE,
               type="source")
```
----------

### MiDNE shiny

The MiDNE pipeline is also implemented as a user-friendly shiny app, available both in MiDNE package and as docker image (https://hub.docker.com/r/bioinfouninascala/midne). <br/>

#### MiDNE function
In R, you can run:
```r
library(MiDNE)
MiDNEshiny()
```

#### MiDNE docker image
To pull the docker image, run in the terminal:
```
$ docker pull bioinfouninascala/midne
```

To run the docker image, map http port 3838 on the host port 8585:
```
$ docker run --rm -p 8585:3838 bioinfouninascala/midne
```

Finally, use MiDNE shiny in a browser by typing http://localhost:8585 (if you are on a local machine) or https://server_name:8585 (if you are on a server machine).
<br/>

----------

## Code and data
### 📂 `R/` — Code Directory
**0. Load example files and networks**

MiDNE includes functions to load example datasets used as a case study, based on TCGA-BRCA omics data and drug information from DrugBank.
- `loadOmicsMat`: loads example omics matrices (e.g., gene expression, CNV, methylation) used in the MiDNE case study.
- `loadDrugs`: loads a single-column data frame cointaning drug IDs obtained from the DrugBank database.
- `loadOmicsNet`: loads example omics-based networks inferred by MiDNE.  
- `loadDrugNet`: loads an isolated drug network constructed by MiDNE.
- `loadBipartiteNet`: loads bipartite networks connecting drugs to genes based on their known physical interactions retrieved from the DrugBank database.
- `loadAnnotation`: loads annotation data depending on the selected mode. If "sample" is selected, it returns sample-level metadata from breast cancer datasets used as case study. This can be useful for users who wish to build gene-centered biological networks starting from specific subtypes. If "gene-drug" is selected, it loads an annotation matrix of known genes and drugs.
  
**1. Pre-processing omics matrices**
- `remove_zero_rows`: removes rows from a matrix - assumed to be features x samples - where all the elements are equal to zero.
- `normalize_omics`: applies normalization to omics datasets to make them comparable across samples.
- `get_intersection_matrices`: takes a list of matrices and returns a new list based on their shared row names.

**2. Network inference**
- `gen_coExpressionNet`: constructs a gene co-expression network based on Pearson's correlations across expression profiles.
- `gen_coAbundanceNet`: constructs a gene co-abundance network based on Spearman's correlations across samples.
- `gen_coDNAmethNet`: constructs a gene network based on co-variation in DNA methylation patterns across samples.
- `gen_coCNVnet`: constructs a gene network based on co-occurrence patterns in copy number variation (CNV) data across samples.
- `gen_isolatedDrugNet`: constructs a drug network where each node represents a drug and no prior connections are assume
- `omics_network_inference`: generic function to infer biological networks from omics data. Supports different input types (e.g., expression, CNV) and inference methods (e.g., correlation, Fisher's exact test).
  
**3. Embedding**
- `gen_sim_mat_M`: applies Random Walk with Restart (RWR) to a multiplex network composed only of omics layers, integrating the information into a gene-by-gene similarity matrix. Each column of the resulting matrix represents the association scores between a given gene (the seed node) and all other genes in the multiplex network. The matrix is column-wise normalized.      
- `gen_sim_mat_MH`: applies Random Walk with Restart (RWR) to a heterogeneous multiplex network that includes both omics and drug layers, returning a (gene + drug)-by-(gene + drug) similarity matrix. Each column represents the association scores between a specific node (gene or drug) and all other nodes in the network. The matrix is column-wise normalized.       

- `get_embedding`: computes a low-dimensional embedding of the RWR similarity matrix by appling the [MultiVERSE algorithm](https://github.com/Lpiol/MultiVERSE) described by [Léo Pio-Lopez, et al.](https://arxiv.org/abs/2008.10085). 
- `get_parallel_umap_embedding`: applies UMAP in parallel for dimensionality reduction on (embedded) similarity matrix. 
- `get_pca_embedding`: applies PCA for dimensionality reduction on (embedded) similarity matrix.
- `get_tsne_embedding`: applies t-SNE in parallel for dimensionality reduction on (embedded) similarity matrix.

**4. Plotting**
- `plot_2D_matrix`: plots a 2D representation of a matrix (e.g., similarity or embedding matrix) as a scatterplot where point color and shape can be customized by providing an annotation data frame.
- `plot_traceback`: generates a bar plot to visualize the contribution of each layer to the proximity score between pairs of nodes as measured by RWR.
- `gen_traceback_net`: generates an interactive network visualization showing the connections between selected node pairs across different layers. Edge widths are proportional to proximity scores, and edges are color-coded by layer.

**Other**
- `fisher_test_post_hoc.R`: performs Fisher's exact test on a contingency table representing the status of a gene pair, and applies post hoc analysis to assess whether the observed imbalance may have biological significance. This function is internally used in `gen_coDNAmethNet.R` and `gen_coCNVnet.R`.
- `get_filtered_corMat_by_adj_pval.R`: filters a correlation matrix based on associated adjusted p-values, retaining only statistically significant correlations. This function is internally used in `gen_coExpressionNet.R` and `gen_coAbundanceNet.R`.
- `create_multiplex`: combines multiple omics-based networks into a multiplex structure, where each layer represents a different omics view of the same set of entities (e.g., genes). 
- `prune_multiplex_network`: prunes the multiplex network by removing edges whose weights fall below a specified threshold. For instance, in a correlation-based network, setting the threshold to 0.5 will remove all edges with an absolute correlation value lower than 0.5.
- `create_layer_transition_matrix`: builds a layer transition matrix for the multiplex network, defining transition probabilities between layers for random walk-based embedding. These transition probabilities are computed using the Jaccard index between layers — the higher the number of shared edges between two layers, the higher the transition probability.
- `traceback_link`: extracts layer-specific proximity values for selected node pairs from an extended RWR matrix. Optionally includes reverse direction (target → source) associations.
- `MiDNEshiny`: launches the Shiny app for summarizing the MiDNE workflow and exploring results interactively, including clustering and enrichment analyses.
- `utils`: contains utility functions used throughout the MiDNE package, including checks, data formatting, and helper functions.

&nbsp;  

### 📂 `inst/extdata/` — Data Directory
**/data**
- `/pharmacological/FDAdrugs.RDS`: a list of drugs approved by the FDA, used in the drug network analysis.
- `/pharmacological/ALLdrugs.RDS`: a comprehensive list of drugs included in the DrugBank reference set.
- `/biological/BRCA_expr_HiSeq.RDS`: expression data from TCGA-BRCA samples measured with the Illumina HiSeq technology.
- `/biological/BRCA_Methylation_Meth450.RDS`: methylation data from TCGA-BRCA samples measured with the Illumina 450K platform.
- `/biological/BRCA_proteome_CDAP.RDS`: proteomics data for TCGA-BRCA samples obtained from the CDAP pipeline.
- `/biological/BRCA_SCNA.RDS`: Somatic copy number alteration (SCNA) data for TCGA-BRCA samples.

**/annotation**
- `TCGA_BRCA_01_28_2016_ClinicalFirehose.tsi`: clinical metadata for TCGA-BRCA samples from the Broad Firehose pipeline.
- `all_genes_drugs_annotation.RDS`: an annotation table that includes both genes and drugs. For genes, the table includes information such as associated biological processes, protein complexes, and transcription factor status. For drugs, it includes annotations such as FDA approval status.

**/networks**
- `/biological/BRCA_filt_codel_network.RDS`: filtered co-deletion network of genes based on SCNA data.
- `/biological/BRCA_filt_coamp_network.RDS`: filtered co-amplification network of genes based on SCNA data.
- `/biological/BRCA_filt_cometh_network.RDS`: filtered co-methylation network based on DNA methylation profiles.
- `/biological/BRCA_filt_prot_net.RDS`: filtered co-abundance network based on proteomics data.
- `/biological/BRCA_filt_coexpr_network.RDS`: filtered gene co-expression network based on transcriptomic data.
- `/pharmacological/FDAdrugs_net.csv`: isolated drug network for FDA-approved drugs.
- `/bipartite/FDA_active_DRUGBANK_bnet.RDS`: bipartite network connecting FDA-approved drugs to their known gene targets.

**/similarity_matrices**
- `emb_FDA_active_drug_5omics_uRWRMHmat.RDS`: embedded similarity matrix computed usign Random Walk with Restart (RWR) on the BRCA heterogeneous network, followed by dimensionality reduction via the `get_embedding` function.
- `umap_emb_FDA_active_drug_5omics_uRWRMHmat.RDS`:  UMAP-based low-dimensional embedding of the above embedded similarity matrix for visualization.

---

## Contacts
If you have any questions or comments, please feel free to email Aurora Brandi (aurora.brandi@unina.it).


