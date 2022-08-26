**High level description**

Model for correcting technical variation in Single-Cell gene expression data
Iterative normalization and clustering method for single-cell gene expression data.

**What is single cell seq?**

The emerging technology of single-cell RNA-seq gives access to gene expression measurements for thousands of cells, allowing discovery and characterization of cell types. Single cell sequqncing presents exciting opportunities to study heterogeneity of expression and characterize unknown cell types. This contrasts traditional bulk gene expression data where the gene expression is measured by an average readout across a bulk of cells.

**Problem with algorithm:**

The data is confounded by technical variation emanating from experimental errors and cell type-specific biases.
Current approaches perform a global normalization prior to analyzing biological signals, which does not resolve missing data or variation dependent on latent cell types.

**Problem with data:**

* One primary reason that makes single-cell RNA-seq analysis challenging is dropouts, where the data only captures a small fraction of the transcriptome of each cell. Almost all computational algorithms developed for single-cell RNA-seq adopted gene selection, dimension reduction or imputation to address the dropouts.
* Biases in cell sampling
* Significant differences in total number of mRNA molecules, as well as variation in library size, defined as sum of amplified mRNA molecules per cell.
* Signals become more stable when individual signals are summarized (such as in a bulk experiment); thus, the increase in resolution due to sc-seq also means a reduction of the stability of the supporting signals.

**How do we solve it:**

Our model is formulated as a hierarchical Bayesian mixture model with cell-specific scalings that aid the iterative normalization and clustering of cells, teasing apart technical variation from biological signals. We demonstrate that this approach is superior to global normalization followed by clustering.

**How do we measure the impact:**

We show identifiability and weak convergence guarantees of our method and present a scalable Gibbs inference algorithm.
This method improves cluster inference in both synthetic and real single-cell data compared with previous methods, and allows easy interpretation and recovery of the underlying structure and cell types.

**Data Collection:**

The sample data is collected in form of (cells * genes), where cells are rows and genes are represented as columns. The values are the gene counts for each cell and gene combination.
The gene data was from the cell type in mouse cortex and hippocampus(http://linnarssonlab.org/cortex/).


**Model:**

We used Dirichlet Process Mixture Model for (https://towardsdatascience.com/tl-dr-dirichlet-process-gaussian-mixture-models-made-easy-12b4d492e5f9) for clustering of genes.

**Why:**

When we have limited prior belief over the cluster number or mixing probabilities, we can turn non-parametric and consider infinitely many of them. Naturally many of these clusters will be redundant, and have mixing probabilities so close to 0 that we can just ignore them. Incredibly, this framework lets the data determine the most likely number of clusters.

* Data-smoothing methods define a “similarity” between cells (e.g., cells that are neighbors in a graph or occupy a small region in a latent space) and adjust expression values for each cell based on expression values in similar cells. These methods usually adjust all expression values, including technical zeros, biological zeros, and observed non-zero values.

* Local normalization techniques used for normlaization based on clusters of cell data. Normalization was required for significant diffrences in library size.

* After Data normalization, we used PCA : PRINICIAPL COMPONENT ANALYSIS for dimensionality reduction of data.


**Output:**

We evaluated BISCUIT’s performance on real world data using mouse cortex cells from Zeisel et al. (2015) that include ground truth labels for 7 known cell types. For computational speed we chose d = 558 genes with largest standard deviation across n = 3005 cells. Figure S7 shows the confusion matrix for inferred classes and Figure 8 shows the mode of inferred classes across 500 Gibbs sweeps post burn-in of 1500 sweeps compared to their actual cell type
labels. Cells are visualized using t-SNE dimensionality reduction (Van der Maaten & Hinton, 2008), as this was
shown to be an effective visualization that captures the cell type structure in single-cell data (Amir et al., 2013).

=======
This repository contains source code for the paper: [Dirichlet Process Mixture Model for Correcting Technical Variation in Single-Cell Gene Expression Data by Sandhya Prabhakaran*, Elham Azizi*, Ambrose J Carr, and Dana Pe’er](http://jmlr.org/proceedings/papers/v48/prabhakaran16.pdf)  in ICML 2016

Software provided in R.

**Dependencies**

1. R - install R from https://cran.r-project.org/. A version higher than 2.12.0 is recommended.
    sudo yum install R
2. Once installed, open a terminal and at the command prompt, type R. 
3. At the R prompt: Install the following R packages by issuing command:

>install.packages(c("MCMCpack","mvtnorm","ellipse","coda","Matrix","Rtsne","gtools","foreach","doParallel","doSNOW","snow","lattice","MASS","bayesm","robustbase","chron","mnormt","schoolmath","devtools","RColorBrewer"))


> **Note:**
Until here, this is a one-time activity. 


**How to run code**


1. Clone/Download code repository.
2. At the R prompt: issue command
>rm(list=ls());
>graphics.off()
3. Issue setwd() to point to the path where the code repository resides. 
eg. if the code is downloaded at “/User/Downloads/BISCUIT/", then type
>working_path <- "/Users/Downloads/BISCUIT/"; 
>setwd(working_path);

4. In start_file.R:

	a. input_file_name: is the name of your input data available as counts. File must be of the form cells x genes. You can also download data used in the ICML paper from [here](https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt).  (Data source: http://linnarssonlab.org/cortex/). 
	>input_file_name <- ‘expression_mRNA_17-Aug-2014.txt’

	b. input_data_tab_delimited: Set to TRUE if the input data is tab-delimited
	
	c. is_format_genes_cells: Set to TRUE if the input data has rows as genes and columns as cells.
	
	d. choose_cells: choose the number of cells or comment out to use all the cells in the input dataset.  
	
   	e. choose_genes: choose the number of genes or comment out to use all the genes in the input dataset.  This will select genes based on the ordered Fiedler vector. For eg., choose_genes <- 20, for the top 20 genes ordered by the magnitude of the Fiedler vector.

	f. gene_batch: set it such that 20 <= gene_batch <= 150. (This will create gene batches across which the Infinite Mixture model will run in parallel)
	
	g. num_iter: Maximum number of MCMC iterations for convergence. Set this such that 15 <= num_iter <= ~50. 

   	h. num_cores: Set this to a value lesser than the total number of cores in your device. For eg, in R, type detectCores(). If this returns a value greater than 1 then set num_cores <- detectCores() - 1, else set num_cores <- 1.
    
    i. z_true_labels_avl: set to TRUE if true labels of cells are available, else set to FALSE.

    j. num_cells_batch: required to split the data in feasible sets for parallel processing the confusion matrix. Set this to 1000 if input number of cells is in the 1000s, else set it to 100. 
	
	k. alpha: DPMM dispersion parameter. A higher value spins more clusters whereas a lower value spins lesser clusters. For the Zeisel et al. data, alpha = 0.005 or less. 

	l. output_folder_rename: give a prefix to rename your existing /output/ folder, if any.

5. At the R prompt: Issue command
> source("start_file.R")

**Output**
______

1. Output folder gets created which has the run log.txt and log_CM.txt files. For resolving any discrepancies while parallel processing of your data, check debug.txt and debug_CM.txt.
2. output/plots/ folder gives the various plots
3. In the R terminal, the latent variables of interest can be obtained by issuing:

>- z_inferred_final               ##  *(class assignment of the cells)*
>- alpha_inferred_final        ##   *(inferred cell-specific alphas)* 
>- beta_inferred_final          ##  *(inferred cell-specific betas)*
>- mu_final                         ##  *(inferred mus. dim(mu_final) <- numgenes x K)*
>- Sigma_final                    ##  *(inferred Sigmas)*
>- total_clusters                ## *(total inferred clusters)*

4. The following are also saved as .txt files for further analysis

>- ~/output/plots/extras/pre_imputed_tSNE_coord.txt  (tSNE coordinates of data pre-imputation)
>- ~/output/plots/extras/post_imputed_tSNE_coord.txt (tSNE coordinates of data post-imputation)
>- ~/output/plots/Inferred_Sigmas/Sigma_final_k.txt (Inferred Sigma matrices per cluster k)
>- ~/output/plots/Inferred_means/mu_final (Inferred mean matrix; rows are numgenes, columns are number of total_clusters)
>- ~/output/plots/Inferred_Sigmas/Genes_selected.csv (Genes selected based on the global co-expression/gene disparity/Fiedler vector)
>- ~/output/plots/Inferred_means/Genes_selected.csv (Genes selected based on the global co-expression/gene disparity/Fiedler vector)
>- ~/output/plots/Inferred_alphas_betas/Final_alphas.csv (Inferred alpha values)
>- ~/output/plots/Inferred_alphas_betas/Final_betas.csv (Inferred beta values)
>- ~/output/plots/extras/Imputed_Y_countspace.txt (Imputed data in count space)
>- ~/output/plots/extras/Imputed_Y_logspace.txt (Imputed data in log space)
>- ~/output/plots/extras/Input_parms_used.txt 
>- ~/output/plots/extras/Output_values.txt (total number of clusters, cluster proportions)
>- ~/output/plots/extras/cluster_posterior_probabilities.csv (posterior probabilities of every cell to every inferred cluster)


**Opening an issue/contacting the developers**
---------------------
1. Under the repository, click Issues.
2. Click New issue.
3. Type a title and description for your issue, new feature you wish to see added etc.
4. When you are finished, click Submit new issue.

