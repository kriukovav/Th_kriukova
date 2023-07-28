# Th_kriukova

This is an analysis pipeline used to   
a) extract TCR repertoires from raw sequncing reads   
b) analyze TCR reperoires and draw all TCR-repertoires related figures from Kriukova et al.    
c) draw scRNA-Seq related figures from Kriukova et al.  

For easier navigation through the script please refer to main.R script. 
You may run the pipeline starting from the raw reads (mixcr and slurm software; raw reads from NCBI SRA bioproject: PRJNA995237 are required), or starting from the processed TCR data (please check comments in main.R for details).

## What data files are here?
For an easier usage, I added to the data folder:  
a) vdjdb-2022-03-30 folder, containing [vdjdb](https://vdjdb.cdr3.net) database   
b) TcellAssay_fdr_fc_filtered folder, containing clonotypes expanded in the T cell fucntional assay (data analysis was performed by Ksenia Lupyr)   
c) gating_V2 folder, containing theresults of in silico flow-cytometry like gating of scRNA-Seq cells    
d) patient_hla.txt - list of HLA alleles, D11   
e) pogorelyy2022.txt - supplementary file from [Pogorelyy et al.](10.1016/j.xcrm.2022.100697)   
f) sample_dates_V2.txt - list of samples from D11.

I added to the outs folder:  
a) ds45k folder - contains processed downsampled to the equal number of UMIs (45,000) bulk TCR beta repertoires from PBMC (required to run the pipeline starting from the processed data)   

## Before you start
a) Please install all the needed R packages

```(r)
install.packages("tidyverse", "Seurat", "here", "pals", "patchwork", "ggrepel", "stringdist", "scatterpie", "ggVennDiagram", "igraph", "ggnetwork, lubridate")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
```
b) Please download our reference peripheral blood CD4 T lymphocytes scRNA-Seq dataset (link). This is an rds file containing an integrated Seurat object (full_reference_return_model.rds). Put it into outs/scRNAseq folder.    

c) Please download the processed scRNA-Seq datasets from effector-memory CD4 T cells in D11 (links). These are rds files containing Seurat objects (full_reference_return_model.rds). Put it into outs/scRNAseq folder.    

d) Always work from the project folder (Th_kriukova)




### In case you prefer to run from raw reads:
Install [mixcr](https://mixcr.com/mixcr/getting-started/installation/) (v4.1.0)
Check if slurm is supported by your administrators. Maybe you will need to add corrections to my shell scripts.  
I recommend to create an ngs folder as a subdirectory of data folder, and to put the files from SRA to data/ngs. 




