# This is a pipeline script which calls the analysis scripts
library(here)
library(pals)


#1 (Start from here to analyze raw data)
ngs_directory <- here::here("data", "ngs") #locate the folder with ngs data
mixcr_directory <- "/home/vkriukova/software/mixcr/mixcr" #v4.1.0 #locate your mixcr command

source("01_preprocessing.R")
rm(list = ls(all.names = TRUE)) #will clear all objects 

#2
ngs_directory <- here::here("data", "ngs") #locate the folder with ngs data
mixcr_directory <- "/home/vkriukova/software/mixcr/mixcr" #v4.1.0 #locate your mixcr command

source("02_preprocessing_5primeRACE.R")
rm(list = ls(all.names = TRUE)) #will clear all objects 

#3 (Start from here to analyze processed TCR repertoires)
path_to_file <- here::here("data", "TcellAssay_fdr_fc_filtered", "clones_TcellAssay.txt")
table_of_responding_clones <- readr::read_tsv(path_to_file) #should contain 'cdr3aa', 'cdr3nt', 'v', 'j' columns, and any number of additional columns starting with 'meta.'
assay_name <- "TcellAssay" #will be used in the output naming
intersect_type <- "cdr3aa-cdr3nt-v-j" # what is the rule to intersect table_of_responding_clones with clonesets2? Only affects naming of the output files 
facetting_vars <- c("meta.subset", "meta.antigen")

source("03_overlap_strict.R")
rm(list = ls(all.names = TRUE)) #will clear all objects

#4 
source("04_PBMC_tracking_edgeR.R")
rm(list = ls(all.names = TRUE)) #will clear all objects

#5
source("05_reformat_pogorelyy2022.R")
rm(list = ls(all.names = TRUE)) #will clear all objects

#6
source("06_reformat_vdjdb.R")
rm(list = ls(all.names = TRUE)) #will clear all objects

#7a
path_to_file <- here::here("data", "vdjdb-2022-03-30", "vdjdb_full_reformated.txt")
table_of_responding_clones <- readr::read_tsv(path_to_file) #should contain 'cdr3aa', 'v', 'j' columns, and any number of additional columns starting with 'meta.'
assay_name <- "vdjdb" #will be used in the output naming
intersect_type <- "cdr3aa-v-j-1mm" # what is the rule to intersect table_of_responding_clones with clonesets2? Only affects naming of the output files
facetting_vars <- c("meta.antigen.species", "meta.mhc.class") # used only for plotting and the generation of statistics file. Please eneter any number of variables from table_of_responding_clones table to allow plot facetting by these variables

source("07_overlap_1MMcdr3aaVJ.R")
rm(list = ls(all.names = TRUE)) #will clear all objects

#7b
path_to_file <- here::here("data", "pogorelyy2022_reformated_refined.txt")
table_of_responding_clones <- readr::read_tsv(path_to_file) #should contain 'cdr3aa', 'v', 'j' columns, and any number of additional columns starting with 'meta.'
assay_name <- "pogorelyy2022" #will be used in the output naming
intersect_type <- "cdr3aa-v-j-1mm" # what is the rule to intersect table_of_responding_clones with clonesets2? Only affects naming of the output files
facetting_vars <- c("meta.MIRA_pool_source", "meta.best_HLA_assoc") # used only for plotting and the generation of statistics file. Please eneter any number of variables from table_of_responding_clones table to allow plot facetting by these variables

source("07_overlap_1MMcdr3aaVJ.R")
rm(list = ls(all.names = TRUE)) #will clear all objects

#8
source("08_figures_pogorelyy_vdjdb.R")
rm(list = ls(all.names = TRUE)) #will clear all objects

#9
source("09_draw_1MM_clusters.R")
rm(list = ls(all.names = TRUE)) #will clear all objects

#10
source("10_draw_1MM_clusters_with_specificities.R")
rm(list = ls(all.names = TRUE)) #will clear all objects

#11
source("11_overlap_clones.R")
rm(list = ls(all.names = TRUE)) #will clear all objects

#12
source("12_overlap_clusters.R")
rm(list = ls(all.names = TRUE)) #will clear all objects

#13
source("13_plot_timeline.R")
rm(list = ls(all.names = TRUE)) #will clear all objects


#14
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")

color_code <- data.frame("number" = as.character(0:15), "colors" = alphabet2()[1:16]) %>%
  mutate(colors = case_when(number == "8" ~ "red",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive_activated", "EffMem_Th22", "Naive",
                                             "CentMem", "Tfh", "Naive_IFN_response", "SelfStopping",
                                             "EffMem_Th2a", "EffMem_IFN_response", "Temra_cytotoxic_Th1", "Cycling"))



source("14_plot_scRNAseq_heatmap.R")
rm(list = ls(all.names = TRUE)) #will clear all objects

#15
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
folder_with_gating_results <- here("data", "gating_V2")
color_code <- data.frame("number" = as.character(0:15), "colors" = alphabet2()[1:16]) %>%
  mutate(colors = case_when(number == "8" ~ "red",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive_activated", "EffMem_Th22", "Naive",
                                             "CentMem", "Tfh", "Naive_IFN_response", "SelfStopping",
                                             "EffMem_Th2a", "EffMem_IFN_response", "Temra_cytotoxic_Th1", "Cycling"))


source("15_dot_plot_comparison_sort_cite.R")
rm(list = ls(all.names = TRUE)) #will clear all objects


#16
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
folder_with_gating_results <- here("data", "gating_V2")

source("16_cardiff_coloring.R")
rm(list = ls(all.names = TRUE)) #will clear all objects


#17
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
color_code <- data.frame("number" = as.character(0:15), "colors" = alphabet2()[1:16]) %>%
  mutate(colors = case_when(number == "8" ~ "red",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive_activated", "EffMem_Th22", "Naive",
                                             "CentMem", "Tfh", "Naive_IFN_response", "SelfStopping",
                                             "EffMem_Th2a", "EffMem_IFN_response", "Temra_cytotoxic_Th1", "Cycling"))


source("17_plot_clonality.R")
rm(list = ls(all.names = TRUE)) #will clear all objects



#18
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
color_code <- data.frame("number" = as.character(0:15), "colors" = alphabet2()[1:16]) %>%
  mutate(colors = case_when(number == "8" ~ "red",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive_activated", "EffMem_Th22", "Naive",
                                             "CentMem", "Tfh", "Naive_IFN_response", "SelfStopping",
                                             "EffMem_Th2a", "EffMem_IFN_response", "Temra_cytotoxic_Th1", "Cycling"))


source("18_annotation.R") 
rm(list = ls(all.names = TRUE)) #will clear all objects



#19
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
path_to_query <- c(here("outs", "scRNAseq", "D11_COVID1.rds"),
                   here("outs", "scRNAseq", "D11_COVID2_rep1.rds"),
                   here("outs", "scRNAseq", "D11_COVID2_rep2.rds"))

color_code <- data.frame("number" = as.character(0:15), "colors" = alphabet2()[1:16]) %>%
  mutate(colors = case_when(number == "8" ~ "red",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive_activated", "EffMem_Th22", "Naive",
                                             "CentMem", "Tfh", "Naive_IFN_response", "SelfStopping",
                                             "EffMem_Th2a", "EffMem_IFN_response", "Temra_cytotoxic_Th1", "Cycling"))


source("19_label_transfer.R") 
rm(list = ls(all.names = TRUE)) #will clear all objects


#20
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
path_to_query <- here("outs", "scRNAseq", "query_scRNA2.rds")

color_code <- data.frame("number" = as.character(0:15), "colors" = alphabet2()[1:16]) %>%
  mutate(colors = case_when(number == "8" ~ "red",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive_activated", "EffMem_Th22", "Naive",
                                             "CentMem", "Tfh", "Naive_IFN_response", "SelfStopping",
                                             "EffMem_Th2a", "EffMem_IFN_response", "Temra_cytotoxic_Th1", "Cycling"))

source("20_scatterpie.R") 
rm(list = ls(all.names = TRUE)) #will clear all objects

#21
source("21_collect_table_of_specific_clones.R") 
rm(list = ls(all.names = TRUE)) #will clear all objects

#22
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")

color_code <- data.frame("number" = as.character(0:15), "colors" = alphabet2()[1:16]) %>%
  mutate(colors = case_when(number == "8" ~ "red",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive_activated", "EffMem_Th22", "Naive",
                                             "CentMem", "Tfh", "Naive_IFN_response", "SelfStopping",
                                             "EffMem_Th2a", "EffMem_IFN_response", "Temra_cytotoxic_Th1", "Cycling"))


