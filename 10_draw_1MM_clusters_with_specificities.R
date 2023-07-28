library(stringdist)
library(here)
library(patchwork)
library(igraph)
library(ggnetwork)
library(tidyverse)

set.seed(11)
# data input (dataset ds45k) ----------------------------------------------

filenames <- list.files(here("outs", "ds45k"), full.names = TRUE)
filenames_short <- str_remove(list.files(here("outs", "ds45k"), full.names = F), pattern = ".txt")

clonesets <- map(filenames, read_tsv) 
names(clonesets) <- filenames_short

clonesets2 <- bind_rows(clonesets, .id = "sample_id")
clonesets2 %>% group_by(sample_id) %>% summarize(cdna = sum(count), total_freq = sum(freq), clonotypes = n()) #Check that everything is correct
rm(filenames, filenames_short, clonesets) #remove unused variables from the memory


# data input - collect all the sets of responding clones -------------------

patient_hla <- read_lines(here("data", "patient_hla.txt"))
pattern_hla <-  patient_hla %>% str_replace("\\-", "\\\\-") %>% str_replace("\\*", "\\\\*") %>% str_replace("\\:", "\\\\:") %>% paste(collapse = "|")


read_tsv(here("outs", "identified_clones", "filtered_clones", "pogorelyy2022.cdr3aa-v-j-1mm.txt")) %>% 
  filter(str_detect(.$meta.best_HLA_assoc, pattern = pattern_hla)) %>%
  filter((PBMCtp01_1 + PBMCtp01_2)/2 >= 3 |
           (PBMCtp03_1 + PBMCtp03_2)/2 >= 3 |
           (PBMCtp08_1 + PBMCtp08_2)/2 >= 3) %>% #filter only clones having at least 3 UMIs in mean_count in vax / covid1 /covid2 timepoints
  select(cdr3aa, v, j) %>% distinct -> pogorelyy2022

read_tsv(here("outs", "identified_clones", "filtered_clones", "TcellAssay.cdr3aa-cdr3nt-v-j.txt")) %>% 
  select(cdr3aa, v, j) %>% distinct -> TcellAssay

read_tsv(here("outs", "identified_clones", "filtered_clones", "vdjdb.cdr3aa-v-j-1mm.txt")) %>%
  filter(meta.species == "HomoSapiens") %>%
  filter(meta.antigen.species == "SARS-CoV-2") %>% 
  filter(str_detect(.$`meta.mhc.a`, pattern = pattern_hla) | str_detect(.$`meta.mhc.b`, pattern = pattern_hla)) %>% #filter only patient's mhc
  filter((PBMCtp01_1 + PBMCtp01_2)/2 >= 3 |
           (PBMCtp03_1 + PBMCtp03_2)/2 >= 3 |
           (PBMCtp08_1 + PBMCtp08_2)/2 >= 3) %>% #filter only clones having at least 3 UMIs in mean_count in vax / covid1 /covid2 timepoints
  filter(cdr3aa != "CASSETSGSTDTQYF") %>% # filter out irrelevant clone (selected based on very high frequency even before covid at tp00)
  select(cdr3aa, v, j) %>% distinct -> vdjdb

read_tsv(here("outs", "identified_clones", "filtered_clones", "edgeRhits.txt")) %>% 
  split(f = .$meta.comparison) %>%
  map(function(x) select(x, cdr3aa, v, j) %>% distinct) -> list_of_clones1


list_of_clones1 <- list_of_clones1[!str_detect(names(list_of_clones1), pattern = "control")]

list_of_clones2 <- list("pogorelyy2022" = pogorelyy2022, "TcellAssay" = TcellAssay, "vdjdb" = vdjdb)

list_of_clones <- c(list_of_clones1, list_of_clones2)

all_clones <- list_of_clones %>%
  bind_rows() %>%
  distinct


# find neighbors ----------------------------------------------------------

# we create a function which takes two sets of TCRs as an input (the query set - which might contain some specific interesting TCRs and the full set - which might contain the total TCR repertoire of the donor.)
# this function finds 0-1 amino acid mismatch neighbors of the query TCR set members in the full TCR set. The search of neighbors is restricted to the CDR3 clones having exactly the same V-J genes.
find_neighbors <- function(full_repertoire, query_repertoire) {
  list_vj_full_repertoire <- full_repertoire %>%
    select(cdr3aa, v, j) %>%
    distinct() %>%
    split(f = paste(.$v, .$j, sep = "_")) # split full repertoire into the list where each element is a repertoire of a certain v-j pair
  
  list_vj_query_repertoire <- query_repertoire %>%
    select(cdr3aa, v, j) %>%
    distinct() %>%
    split(f = paste(.$v, .$j, sep = "_")) # split the query repertoire into the list where each element is a repertoire of a certain v-j pair
  
  filtered_list_vj_full_repertoire <- list_vj_full_repertoire[names(list_vj_full_repertoire) %in% names(list_vj_query_repertoire)] # filter full repertoire in such a way that it only contains vj pairs from the query (for speed up)
  
  list_vj_query_repertoire <- list_vj_query_repertoire[names(list_vj_query_repertoire) %in% names(filtered_list_vj_full_repertoire)] # filter query repertoire in such a way that it only contains vj pairs from the full repertoire (most probably this step is unneccessary)
  
  distances <- map2(filtered_list_vj_full_repertoire, list_vj_query_repertoire, function(x, y) stringdistmatrix(y$cdr3aa, x$cdr3aa, method = "hamming", useNames = T))
  
  result <- map(distances, function(x){
    x %>%
      as.data.frame() %>%
      rownames_to_column(var = "query") %>%
      pivot_longer(names_to = "full", values_to = "distance", cols = starts_with("C")) %>%
      filter(distance <= 1)
  })
  
  result_table <- bind_rows(result, .id = "vj") %>%
    separate(vj, into = c("v", "j"), sep = "_") %>%
    select(full, v, j, query, distance)
  
  colnames(result_table) <- c("cdr3aa", "v", "j", "meta.query", "meta.distance")
  
  result_table
}


neighbors <- find_neighbors(full_repertoire = clonesets2, query_repertoire = all_clones) # find neighbors of the responding clones in the full TCR repertoire of the donor
write_tsv(neighbors, here("outs", "identified_clones", "filtered_clones","neighbors.txt"))

list_vj <- neighbors %>%
  select(cdr3aa, v, j) %>%
  distinct() %>%
  split(f = paste(.$v, .$j, sep = "_")) 

create_graph_dfs_from_neighbors_df <- function(neighbors_df, vj_prefix) {
  neighbors_df %>%
    select(cdr3aa) %>%
    pull() %>%
    stringdistmatrix(a = ., b = ., method = "hamming", useNames = T) %>%
    as.data.frame() %>%
    rownames_to_column(var = "var1") %>%
    pivot_longer(names_to = "var2", values_to = "distance", cols = starts_with("C")) %>%
    filter(distance <= 1) %>%
    select(var1, var2) %>%
    mutate(var1 = paste(var1, vj_prefix, sep = "_"),
           var2 = paste(var2, vj_prefix, sep = "_"))
}


graph_dfs_vj <- imap(list_vj, function(x, y) create_graph_dfs_from_neighbors_df(neighbors_df = x, vj_prefix = y))

graphs_vj <- graph_from_data_frame(bind_rows(graph_dfs_vj), directed = F)


ggnetwork_vj <- graphs_vj %>% ggnetwork(layout = igraph::with_graphopt()) %>% distinct()

ggnetwork_combined <- ggnetwork_vj %>%
  mutate(name2 = name) %>%
  separate(name2, into = c("cdr3aa", "v", "j"), sep = "_")


# add info about types of the clones --------------------------------------

antigen_names <- c("surface glycoprotein", "Spike", "S", "membrane glycoprotein",  "Matrix", "M", "nucleocapsid phosphoprotein", "N", "Nucleocapsid", "MN",
                   "ORF1ab", "ORF3", "ORF6", "NSP3", "ORF7a", "RNP", "ORF14", "ORF9b","Envelope", "ORF8", "ORF10")
meta.antigen <- c("S", "S", "S", "M", "M", "M", "N", "N", "N", "MN", rep("other", 11))

meta.colors <- c(rep("#0191FF", 3), rep("#FF009E", 3), rep("#FF9201", 3), "#50D6A1", rep("grey30", 11))

df_antigen <- data.frame(meta.antigen, antigen_names, meta.colors)



read_tsv("outs/identified_clones/filtered_clones/pogorelyy2022.cdr3aa-v-j-1mm.txt") %>% 
  filter(str_detect(.$meta.best_HLA_assoc, pattern = pattern_hla)) %>% 
  select(cdr3aa, v, j, meta.MIRA_pool_source) %>% 
  left_join(df_antigen, by = c("meta.MIRA_pool_source" = "antigen_names")) %>%
  select(cdr3aa, v, j, meta.antigen, meta.colors) %>%
  distinct -> pogorelyy2022

read_tsv("outs/identified_clones/filtered_clones/TcellAssay.cdr3aa-cdr3nt-v-j.txt") %>% 
  select(cdr3aa, v, j, meta.antigen) %>% 
  left_join(df_antigen, by = "meta.antigen") %>%
  select(cdr3aa, v, j, meta.antigen, meta.colors) %>%
  distinct -> TcellAssay

read_tsv("outs/identified_clones/filtered_clones/vdjdb.cdr3aa-v-j-1mm.txt") %>% 
  filter(meta.species == "HomoSapiens") %>%
  filter(meta.antigen.species == "SARS-CoV-2") %>% 
  filter(str_detect(.$`meta.mhc.a`, pattern = pattern_hla) | str_detect(.$`meta.mhc.b`, pattern = pattern_hla)) %>% #filter only patient's mhc
  select(cdr3aa, v, j, meta.antigen.gene) %>% 
  left_join(df_antigen, by = c("meta.antigen.gene" = "antigen_names")) %>% 
  select(cdr3aa, v, j, meta.antigen, meta.colors) %>%
  distinct -> vdjdb

neighbors %>%
  #left_join(all_clones, by = c("meta.query" = "cdr3aa", "v", "j")) %>%
  left_join(pogorelyy2022, by = c("cdr3aa", "v", "j")) %>%
  left_join(TcellAssay, by = c("cdr3aa", "v", "j"), suffix = c(".pogorelyy2022", ".TcellAssay")) %>%
  left_join(vdjdb, by = c("cdr3aa", "v", "j")) %>%
  rename("meta.antigen.vdjdb" = "meta.antigen",
         "meta.colors.vdjdb" = "meta.colors") %>%
  select(-meta.query, -meta.distance) %>%
  distinct() -> clone_meta
  
ggnetwork_combined %>%
  left_join(clone_meta, by = c("cdr3aa", "v", "j")) -> plot_data

plot_data$meta.colors.pogorelyy2022[is.na(plot_data$meta.colors.pogorelyy2022)] <- "grey90"
plot_data$meta.colors.vdjdb[is.na(plot_data$meta.colors.vdjdb)] <- "grey90"
plot_data$meta.colors.TcellAssay[is.na(plot_data$meta.colors.TcellAssay)] <- "grey90"

plot_data %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges()+
  geom_nodes(aes(color = meta.colors.vdjdb), shape = 21, stroke = 1, alpha = 0.7,  size = 2) +
  geom_nodes(aes(color = meta.colors.pogorelyy2022), shape = 21, stroke = 1, alpha = 0.7,  size = 1.3) +
  geom_nodes(aes(color = meta.colors.TcellAssay), alpha = 0.7,  size = 0.9) +
  scale_color_identity(labels = unique(df_antigen$meta.antigen), breaks = unique(df_antigen$meta.colors), guide = "legend") +
  theme_void() +
  coord_fixed() -> plots


ggsave(plot = plots, filename = paste("all_clustersSpecificities", "pdf", sep = "."), path = here("outs", "identified_clones", "figures"), height = 4.5, width = 4.5)
write_rds(plots, paste(here("outs", "identified_clones", "figures"), "all_clustersSpecificities.rds", sep = "/"))

