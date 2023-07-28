
# packages ----------------------------------------------------------------
library(here)
library(tidyverse)


# data input --------------------------------------------------------------

read_tsv(here("data", "pogorelyy2022.txt")) -> pogorelyy2022


# reformat data -----------------------------------------------------------

pogorelyy2022$best_HLA_assoc -> hla
hla <- unique(hla)

as.data.frame(hla) %>%
  filter(str_detect(hla, "HLA-D.AB\\*..\\:.._..\\:..")) %>%
  tidyr::extract(hla, into = c("HLA", "A","B"), "HLA\\-D([A-Z])AB\\*([0-9][0-9]\\:[0-9][0-9])\\_([0-9][0-9]\\:[0-9][0-9])", remove = F) %>%
  mutate(hla_clean = paste0("HLA-D", .$HLA, "A1", "*", .$A, "_", "HLA-D", .$HLA, "B1", "*", .$B)) %>%
  select(hla, hla_clean) -> part1 # I expand the HLA data for HLA-DQ/DP alpha and beta chains

as.data.frame(hla) %>%
  filter(str_detect(hla, "HLA-DRDQ\\*..\\:.._..\\:.._..\\:..")) %>%
  tidyr::extract(hla,into = c("HLA", "HLA2", "B","A", "B2"), "HLA\\-D([A-Z])D([A-Z])\\*([0-9][0-9]\\:[0-9][0-9])\\_([0-9][0-9]\\:[0-9][0-9])\\_([0-9][0-9]\\:[0-9][0-9])", remove = F) %>%
  mutate(hla_clean = paste0("HLA-D", .$HLA, "B1", "*", .$B, "_", "HLA-D", .$HLA2, "A1", "*", .$A, "_", "HLA-D", .$HLA2, "B1", "*", .$B2)) %>%
  select(hla, hla_clean) -> part2 # I expand the HLA data for HLA-DRB1/DQ haplotypes

HLA_labels_clean <- bind_rows(part1, part2) %>%
  full_join(as.data.frame(hla)) %>% # I add those HLA labels, that are already written in the expanded form
  mutate(hla_clean = case_when(is.na(hla_clean) ~ hla,
                               TRUE ~ hla_clean))

pogorelyy2022 <- pogorelyy2022 %>%
  left_join(HLA_labels_clean, by = c("best_HLA_assoc" = "hla")) %>%
  mutate(best_HLA_assoc = hla_clean) %>%
  select(-hla_clean)


pogorelyy2022 %>% 
  mutate(vb = str_remove(.$vb, "\\*.*")) %>%
  mutate(jb = str_remove(.$jb, "\\*.*")) %>%
  rename(v = vb,
         j = jb,
         cdr3aa = cdr3b) %>%
  relocate(cdr3aa, v, j) -> pogorelyy2022_reformated

names(pogorelyy2022_reformated)[4:length(names(pogorelyy2022_reformated))] <- paste0("meta.", names(pogorelyy2022_reformated)[4:length(names(pogorelyy2022_reformated))])

# summarize cluster info  -------------------------------------------------

pogorelyy2022_reformated %>% 
  distinct(meta.tcrdist120_cluster_id, meta.va) %>% 
  group_by(meta.tcrdist120_cluster_id) %>%
  summarize(unique_va_count = n()) -> unique_va_count # create a df with number of unique Va segments in the clusters

pogorelyy2022_reformated %>%
  distinct(meta.tcrdist120_cluster_id, v) %>% 
  group_by(meta.tcrdist120_cluster_id) %>%
  summarize(unique_vb_count = n()) -> unique_vb_count # create a df with number of unique Vb segments in the clusters

pogorelyy2022_reformated %>%
  group_by(meta.tcrdist120_cluster_id) %>%
  summarize(cluster_size = n()) -> cluster_size # create a df with cluster sizes

cluster_size %>%
  left_join(unique_va_count, by = "meta.tcrdist120_cluster_id") %>%
  left_join(unique_vb_count, by = "meta.tcrdist120_cluster_id") %>%
  mutate(ab_cluster_type = case_when((unique_va_count - unique_vb_count) >= 4 ~ "b",
                                     (unique_va_count - unique_vb_count) <= -4 ~ "a",
                                     TRUE ~ "ab")) %>% 
  filter(ab_cluster_type == "a") %>%
  select(meta.tcrdist120_cluster_id) %>% pull -> a_cluster # cluster, which specificity is determined by alpha chain - I will exclude this cluster from analysis, as we only study TCRbeta repertoires
  

pogorelyy2022_reformated %>%
  distinct(meta.tcrdist120_cluster_id, meta.best_HLA_assoc) %>%
  group_by(meta.tcrdist120_cluster_id) %>%
  summarize(hla_hits_number = n()) %>%
  filter(hla_hits_number > 1) %>%
  select(meta.tcrdist120_cluster_id) %>% pull -> clusters_with_nonunique_hla # I found cluster IDs which are composed of TCRs with different predicted HLA restriction

pogorelyy2022_reformated %>%
  filter(meta.tcrdist120_cluster_id %in% clusters_with_nonunique_hla) %>% 
  select(1:10, meta.best_HLA_assoc) %>% 
  split(f = .$meta.tcrdist120_cluster_id) %>%
  map(function(x) x %>% 
        group_by(meta.best_HLA_assoc) %>% 
        summarise(n = n()) %>% # count how many cluster-members are restricted to the given HLA
        arrange(desc(n)) %>% # arrange data so that the most represented in cluster HLAs are on top of the table for the given cluster
        filter(n >= mean(n)) %>% # I only take overepresented in the cluster HLAs. I am going to assign the restriction to these HLA to the whole cluster
        select(meta.best_HLA_assoc) %>% pull %>% paste(collapse = "_") %>%
        data.frame("meta.clusterHLA" = .)) %>%
  bind_rows(.id = "meta.tcrdist120_cluster_id") %>%
  mutate(meta.tcrdist120_cluster_id = as.numeric(meta.tcrdist120_cluster_id)) -> expected_hla_for_clusters

pogorelyy2022_reformated %>%
  left_join(expected_hla_for_clusters) %>% # add HLA info for clusters with nonunique initial best HLA assoc
  mutate(meta.clusterHLA = case_when(is.na(meta.clusterHLA) ~ meta.best_HLA_assoc,
                                          TRUE ~ meta.clusterHLA)) %>%
  mutate(meta.best_HLA_assoc = meta.clusterHLA) %>%
  filter(! meta.tcrdist120_cluster_id %in% a_cluster) %>% #filter out alpha clusters
  select(- meta.clusterHLA) -> pogorelyy2022_reformated_refined
  
  



# save data ---------------------------------------------------------------

write_tsv(pogorelyy2022_reformated, here("data", "pogorelyy2022_reformated.txt"))
write_tsv(pogorelyy2022_reformated_refined, here("data", "pogorelyy2022_reformated_refined.txt"))

