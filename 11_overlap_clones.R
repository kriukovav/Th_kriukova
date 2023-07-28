
# packages ----------------------------------------------------------------
library(here)
library(ggVennDiagram)
library(tidyverse)

set.seed(42)

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

vd_data_1 <- list_of_clones1[!str_detect(names(list_of_clones1), pattern = "control")] %>%
  map(function(x) {
    x %>%
      mutate(clone_id = paste(cdr3aa, v, j, sep = "_")) %>%
      select(clone_id) %>%
      pull()
  })

vd_data_2 <- list("pogorelyy2022" = pogorelyy2022, "TcellAssay" = TcellAssay, "vdjdb" = vdjdb) %>%
  map(function(x) {
    x %>%
      mutate(clone_id = paste(cdr3aa, v, j, sep = "_")) %>%
      select(clone_id) %>%
      pull()
  })

vd_data_3 <- c(vd_data_1, vd_data_2)

vd_data_1_pooled <- list("timetracking" = unique(unlist(vd_data_1, use.names = F)))

vd_data_4 <- c(vd_data_1_pooled, vd_data_2)

vd_data_4 <- vd_data_4[sort(names(vd_data_4))]

plot <- ggVennDiagram(vd_data_4, label = "count") + 
  scale_x_continuous(expand = expansion(mult = .2)) + 
  scale_y_continuous(expand = expansion(mult = .2)) + 
  scale_fill_gradient(low="white",high = "pink") +
  ggtitle("CDR3aa + v + j overlaps") +
  theme(legend.position = "none")

ggsave(plot = plot, filename = "venn_four_groups.pdf", path = here("outs", "identified_clones", "figures"), height = 7, width = 7)
write_rds(plot, paste(here("outs", "identified_clones", "figures"), "venn_four_groups.rds", sep = "/"))


