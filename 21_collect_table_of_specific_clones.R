library(tidyverse)

# data input - collect all the sets of responding clones -------------------

patient_hla <- read_lines(here("data", "patient_hla.txt"))
pattern_hla <-  patient_hla %>% str_replace("\\-", "\\\\-") %>% str_replace("\\*", "\\\\*") %>% str_replace("\\:", "\\\\:") %>% paste(collapse = "|")


read_tsv(here("outs", "identified_clones", "filtered_clones", "pogorelyy2022.cdr3aa-v-j-1mm.txt")) %>% 
  filter(str_detect(.$meta.best_HLA_assoc, pattern = pattern_hla)) %>%
  filter((PBMCtp01_1 + PBMCtp01_2)/2 >= 3 |
           (PBMCtp03_1 + PBMCtp03_2)/2 >= 3 |
           (PBMCtp08_1 + PBMCtp08_2)/2 >= 3) %>% #filter only clones having at least 3 UMIs in mean_count in vax / covid1 /covid2 timepoints
  select(cdr3aa, v, j, meta.MIRA_pool_source) %>% distinct() -> pogorelyy2022

read_tsv(here("outs", "identified_clones", "filtered_clones", "TcellAssay.cdr3aa-cdr3nt-v-j.txt")) %>%
  select(cdr3aa, v, j, cdr3nt, meta.antigen, meta.subset) %>% distinct -> TcellAssay

read_tsv(here("outs", "identified_clones", "filtered_clones", "vdjdb.cdr3aa-v-j-1mm.txt")) %>%
  filter(meta.species == "HomoSapiens") %>%
  filter(meta.antigen.species == "SARS-CoV-2") %>% 
  filter(str_detect(.$`meta.mhc.a`, pattern = pattern_hla) | str_detect(.$`meta.mhc.b`, pattern = pattern_hla)) %>% #filter only patient's mhc
  filter((PBMCtp01_1 + PBMCtp01_2)/2 >= 3 |
           (PBMCtp03_1 + PBMCtp03_2)/2 >= 3 |
           (PBMCtp08_1 + PBMCtp08_2)/2 >= 3) %>% #filter only clones having at least 3 UMIs in mean_count in vax / covid1 /covid2 timepoints
  filter(cdr3aa != "CASSETSGSTDTQYF") %>% # filter out irrelevant clone (selected based on very high frequency even before covid at tp00)
  select(cdr3aa, v, j, meta.antigen.gene) %>% distinct -> vdjdb

read_tsv(here("outs", "identified_clones", "filtered_clones", "edgeRhits.txt")) %>%
  filter(!str_detect(.$meta.comparison, "control")) %>%
  select(cdr3aa, v, j, cdr3nt, meta.comparison) %>%
  distinct() -> timetracking


TcellAssay <- TcellAssay %>%
  rename("TcellAssay.cdr3nt" = "cdr3nt", "TcellAssay.meta.antigen" = "meta.antigen", "TcellAssay.meta.subset" = "meta.subset") %>%
  group_by(cdr3aa, v, j) %>%
  summarize(TcellAssay.cdr3nt = paste(TcellAssay.cdr3nt),
            TcellAssay.meta.antigen = paste(TcellAssay.meta.antigen),
            TcellAssay.meta.subset = paste(TcellAssay.meta.subset)) 

pogorelyy2022 <- pogorelyy2022 %>%
  rename("pogorelyy2022.meta.antigen" = "meta.MIRA_pool_source") %>%
  mutate(pogorelyy2022.meta.antigen = case_when(pogorelyy2022.meta.antigen == "membrane glycoprotein" ~ "M",
                                                pogorelyy2022.meta.antigen == "surface glycoprotein" ~ "S",
                                                TRUE ~ pogorelyy2022.meta.antigen))


vdjdb <- vdjdb %>%
  rename("vdjdb.meta.antigen" = "meta.antigen.gene") %>%
  mutate(vdjdb.meta.antigen = case_when(vdjdb.meta.antigen == "Spike" ~ "S",
                                        vdjdb.meta.antigen == "Nucleocapsid" ~ "N",
                                        vdjdb.meta.antigen == "Matrix" ~ "M"))

timetracking <- timetracking %>%
  rename("timetracking.meta.comparison" = "meta.comparison",
         "timetracking.cdr3nt" = "cdr3nt")

list_of_SarsCoV2_specific_clones <- TcellAssay %>%
  full_join(timetracking, by = c("cdr3aa", "v", "j")) %>%
  full_join(vdjdb, by = c("cdr3aa", "v", "j")) %>%
  full_join(pogorelyy2022, by = c("cdr3aa", "v", "j"))

distinct_clonotypes <- list_of_SarsCoV2_specific_clones %>%
  select(cdr3aa, v, j) %>%
  distinct() 

distinct_clonotypes$clonotype_id_number <- 1:length(distinct_clonotypes$cdr3aa)

distinct_clonotypes %>%
  full_join(list_of_SarsCoV2_specific_clones, by = c("cdr3aa", "v", "j")) %>% 
  mutate(donor_id = "D11") %>%
  relocate(clonotype_id_number, .before = cdr3aa) -> list_of_SarsCoV2_specific_clones_v2

write_tsv(list_of_SarsCoV2_specific_clones_v2, here::here("outs", "identified_clones", "filtered_clones", "supplementary_table_for_paper.tsv"))
  



