
# packages ----------------------------------------------------------------
library(here)
library(tidyverse)


# data input --------------------------------------------------------------

read_tsv(here("data", "vdjdb-2022-03-30", "vdjdb_full.txt")) -> vdjdb


# reformat data -----------------------------------------------------------

vdjdb %>%
  mutate(v.beta = str_remove(.$v.beta, "\\*.*")) %>%
  mutate(j.beta = str_remove(.$j.beta, "\\*.*")) %>%
  rename(v = v.beta,
         j = j.beta,
         cdr3aa = cdr3.beta) %>%
  relocate(cdr3aa, v, j) -> vdjdb_reformated
  
names(vdjdb_reformated)[4:length(names(vdjdb_reformated))] <- paste0("meta.", names(vdjdb_reformated)[4:length(names(vdjdb_reformated))])


# save data ---------------------------------------------------------------

write_tsv(vdjdb_reformated, here("data", "vdjdb-2022-03-30", "vdjdb_full_reformated.txt"))
  