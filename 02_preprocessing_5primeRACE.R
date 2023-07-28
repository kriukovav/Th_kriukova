
# packages ----------------------------------------------------------------

if (!require("here", quietly = TRUE))
  install.packages("here")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

library(here)
library(tidyverse) #v1.3.2

# ngs_directory <- here("data", "ngs") #locate the folder with ngs data
# 
# mixcr_directory <- "/home/vkriukova/software/mixcr/mixcr" #v4.1.0
# 
slurm_submit_mixcr_5RACE.sh_path <- here("slurm_submit_mixcr_5RACE.sh")

# data input --------------------------------------------------------------

filenames <- list.files(ngs_directory, pattern = "\\_CD4\\_.*TRB|\\_CD8\\_.*TRB", full.names = T) #full paths to fastq files
filenames_short <- str_remove(list.files(ngs_directory, pattern = "\\_CD4\\_.*TRB|\\_CD8\\_.*TRB", full.names = F), pattern = "_S.*") #sample names

data.frame(filenames, "sample_id" = filenames_short) %>%
  mutate(R = case_when(str_detect(.$filenames, pattern = "_R1_") ~ "R1",
                       str_detect(.$filenames, pattern = "_R2_") ~ "R2")) %>%
  pivot_wider(names_from = "R", values_from = "filenames") %>%
  relocate(sample_id, .after = R2) %>%
  write_tsv(here("filenames.txt"), col_names = F) #generate the filenames.txt file containing full paths to the ngs files and sample_ids


# mixcr analysis ----------------------------------------------------------

paste(slurm_submit_mixcr_5RACE.sh_path, mixcr_directory) -> slurm_submit_mixcr_5RACE.sh #create a variable - shell command, which runs the slurm_submit_mixcr_5RACE.sh script. This script further executes run_mixcr_5RACE.sh script in slurm batches
system(slurm_submit_mixcr_5RACE.sh) #run slurm_submit_mixcr_5RACE.sh script


# mixcr alignment report --------------------------------------------------

filenames <- list.files(here("outs", "mixcr", "stdout_files"), pattern = "\\_CD4\\_.*TRB|\\_CD8\\_.*TRB", full.names = T) #full paths to mixcr stdout files
filenames_short <- str_remove(list.files(here("outs", "mixcr", "stdout_files"), pattern = "\\_CD4\\_.*TRB|\\_CD8\\_.*TRB", full.names = F), pattern = ".txt") #sample names

t <- map(filenames, function(x) read_lines(x))
names(t) <- filenames_short

all_data <- imap(t, function(x, y) data.frame("sample_id" = y, "text" = x)) %>%
  bind_rows

all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "Total sequencing reads:")) %>%
  mutate(text = str_remove_all(text, "Total sequencing reads: "),
         text = as.double(text)) %>%
  rename(total_sequencing_reads = text) ->
  total_sequencing_reads

all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "Successfully aligned reads:")) %>%
  mutate(text = str_remove_all(text, "Successfully aligned reads: | \\(.*"),
         text = as.double(text)) %>%
  rename(total_aligned_reads = text) ->
  total_aligned_reads


all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "TRA chains:")) %>%
  slice_head() %>%
  mutate(text = str_remove_all(text, "TRA chains: ")) %>%
  separate(col = text, into = c("tra_reads_count", "tra_reads_fraction"), " ") %>%
  mutate(tra_reads_fraction = str_remove_all(tra_reads_fraction, "[\\(\\)%]"),
         tra_reads_fraction = as.double(tra_reads_fraction)) ->
  tra_reads

all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "TRB chains:")) %>%
  slice_tail() %>%
  mutate(text = str_remove_all(text, "TRB chains: ")) %>%
  separate(col = text, into = c("trb_reads_count", "trb_reads_fraction"), " ") %>%
  mutate(trb_reads_fraction = str_remove_all(trb_reads_fraction, "[\\(\\)%]"),
         trb_reads_fraction = as.double(trb_reads_fraction)) ->
  trb_reads

all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "  Number of groups accepted: ")) %>%
  mutate(text = str_remove_all(text, "  Number of groups accepted: | \\(.*"),
         text = as.double(text)) %>%
  rename(accepted_UMI = text) ->
  accepted_UMI

all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "    Effective threshold: ")) %>%
  mutate(text = str_remove_all(text, "    Effective threshold: "),
         text = as.double(text)) %>%
  rename(UMI_coverage_threshold_auto = text) ->
  UMI_coverage_threshold_auto

all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "    minRecordsPerConsensus: ")) %>%
  slice_tail(n = 1) %>%
  mutate(text = str_remove_all(text, "    minRecordsPerConsensus: "),
         text = as.double(text)) %>%
  rename(UMI_coverage_threshold_used = text) ->
  UMI_coverage_threshold_used

all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "Final clonotype count: ")) %>%
  slice_tail(n = 1) %>%
  mutate(text = str_remove_all(text, "Final clonotype count: "),
         text = as.double(text)) %>%
  rename(final_clonotype_count = text) ->
  final_clonotype_count

all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "Average number of reads per clonotype: ")) %>%
  slice_tail(n = 1) %>%
  mutate(text = str_remove_all(text, "Average number of reads per clonotype: "),
         text = as.double(text)) %>%
  rename(average_number_reads_per_clonotype = text) ->
  average_number_reads_per_clonotype

alignment_report <- total_sequencing_reads %>%
  left_join(total_aligned_reads) %>%
  left_join(tra_reads) %>%
  left_join(trb_reads) %>%
  left_join(accepted_UMI) %>%
  left_join(UMI_coverage_threshold_auto) %>%
  left_join(UMI_coverage_threshold_used) %>%
  left_join(final_clonotype_count) %>%
  left_join(average_number_reads_per_clonotype)

write_tsv(alignment_report, here("outs", "mixcr", "alignment_report_5RACE.txt"))


# collect final clonesets -------------------------------------------------

filenames <- list.files(here("outs", "mixcr"), pattern = "D11\\_CD4.*TRB\\.clones\\_TRB\\.tsv|D11\\_CD8.*TRB\\.clones\\_TRB\\.tsv", full.names = T)
filenames_short <- str_remove(list.files(here("outs", "mixcr"), pattern = "D11\\_CD4.*TRB\\.clones\\_TRB\\.tsv|D11\\_CD8.*TRB\\.clones\\_TRB\\.tsv", full.names = F), pattern = "\\.clones.*")

clonesets <- map(filenames, function(x) read_tsv(x))
names(clonesets) <- filenames_short

clonesets %>%
  map(function(x){
    select(x, uniqueUMICount, nSeqCDR3, aaSeqCDR3, allVHitsWithScore, allDHitsWithScore, allJHitsWithScore) %>%
      mutate(allVHitsWithScore = str_remove(.$allVHitsWithScore, pattern = "\\*.*")) %>%
      mutate(allJHitsWithScore = str_remove(.$allJHitsWithScore, pattern = "\\*.*")) %>%
      mutate(allDHitsWithScore = str_remove(.$allDHitsWithScore, pattern = "\\*.*")) -> temp_table
    names(temp_table) <- c("count", "cdr3nt", "cdr3aa", "v", "d", "j")
    group_by(temp_table, cdr3nt, cdr3aa, v, d, j) %>%
      summarise(count = sum(count), .groups = "drop") %>%
      filter(str_detect(.$cdr3aa, pattern = "^C.*F$")) %>% #remove CDR3 clonotypes not beginning with conserved C or ending with conserved F
      filter(!str_detect(.$cdr3aa, pattern = "\\_|\\*")) %>%  #remove CDR3 clonotypes with frameshifts or stopcodons
      mutate(freq = count/sum(count)) %>%
      relocate(count, freq, .before = "cdr3nt") %>%
      arrange(desc(count))
  }) -> clonesets_clear

dir.create(path = here("outs", "clonesets"))
iwalk(clonesets_clear, ~ write_tsv(.x, here("outs", "clonesets", paste0(.y, ".txt"))))

clonesets_clear %>%
  bind_rows(.id = "sample_id") %>%
  group_by(sample_id) %>%
  summarise(umi_sum = sum(count),
            freq_sum = sum(freq)) # check total filtered UMI count before downsample procedure




# de novo timepoint zero collection ---------------------------------------

filenames <- list.files(here("outs", "clonesets"), pattern = "D11\\_CD4.*TRB\\.txt|D11\\_CD8.*TRB\\.txt", full.names=TRUE)

filenames_short <- str_remove(list.files(here("outs", "clonesets"), pattern = "D11\\_CD4.*TRB\\.txt|D11\\_CD8.*TRB\\.txt", full.names=FALSE), pattern = ".txt")

clonesets <- map(filenames, read_tsv)
names(clonesets) <- filenames_short


# My task is to produce two pseudo-bulk PBMC files with 33k UMIs.
# The problem is I do not remember what are the replicates in this study (technical or biological?)
# Then I just remove everything that contains "repl2"

clonesets <- clonesets[!str_detect(names(clonesets), "2repl")]
clonesets_CD4_set1 <- clonesets[str_detect(names(clonesets), "CD4.*._1_|CD4.*._4_|CD4.*._5_")]
clonesets_CD4_set2 <- clonesets[str_detect(names(clonesets), "CD4.*._2_|CD4.*._3_|CD4.*._6_|CD4.*._7_")]
clonesets_CD8_set1 <- clonesets[str_detect(names(clonesets), "CD8.*._1_|CD8.*._3_|CD8.*._6_")]
clonesets_CD8_set2 <- clonesets[str_detect(names(clonesets), "CD8.*._2_|CD8.*._4_|CD8.*._5_|CD8.*._7_")]


sets <- list(clonesets_CD4_set1, clonesets_CD4_set2, clonesets_CD8_set1, clonesets_CD8_set2)
names(sets) <- c("clonesets_CD4_set1", "clonesets_CD4_set2", "clonesets_CD8_set1", "clonesets_CD8_set2")

pooled <- sets %>%
  map(bind_rows) %>%
  map(function(x) {
    group_by(x, cdr3nt, cdr3aa, v, d, j) %>%
      summarise(count = sum(count), .groups = "drop") %>%
      mutate(freq = count / sum(count)) %>%
      relocate(count, freq, .before = "cdr3nt") %>%
      arrange(desc(count))
  })

# some stats across pooled samples
pooled %>%
  map(function(x) data.frame(UMIs = sum(x$count), total_freq = sum(x$freq), n_clonotypes = length(x$count)))

DownSample <- function(data, UMI) { #data should contain a "count" column with the UMI count of each clonotype, as well as "cdr3aa", "cdr3nt", "v", "j" columns
  set.seed(42)
  data %>%
    slice(rep(1:n(), times = data$count)) %>%
    sample_n(size = UMI, replase = FALSE) %>%
    group_by(cdr3nt, cdr3aa, v, d, j) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(freq = count/sum(count)) %>%
    relocate(count, freq, .before = "cdr3nt") %>%
    arrange(desc(count))
}



ds_CD4_set1 <- DownSample(pooled$clonesets_CD4_set1, UMI = 45000*2/3)
ds_CD4_set2 <- DownSample(pooled$clonesets_CD4_set2, UMI = 45000*2/3)
ds_CD8_set1 <- DownSample(pooled$clonesets_CD8_set1, UMI = 45000/3)
ds_CD8_set2 <- DownSample(pooled$clonesets_CD8_set2, UMI = 45000/3)

#check overlaps
sum(paste(ds_CD4_set1$cdr3nt, ds_CD4_set1$cdr3aa, ds_CD4_set1$v, ds_CD4_set1$j) %in% paste(ds_CD8_set1$cdr3nt, ds_CD8_set1$cdr3aa, ds_CD8_set1$v, ds_CD8_set1$j))
sum(paste(ds_CD4_set2$cdr3nt, ds_CD4_set2$cdr3aa, ds_CD4_set2$v, ds_CD4_set2$j) %in% paste(ds_CD8_set2$cdr3nt, ds_CD8_set2$cdr3aa, ds_CD8_set2$v, ds_CD8_set2$j))


D11_PBMC_PBMCtp00_1_TRB <- bind_rows(ds_CD4_set1, ds_CD8_set1) %>%
  group_by(cdr3nt, cdr3aa, v, d, j) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(freq = count/sum(count)) %>%
  arrange(desc(count)) %>%
  relocate(count, freq, .before = cdr3nt)

D11_PBMC_PBMCtp00_2_TRB <- bind_rows(ds_CD4_set2, ds_CD8_set2) %>%
  group_by(cdr3nt, cdr3aa, v, d, j) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(freq = count/sum(count)) %>%
  arrange(desc(count)) %>%
  relocate(count, freq, .before = cdr3nt)

#save the de novo created ds45k tp00 PBMC
write_tsv(D11_PBMC_PBMCtp00_1_TRB, here("outs", "ds45k", "D11_PBMC_PBMCtp00_1_TRB.txt"))
write_tsv(D11_PBMC_PBMCtp00_2_TRB, here("outs", "ds45k", "D11_PBMC_PBMCtp00_2_TRB.txt"))
