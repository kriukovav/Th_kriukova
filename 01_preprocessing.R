
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
slurm_submit_mixcr.sh_path <- here("slurm_submit_mixcr.sh")

# data input --------------------------------------------------------------

filenames <- list.files(ngs_directory, pattern = "PBMC", full.names = T) #full paths to fastq files
filenames_short <- str_remove(list.files(ngs_directory, pattern = "PBMC", full.names = F), pattern = "_S.*") #sample names

data.frame(filenames, "sample_id" = filenames_short) %>%
  mutate(R = case_when(str_detect(.$filenames, pattern = "_R1_") ~ "R1",
                       str_detect(.$filenames, pattern = "_R2_") ~ "R2")) %>%
  pivot_wider(names_from = "R", values_from = "filenames") %>%
  relocate(sample_id, .after = R2) %>%
  write_tsv(here("filenames.txt"), col_names = F) #generate the filenames.txt file containing full paths to the ngs files and sample_ids

# mixcr analysis ----------------------------------------------------------

paste(slurm_submit_mixcr.sh_path, mixcr_directory) -> slurm_submit_mixcr.sh #create a variable - shell command, which runs the slurm_submit_mixcr.sh script. This script further executes run_mixcr.sh script in slurm batches
system(slurm_submit_mixcr.sh) #run slurm_submit_mixcr.sh script


# mixcr alignment report --------------------------------------------------

filenames <- list.files(here("outs", "mixcr", "stdout_files"), pattern = "PBMC", full.names = T) #full paths to mixcr stdout files
filenames_short <- str_remove(list.files(here("outs", "mixcr", "stdout_files"), pattern = "PBMC", full.names = F), pattern = ".txt") #sample names

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

write_tsv(alignment_report, here("outs", "mixcr", "alignment_report.txt"))


# collect final clonesets -------------------------------------------------

filenames <- list.files(here("outs", "mixcr"), pattern = "TRB\\.clones\\_TRB\\.tsv", full.names = T)
filenames_short <- str_remove(list.files(here("outs", "mixcr"), pattern = "TRB\\.clones\\_TRB\\.tsv", full.names = F), pattern = "\\.clones.*")

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

ds45k <- map(clonesets_clear, function(x) DownSample(data = x, UMI = 45000))

dir.create(path = here("outs", "ds45k"))
iwalk(ds45k, ~ write_tsv(.x, here("outs", "ds45k", paste0(.y, ".txt"))))











