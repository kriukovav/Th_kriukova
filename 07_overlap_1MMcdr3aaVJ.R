
# packages ----------------------------------------------------------------

library(stringdist)
library(here)
library(patchwork)
library(tidyverse)

dir.create(here("outs", "identified_clones", "statistics"), recursive = T)
dir.create(here("outs", "identified_clones", "figures"), recursive = T)
dir.create(here("outs", "identified_clones", "filtered_clones"), recursive = T)

# data input (dataset ds45k) ----------------------------------------------

filenames <- list.files(here("outs", "ds45k"), full.names = TRUE)
filenames_short <- str_remove(list.files(here("outs", "ds45k"), full.names = F), pattern = ".txt")

clonesets <- map(filenames, read_tsv) 
names(clonesets) <- filenames_short

clonesets2 <- bind_rows(clonesets, .id = "sample_id")
clonesets2 %>% group_by(sample_id) %>% summarize(cdna = sum(count), total_freq = sum(freq), clonotypes = n()) #Check that everything is correct
rm(filenames, filenames_short, clonesets) #remove unused variables from the memory

# data input (table of responding clones) ---------------------------------------

# path_to_file <- here("data", "vdjdb-2022-03-30", "vdjdb_full_reformated.txt")
# 
# table_of_responding_clones <- read_tsv(path_to_file) #should contain 'cdr3aa', 'v', 'j' columns, and any number of additional columns starting with 'meta.'
# 
# assay_name <- "vdjdb" #will be used in the output naming
# 
# intersect_type <- "cdr3aa-v-j-1mm" # what is the rule to intersect table_of_responding_clones with clonesets2? Only affects naming of the output files 
# 
# facetting_vars <- c("meta.antigen.species", "meta.mhc.class") # used only for plotting and the generation of statistics file. Please eneter any number of variables from table_of_responding_clones table to allow plot facetting by these variables

# convert ds45k dataset into wide format ----------------------------------

clonesets2_wide <- clonesets2 %>%
  mutate(timepoint_replicate = str_extract(.$sample_id, pattern = "PBMCtp[0-9][0-9]\\_[0-9]")) %>% # use regex to grep info about timepoint and replicate from the sample name
  pivot_wider(names_from = timepoint_replicate, values_from = count, id_cols = c("cdr3aa", "v", "j"), values_fill = 0, values_fn = sum) # I sum counts of convergent TCR-CDR3 clonotypes

table_of_responding_clones <- table_of_responding_clones %>%
  select(cdr3aa, v, j, starts_with("meta.")) %>%
  distinct() # remove multiple entries of convergent clonotypes in the table_of_responding_clones before inner_join


# intersect by cdr3aa+v+j +1mm --------------------------------------------

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


neighbors <- find_neighbors(full_repertoire = clonesets2_wide, query_repertoire = table_of_responding_clones) # find neighbors of the responding clones in the full TCR repertoire of the donor

clonesets2_wide_filtered <- neighbors %>%
  inner_join(clonesets2_wide, by = c("cdr3aa", "v", "j")) %>% # filter only those TCR clonotypes, which are present within neighbors set
  left_join(table_of_responding_clones, by = c("meta.query" = "cdr3aa", "v", "j"))  %>% # add some missing metadata related to the query responding TCRs - so that neighbors automatically inherit the same metadata
  relocate(starts_with("meta."), .after = "j") # relocate columns

write_tsv(clonesets2_wide_filtered, here("outs", "identified_clones", "filtered_clones", paste(assay_name, intersect_type, "txt", sep = ".")))

clonesets2_wide_filtered %>%
  group_by(across({{facetting_vars}})) %>%
  summarize(count = n(), .groups = "drop") %>%
  write_tsv(here("outs", "identified_clones", "statistics", paste(assay_name, intersect_type, "txt", sep = "."))) # counts of identified clones 


# prepare data for plotting -----------------------------------------------

clonesets2_long_filtered <- clonesets2_wide_filtered %>%
  pivot_longer(names_to = "timepoint_replicate", values_to = "count", cols = starts_with("PBMC")) %>%
  separate(timepoint_replicate, into = c("timepoint", "replicate"), sep = "_") %>%
  group_by(cdr3aa, v, j, timepoint) %>%
  mutate(replicate_mean = mean(count)) %>%
  mutate(plotting_line_group = paste(cdr3aa, v, j, sep = "_")) # we prepared data for plotting


data.frame(timepoint = c("PBMCtp01", "PBMCtp03", "PBMCtp08"),
           condition = c("VAX", "COVID", "COVID"),
           color = c("darkorchid3", "#009988", "#009988")) -> plot_labels # data frame with additional info for plotting

read_tsv(here("outs", "identified_clones", "statistics", paste(assay_name, intersect_type, "txt", sep = "."))) -> plot_statistics


ymax = max(clonesets2_long_filtered$count) + max(clonesets2_long_filtered$count) * 0.5
ytext =  ymax * 0.7
ynumbers = ymax * 0.7

top_five_clones <- clonesets2_long_filtered %>%
  select(-count, -replicate) %>%
  distinct() %>%
  group_by(cdr3aa, v, j) %>%
  slice_max(n = 1, order_by = replicate_mean, with_ties = F) %>%
  ungroup() %>%
  group_by(across({{facetting_vars}})) %>%
  slice_max(n = 5, order_by = replicate_mean, with_ties = F) %>%
  ungroup()

line_palette <- c("#0077bb", "#33bbee", "#009988", "#EE7733", "#cc3311")

#if number of lines to be colored on the plots is not multiple of the number of colors in color palette. Adjust the length of color palette vecor by adding N more colors
if ((length(top_five_clones$plotting_line_group)%%length(line_palette)) != 0) top_five_clones$line_color <- c(rep(line_palette, length(top_five_clones$plotting_line_group)%/%length(line_palette)), line_palette[1 : (length(top_five_clones$plotting_line_group)%%length(line_palette))])

#if number of lines to be colored on the plots is multiple of the number of colors in color palette. Just rep the color palette vector desired number of times
if ((length(top_five_clones$plotting_line_group)%%length(line_palette)) == 0) top_five_clones$line_color <- rep(line_palette, length(top_five_clones$plotting_line_group)%/%length(line_palette))


top_five_clones <- top_five_clones %>%
  select({{facetting_vars}}, plotting_line_group, line_color)

clonesets2_long_filtered %>%
  ggplot() +
  geom_vline(data = plot_labels, aes(xintercept = timepoint), color = "grey70", lty = "dashed") +
  geom_line(aes(x = timepoint, y = replicate_mean, group = plotting_line_group), alpha = 0.2, size = 0.4) +
  geom_line(data = inner_join(clonesets2_long_filtered, top_five_clones), aes(x = timepoint, y = replicate_mean, group = plotting_line_group, color = line_color), size = 0.4) +
  #geom_point(aes(x = timepoint, y = count), size = 0.5) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x=element_text(angle = 90),
        panel.grid.major.y = element_line(size = 0.05, linetype = 'solid', colour = "grey90"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste(assay_name, intersect_type, "within ds45k dataset")) +
  ylab("clone size (UMI count)") +
  xlab("time point") +
  facet_wrap(facets = facetting_vars) +
  geom_text(data = plot_labels, aes(x = timepoint, label = condition, color = color),
            y = ytext, size = 4, angle = 90) +
  scale_color_identity() +
  geom_text(data = plot_statistics,
            aes(label = paste("n =", count)),
            x = 5.5, y = ynumbers, size = 4) +
  ylim(0, ymax) -> p


ggsave(plot = p, filename = paste(assay_name, intersect_type, "pdf", sep = "."), path = here("outs", "identified_clones", "figures"), height = 10, width = 10)





