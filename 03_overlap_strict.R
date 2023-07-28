
# packages ----------------------------------------------------------------

if (!require("patchwork", quietly = TRUE))
  install.packages("patchwork")

library(here)
library(patchwork)
library(tidyverse)
library(ggrepel)

dir.create(here("outs", "identified_clones", "statistics"), recursive = T)
dir.create(here("outs", "identified_clones", "figures"), recursive = T)
dir.create(here("outs", "identified_clones", "filtered_clones"), recursive = T)

set.seed(42)

# data input (dataset ds45k) ----------------------------------------------

filenames <- list.files(here("outs", "ds45k"), full.names = TRUE)
filenames_short <- str_remove(list.files(here("outs", "ds45k"), full.names = F), pattern = ".txt")

clonesets <- map(filenames, read_tsv) 
names(clonesets) <- filenames_short

clonesets2 <- bind_rows(clonesets, .id = "sample_id")
clonesets2 %>% group_by(sample_id) %>% summarize(cdna = sum(count), total_freq = sum(freq), clonotypes = n()) #Check that everything is correct
rm(filenames, filenames_short, clonesets) #remove unused variables from the memory

# data input (table of responding clones) ---------------------------------------

# path_to_file <- here("data", "TcellAssay_fdr_fc_filtered/clones_TcellAssay.txt")
# 
# table_of_responding_clones <- read_tsv(path_to_file) #should contain 'cdr3aa', 'cdr3nt', 'v', 'j' columns, and any number of additional columns starting with 'meta.'
# 
# assay_name <- "TcellAssay" #will be used in the output naming
# 
# intersect_type <- "cdr3aa-cdr3nt-v-j" # what is the rule to intersect table_of_responding_clones with clonesets2? Only affects naming of the output files
# 
# facetting_vars <- c("meta.subset", "meta.antigen") # used only for plotting and the generation of statistics file. Please eneter any number of variables from table_of_responding_clones table to allow plot facetting by these variables

# convert ds45k dataset into wide format ----------------------------------

clonesets2_wide <- clonesets2 %>%
  mutate(timepoint_replicate = str_extract(.$sample_id, pattern = "PBMCtp[0-9][0-9]\\_[0-9]")) %>% # use regex to grep info about timepoint and replicate from the sample name
  pivot_wider(names_from = timepoint_replicate, values_from = count, id_cols = c("cdr3aa", "cdr3nt", "v", "j"), values_fill = 0)
  
clonesets2_wide_filtered <- table_of_responding_clones %>%
  inner_join(clonesets2_wide, by = c("cdr3aa", "cdr3nt", "v", "j"))

write_tsv(clonesets2_wide_filtered, here("outs", "identified_clones", "filtered_clones", paste(assay_name, intersect_type, "txt", sep = ".")))

clonesets2_wide_filtered %>%
  group_by(across({{facetting_vars}})) %>%
  summarize(count = n(), .groups = "drop") %>%
  write_tsv(here("outs", "identified_clones", "statistics", paste(assay_name, intersect_type, "txt", sep = "."))) # counts of identified clones 
  

# prepare data for plotting -----------------------------------------------

clonesets2_long_filtered <- clonesets2_wide_filtered %>%
  pivot_longer(names_to = "timepoint_replicate", values_to = "count", cols = starts_with("PBMC")) %>%
  separate(timepoint_replicate, into = c("timepoint", "replicate"), sep = "_") %>%
  group_by(cdr3aa, cdr3nt, v, j, timepoint) %>%
  mutate(replicate_mean = mean(count)) %>%
  mutate(plotting_line_group = paste(cdr3aa, cdr3nt, v, j, sep = "_")) # we prepared data for plotting
  

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
  group_by(cdr3aa, cdr3nt, v, j) %>%
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

data_repel <- inner_join(clonesets2_long_filtered, top_five_clones) %>%
  filter(timepoint == "PBMCtp08") %>%
  select(timepoint, line_color, replicate_mean, cdr3aa, cdr3nt, v, j, {{facetting_vars}}) %>%
  distinct()

clonesets2_long_filtered %>%
  ggplot() +
  geom_vline(data = plot_labels, aes(xintercept = timepoint), color = "grey70", lty = "dashed") +
  geom_line(aes(x = timepoint, y = replicate_mean, group = plotting_line_group), alpha = 0.2, size = 0.4) +
  geom_line(data = inner_join(clonesets2_long_filtered, top_five_clones), aes(x = timepoint, y = replicate_mean, group = plotting_line_group, color = line_color), size = 0.4) +
  geom_label_repel(data = data_repel, aes(x = timepoint, y = replicate_mean, label = cdr3aa, color = line_color),
                   size = 2, fill = "white", alpha = 0.9,
                   direction = "y", xlim = c(11, Inf),
                   min.segment.length = Inf) +
  scale_x_discrete(expand = c(c(0, 0), c(0,4))) +
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





