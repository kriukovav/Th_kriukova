library(here)
library(tidyverse)
library(ggrepel)

set.seed(42)


# data input - collect all the sets of resonding clones -------------------

patient_hla <- read_lines(here("data", "patient_hla.txt"))
pattern_hla <-  patient_hla %>% str_replace("\\-", "\\\\-") %>% str_replace("\\*", "\\\\*") %>% str_replace("\\:", "\\\\:") %>% paste(collapse = "|")


read_tsv(here("outs", "identified_clones", "filtered_clones", "pogorelyy2022.cdr3aa-v-j-1mm.txt")) %>% 
  filter(str_detect(.$meta.best_HLA_assoc, pattern = pattern_hla)) %>%
  filter((PBMCtp01_1 + PBMCtp01_2)/2 >= 3 |
           (PBMCtp03_1 + PBMCtp03_2)/2 >= 3 |
           (PBMCtp08_1 + PBMCtp08_2)/2 >= 3) %>% #filter only clones having at least 3 UMIs in mean_count in vax / covid1 /covid2 timepoints
  select(cdr3aa, v, j, meta.MIRA_pool_source, starts_with("PBMC")) %>% distinct() -> pogorelyy2022

read_tsv(here("outs", "identified_clones", "filtered_clones", "vdjdb.cdr3aa-v-j-1mm.txt")) %>%
  filter(meta.species == "HomoSapiens") %>%
  filter(meta.antigen.species == "SARS-CoV-2") %>% 
  filter(str_detect(.$`meta.mhc.a`, pattern = pattern_hla) | str_detect(.$`meta.mhc.b`, pattern = pattern_hla)) %>% #filter only patient's mhc
  filter((PBMCtp01_1 + PBMCtp01_2)/2 >= 3 |
           (PBMCtp03_1 + PBMCtp03_2)/2 >= 3 |
           (PBMCtp08_1 + PBMCtp08_2)/2 >= 3) %>% #filter only clones having at least 3 UMIs in mean_count in vax / covid1 /covid2 timepoints
  filter(cdr3aa != "CASSETSGSTDTQYF") %>% # filter out irrelevant clone (selected based on very high frequency even before covid at tp00)
  select(cdr3aa, v, j, meta.antigen.gene, starts_with("PBMC")) %>% distinct -> vdjdb

# prepare data for plotting pogorelyy2022 -----------------------------------------------


data.frame(timepoint = c("PBMCtp01", "PBMCtp03", "PBMCtp08"),
           condition = c("VAX", "COVID", "COVID"),
           color = c("darkorchid3", "#009988", "#009988")) -> plot_labels # data frame with additional info for plotting


clonesets2_wide_filtered <- pogorelyy2022
facetting_vars <- colnames(clonesets2_wide_filtered)[str_detect(colnames(clonesets2_wide_filtered), pattern = "meta.")]
assay_name <- "CellRepMed.Pogorelyy.et.al" #will be used in the output naming
intersect_type <- "cdr3aa-v-j-1mm" # what is the rule to intersect table_of_responding_clones with clonesets2? Only affects naming of the output files 


clonesets2_long_filtered <- clonesets2_wide_filtered %>%
  pivot_longer(names_to = "timepoint_replicate", values_to = "count", cols = starts_with("PBMC")) %>%
  separate(timepoint_replicate, into = c("timepoint", "replicate"), sep = "_") %>%
  group_by(cdr3aa, v, j, timepoint) %>%
  mutate(replicate_mean = mean(count)) %>%
  mutate(plotting_line_group = paste(cdr3aa, v, j, sep = "_")) # we prepared data for plotting


clonesets2_wide_filtered %>%
  group_by(across({{facetting_vars}})) %>%
  summarize(count = n(), .groups = "drop") -> plot_statistics

clonesets2_long_filtered %>%
  group_by(across({{facetting_vars}})) %>%
  summarize(max_count = max(count),
            ymax = max(count) + 0.1 * max(count),
            ytext = ymax * 0.9,
            ynumbers = ymax * 0.9) -> plot_y_axis # calculate the values for ymax, and y values to plot the text on plots

merge(plot_y_axis, plot_labels) -> plot_labels # introduce individual text y values for each facet
merge(plot_y_axis, plot_statistics) -> plot_statistics # introduce individual text y values for each facet

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

data_repel <- inner_join(clonesets2_long_filtered, top_five_clones) %>%
  filter(timepoint == "PBMCtp08") %>%
  select(timepoint, line_color, replicate_mean, cdr3aa, v, j, {{facetting_vars}}) %>%
  distinct()

clonesets2_long_filtered %>%
  ggplot() +
  geom_vline(data = plot_labels, aes(xintercept = timepoint), color = "grey70", lty = "dashed") +
  geom_line(aes(x = timepoint, y = replicate_mean, group = plotting_line_group), alpha = 0.2, size = 0.4) +
  geom_line(data = inner_join(clonesets2_long_filtered, top_five_clones), aes(x = timepoint, y = replicate_mean, group = plotting_line_group, color = line_color), size = 0.4) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x=element_text(angle = 90),
        panel.grid.major.y = element_line(size = 0.05, linetype = 'solid', colour = "grey90"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  ggtitle(paste(assay_name, intersect_type, "within ds45k dataset")) +
  ylab("clone size (UMI count)") +
  xlab("time point") +
  facet_wrap(facets = facetting_vars, scales = "free_y") +
  geom_blank(data = plot_y_axis, aes(y = ymax)) + # I need this to introduce individual ylim values for each facet
  geom_text(data = plot_labels, aes(x = timepoint, label = condition, color = color, y = ytext),
            size = 5, angle = 90) +
  scale_color_identity() +
  geom_text(data = plot_statistics,
            aes(label = paste("n =", count), y = ynumbers),
            x = 6, size = 5) +
  geom_label_repel(data = data_repel, aes(x = timepoint, y = replicate_mean, label = cdr3aa, color = line_color),
                   size = 4, fill = "white", alpha = 0.9,
                   direction = "y", xlim = c(11, Inf),
                   min.segment.length = Inf) +
  scale_x_discrete(expand = c(c(0, 0), c(0,7)))

ggsave(plot = last_plot(), filename = paste("HLAmatched", assay_name, intersect_type, "pdf", sep = "."), path = here("outs", "identified_clones", "figures"), height = 4, width = 15)



# prepare data for plotting vdjdb -----------------------------------------



data.frame(timepoint = c("PBMCtp01", "PBMCtp03", "PBMCtp08"),
           condition = c("VAX", "COVID", "COVID"),
           color = c("darkorchid3", "#009988", "#009988")) -> plot_labels # data frame with additional info for plotting



clonesets2_wide_filtered <- vdjdb
facetting_vars <- colnames(clonesets2_wide_filtered)[str_detect(colnames(clonesets2_wide_filtered), pattern = "meta.")]
assay_name <- "vdjdb" #will be used in the output naming
intersect_type <- "cdr3aa-v-j-1mm" # what is the rule to intersect table_of_responding_clones with clonesets2? Only affects naming of the output files 


clonesets2_long_filtered <- clonesets2_wide_filtered %>%
  pivot_longer(names_to = "timepoint_replicate", values_to = "count", cols = starts_with("PBMC")) %>%
  separate(timepoint_replicate, into = c("timepoint", "replicate"), sep = "_") %>%
  group_by(cdr3aa, v, j, timepoint) %>%
  mutate(replicate_mean = mean(count)) %>%
  mutate(plotting_line_group = paste(cdr3aa, v, j, sep = "_")) # we prepared data for plotting

clonesets2_wide_filtered %>%
  group_by(across({{facetting_vars}})) %>%
  summarize(count = n(), .groups = "drop") -> plot_statistics

clonesets2_long_filtered %>%
  group_by(across({{facetting_vars}})) %>%
  summarize(max_count = max(count),
            ymax = max(count) + 0.1 * max(count),
            ytext = ymax * 0.9,
            ynumbers = ymax * 0.9) -> plot_y_axis # calculate the values for ymax, and y values to plot the text on plots

merge(plot_y_axis, plot_labels) -> plot_labels # introduce individual text y values for each facet
merge(plot_y_axis, plot_statistics) -> plot_statistics # introduce individual text y values for each facet

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

data_repel <- inner_join(clonesets2_long_filtered, top_five_clones) %>%
  filter(timepoint == "PBMCtp08") %>%
  select(timepoint, line_color, replicate_mean, cdr3aa, v, j, {{facetting_vars}}) %>%
  distinct()

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
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  ggtitle(paste(assay_name, intersect_type, "within ds45k dataset")) +
  ylab("clone size (UMI count)") +
  xlab("time point") +
  facet_wrap(facets = facetting_vars, scales = "free_y", nrow = 1) +
  geom_blank(data = plot_y_axis, aes(y = ymax)) + # I need this to introduce individual ylim values for each facet
  geom_text(data = plot_labels, aes(x = timepoint, label = condition, color = color, y = ytext),
            size = 5, angle = 90) +
  scale_color_identity() +
  geom_text(data = plot_statistics,
            aes(label = paste("n =", count), y = ynumbers),
            x = 6, size = 5) +
  geom_label_repel(data = data_repel, aes(x = timepoint, y = replicate_mean, label = cdr3aa, color = line_color),
                   size = 4, fill = "white", alpha = 0.9,
                   direction = "y", xlim = c(11, Inf),
                   min.segment.length = Inf) +
  scale_x_discrete(expand = c(c(0, 0), c(0,7)))

ggsave(plot = last_plot(), filename = paste("HLAmatched", assay_name, intersect_type, "pdf", sep = "."), path = here("outs", "identified_clones", "figures"), height = 4, width = 15)
