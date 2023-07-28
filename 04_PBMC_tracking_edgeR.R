
# packages ----------------------------------------------------------------


library(edgeR)
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)

set.seed(42)


# load data ds45k --------------------------------------------------------

filenames <- list.files(here("outs", "ds45k"), pattern = "PBMC", full.names = T)
filenames_short <- str_remove(list.files(here("outs", "ds45k"), pattern = "PBMC", full.names = F), pattern = ".txt")


clonesets <- map(filenames, function(x) read_tsv(x))
names(clonesets) <- filenames_short



# prepare data for edgeR --------------------------------------------------

# prepare data for edgeR - wide-format matrix with samples in columns and cDNA counts for cdr3aa-cdr3nt-v-j clonotypes in rows
clones_matrix <- clonesets %>%
  bind_rows(.id = "sample_id") %>%
  mutate(clone_id = paste(cdr3aa, cdr3nt, v, j, sep = "_")) %>%
  pivot_wider(id_cols = "clone_id", names_from = sample_id, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "clone_id") %>%
  as.matrix()

# create a list with pattern to filter columns in the clones_matrix, and prepare a list of vectors, indicating experimental and reference samples for edgeR. Both lists contain 7 elements - because we do 7 edgeR runs
VAX_response <- c("PBMCtp01", "PBMCtp02", "PBMCtp04", "PBMCtp07") %>% paste(collapse = "|")
groups_VAX_response <- c(rep("expansion", 2), rep("reference", 6))

COVID1_response <- c("PBMCtp02", "PBMCtp03", "PBMCtp04", "PBMCtp07") %>% paste(collapse = "|")
groups_COVID1_response <- c(rep("reference", 2), rep("expansion", 2), rep("reference", 4))

COVID2_response <- c("PBMCtp02", "PBMCtp04", "PBMCtp07", "PBMCtp08") %>% paste(collapse = "|")
groups_COVID2_response <- c(rep("reference", 6), rep("expansion", 2))

control1_response <- c("PBMCtp02", "PBMCtp04", "PBMCtp05", "PBMCtp07") %>% paste(collapse = "|")
groups_control1_response <- c(rep("expansion", 2), rep("reference", 6))

control2_response <- c("PBMCtp02", "PBMCtp04", "PBMCtp05", "PBMCtp07") %>% paste(collapse = "|")
groups_control2_response <- c(rep("reference", 2), rep("expansion", 2), rep("reference", 4))

control3_response <- c("PBMCtp02", "PBMCtp04", "PBMCtp05", "PBMCtp07") %>% paste(collapse = "|")
groups_control3_response <- c(rep("reference", 4), rep("expansion", 2), rep("reference", 2))

control4_response <- c("PBMCtp02", "PBMCtp04", "PBMCtp05", "PBMCtp07") %>% paste(collapse = "|")
groups_control4_response <- c(rep("reference", 6), rep("expansion", 2))

list_timepoints <- list(VAX_response, COVID1_response, COVID2_response, control1_response, control2_response, control3_response, control4_response)
names(list_timepoints) <- c("VAX_response", "COVID1_response", "COVID2_response", "control1_response", "control2_response", "control3_response", "control4_response")

list_edgeRgroups <- list(groups_VAX_response, groups_COVID1_response, groups_COVID2_response, groups_control1_response, groups_control2_response, groups_control3_response, groups_control4_response)
names(list_edgeRgroups) <- c("groups_VAX_response", "groups_COVID1_response", "groups_COVID2_response", "groups_control1_response", "groups_control2_response", "groups_control3_response", "groups_control4_response")

# function which 1) filters only necessary columns in the wide-format clonotype data, assigns them into experimental and reference groups, and runs edgeR pipeline
do_edgeR_analysis <- function(count_matrix, pattern_colnames, vector_edgeRgroups, nTags = 100, filter_FDR = 0.05, filterFC = -1){
  count_matrix_filtered <- count_matrix[, colnames(count_matrix)[str_detect(colnames(count_matrix), pattern = pattern_colnames)]]
  clones_DGE <- DGEList(count_matrix_filtered, group = vector_edgeRgroups)
  clones_DGE <- calcNormFactors(clones_DGE)
  clones_DGE <- estimateDisp(clones_DGE)
  et <- exactTest(clones_DGE)
  result <- topTags(et, n = nTags)$table %>%
    filter(FDR < filter_FDR & logFC < filterFC)
  result
}


# run edgeR and plot results ----------------------------------------------

edgeR_hits <- map2(list_timepoints, list_edgeRgroups, function(x, y) do_edgeR_analysis(count_matrix = clones_matrix, pattern_colnames = x, vector_edgeRgroups = y)) # run edgeR
edgeR_hits <- edgeR_hits[(map(edgeR_hits, function(x) length(row.names(x))) %>% unlist()) > 0] # I filter only those elements of the list, which contain at least one hit

edgeR_hits %>%
  bind_rows(.id = "meta.comparison") %>%
  group_by(meta.comparison) %>%
  summarize(count = n(), .groups = "drop") %>% 
  write_tsv(file = here("outs", "identified_clones", "statistics", paste("edgeRhits", "txt", sep = "."))) # saved edgeR hits statistics

edgeR_hits %>%
  map(function(x) {
    rename_with(x, .fn = ~paste("meta", ., sep = ".")) %>%
    rownames_to_column(var = "clone_id")
    }) %>%
  bind_rows(.id = "meta.comparison") %>%
  separate(clone_id, into = c("cdr3aa", "cdr3nt", "v", "j"), sep = "_") %>%
  write_tsv(file = here("outs", "identified_clones", "filtered_clones", paste("edgeRhits", "txt", sep = "."))) # saved edgeR hits 


data.frame(timepoint = c("PBMCtp01", "PBMCtp03", "PBMCtp08"),
           condition = c("VAX", "COVID", "COVID"),
           color = c("darkorchid3", "#009988", "#009988")) -> plot_labels # data frame with additional info for plotting

read_tsv(here("outs", "identified_clones", "statistics", paste("edgeRhits", "txt", sep = "."))) -> plot_statistics
read_tsv(here("outs", "identified_clones", "filtered_clones", paste("edgeRhits", "txt", sep = "."))) -> edgeR_hits

# prepare data in long format for plotting
clonesets2_long_filtered <- clones_matrix %>%
  as.data.frame %>%
  rownames_to_column(var = "clone_id") %>%
  inner_join(edgeR_hits %>% mutate(clone_id = paste(cdr3aa, cdr3nt, v, j, sep = "_")), by = "clone_id") %>% 
  pivot_longer(names_to = "sample_id", values_to = "count", cols = starts_with("D1")) %>%
  mutate(timepoint = str_extract(.$sample_id, pattern = "PBMCtp[0-9][0-9]")) %>%
  mutate(replicate = str_extract(.$sample_id, pattern = "\\_[0-9]\\_") %>% str_extract("[0-9]")) %>%
  group_by(clone_id, timepoint) %>%
  mutate(replicate_mean = mean(count)) %>%
  mutate(plotting_line_group = clone_id) %>%
  ungroup()# we prepared data for plotting

clonesets2_long_filtered %>%
  group_by(meta.comparison) %>%
  summarize(max_count = max(count),
            ymax = max(count) + 0.1 * max(count),
            ytext = ymax * 0.9,
            ynumbers = ymax * 0.9) -> plot_y_axis # calculate the values for ymax, and y values to plot the text on plots

merge(plot_y_axis, plot_labels) -> plot_labels # introduce individual text y values for each facet
merge(plot_y_axis, plot_statistics) -> plot_statistics # introduce individual text y values for each facet

top_five_clones <- clonesets2_long_filtered %>%
  select(-count, -replicate) %>%
  distinct() %>%
  group_by(clone_id) %>%
  slice_max(n = 1, order_by = replicate_mean, with_ties = F) %>%
  ungroup() %>%
  group_by(meta.comparison) %>%
  slice_max(n = 5, order_by = replicate_mean, with_ties = F) %>%
  ungroup()

line_palette <- c("#0077bb", "#33bbee", "#009988", "#EE7733", "#cc3311")

#if number of lines to be colored on the plots is not multiple of the number of colors in color palette. Adjust the length of color palette vecor by adding N more colors
if ((length(top_five_clones$plotting_line_group)%%length(line_palette)) != 0) top_five_clones$line_color <- c(rep(line_palette, length(top_five_clones$plotting_line_group)%/%length(line_palette)), line_palette[1 : (length(top_five_clones$plotting_line_group)%%length(line_palette))])

#if number of lines to be colored on the plots is multiple of the number of colors in color palette. Just rep the color palette vector desired number of times
if ((length(top_five_clones$plotting_line_group)%%length(line_palette)) == 0) top_five_clones$line_color <- rep(line_palette, length(top_five_clones$plotting_line_group)%/%length(line_palette))


top_five_clones <- top_five_clones %>%
  select(meta.comparison, plotting_line_group, line_color)

data_repel <- inner_join(clonesets2_long_filtered, top_five_clones) %>%
  filter(timepoint == "PBMCtp08") %>%
  select(timepoint, line_color, replicate_mean, clone_id, meta.comparison) %>%
  distinct() %>%
  separate(clone_id, sep = "_", into = c("cdr3aa"), extra = "drop")

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
  theme_minimal() +
  theme(legend.position = "none", axis.text.x=element_text(angle = 90),
        panel.grid.major.y = element_line(size = 0.05, linetype = 'solid', colour = "grey90"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("EdgeR hits") +
  ylab("clone size (UMI count)") +
  xlab("time point") +
  facet_wrap(facets = vars(meta.comparison), scales = "free_y", nrow = 2) +
  geom_blank(data = plot_y_axis, aes(y = ymax)) + # I need this to introduce individual ylim values for each facet
  geom_text(data = plot_labels, aes(x = timepoint, label = condition, color = color, y = ytext),
            size = 4, angle = 90) +
  scale_color_identity() +
  geom_text(data = plot_statistics,
            aes(label = paste("n =", count), y = ynumbers),
            x = 5.5, size = 4)



ggsave(plot = last_plot(), filename = paste("EdgeR", "pdf", sep = "."), path = here("outs", "identified_clones", "figures"), height = 7, width = 18.5)


# The commented code further creates a figure optimized for the paper (including fonts, sizes, etc). It is not necessary to run it

# clonesets2_long_filtered %>%
#   filter(!str_detect(.$meta.comparison, "control")) %>%
#   ggplot() +
#   geom_vline(data = plot_labels %>% filter(!str_detect(.$meta.comparison, "control")), aes(xintercept = timepoint), color = "grey70", lty = "dashed") +
#   geom_line(aes(x = timepoint, y = replicate_mean, group = plotting_line_group), alpha = 0.2, size = 0.4) +
#   geom_line(data = inner_join(clonesets2_long_filtered %>% filter(!str_detect(.$meta.comparison, "control")),
#                               top_five_clones), aes(x = timepoint, y = replicate_mean, group = plotting_line_group, color = line_color), size = 0.4) +
#   geom_label_repel(data = data_repel %>% filter(!str_detect(.$meta.comparison, "control")), aes(x = timepoint, y = replicate_mean, label = cdr3aa, color = line_color),
#                    size = 4, fill = "white", alpha = 0.9,
#                    direction = "y", xlim = c(11, Inf),
#                    min.segment.length = Inf) +
#   scale_x_discrete(expand = c(c(0, 0), c(0,12))) +
#   theme_minimal() +
#   theme(legend.position = "none", axis.text.x=element_text(angle = 90),
#         panel.grid.major.y = element_line(size = 0.05, linetype = 'solid', colour = "grey90"),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         strip.text = element_text(size = 14)) +
#   ggtitle("EdgeR hits") +
#   ylab("clone size (UMI count)") +
#   xlab("time point") +
#   facet_wrap(facets = vars(meta.comparison), nrow = 1) +
#   geom_blank(data = plot_y_axis %>% filter(!str_detect(.$meta.comparison, "control")), aes(y = ymax)) + # I need this to introduce individual ylim values for each facet
#   geom_text(data = plot_labels %>% filter(!str_detect(.$meta.comparison, "control")), aes(x = timepoint, label = condition, color = color),
#             y = 150, size = 5, angle = 90) +
#   scale_color_identity() +
#   geom_text(data = plot_statistics %>% filter(!str_detect(.$meta.comparison, "control")),
#            aes(label = paste("n =", count)),
#            y = 100, x = 6, size = 5)
# 
# 
# 
# ggsave(plot = last_plot(), filename = paste("optimized_for_Paper_", "EdgeR", "pdf", sep = "."), path = here("outs", "identified_clones", "figures"), height = 4, width = 13)
# 
# 
