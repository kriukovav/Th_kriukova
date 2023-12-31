library(Seurat)
library(tidyverse)
library(here)
library(patchwork)
library(pals)

set.seed(42)

# path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
# folder_with_gating_results <- here("data", "gating_V2")"

dir.create(here("outs", "scRNAseq", "figures"), recursive = T)

# import data --------------------------------------------------

intgr <- read_rds(path_to_intgr_seurat) # read seurat object (fully prepared by Daniil Lukyanov)


# downsample data for sortseq -------------------------------------

# subset 20000 cells 
# for sortseq we only include cells from this study
intgr_for_sortseq <- subset(intgr, subset = orig.ident %in% c("DL01_rep1", "DL01_rep2", "DL04_rep1", "DL05_rep1"))
intgr_for_sortseq <- subset(intgr_for_sortseq, cells = sample(1:length(rownames(intgr_for_sortseq[[]])), 20000))

# dot plot sortseq -------------------------------------------------------------

intgr_for_sortseq[[]] %>%
  rownames_to_column(var = "barcode") %>%
  left_join(annotation_table) %>%
  left_join(color_code) %>%
  mutate(number = as.factor(as.numeric(.$number))) -> plot_data

plot_data %>%
  pivot_longer(cols = starts_with("sortseq") & ends_with("trb"), names_to = "Sort_seq_phenotype", values_to = "freq") %>%
  mutate(Sort_seq_phenotype = str_remove_all(.$Sort_seq_phenotype, "sortseq\\_|\\_trb")) %>%
  mutate(Sort_seq_phenotype = str_replace(.$Sort_seq_phenotype, "\\.", "\\-")) -> plot_data_long

plot_data_long %>% 
  mutate(freq = case_when(is.na(.$freq) ~ 0,
                          TRUE ~ 1)) -> plot_data_long


plot_data_long %>%
  select(barcode, Sort_seq_phenotype, freq, number) %>%
  split(f = .$Sort_seq_phenotype) %>%
  map(function(x) {
    x %>%
      group_by(number) %>%
      summarise(n_cells_in_cluster = n(),
                Th_positive_cells = sum(freq),
                `% of cluster stained by Sort-Seq` = mean(freq))
  }) %>%
  bind_rows(.id = "Sort_seq_phenotype") -> simple_heatmap_v1


plot_data_long %>%
  select(barcode, Sort_seq_phenotype, freq) %>%
  distinct() %>%
  group_by(Sort_seq_phenotype) %>%
  summarize(size_of_sort_seq_subset = sum(freq)) -> size_of_sort_seq_subset_table

plot_data_long %>%
  select(barcode, number) %>%
  distinct() %>%
  group_by(number) %>%
  summarize(size_of_seurat_cluster = n()) -> size_of_seurat_cluster

plot_data_long %>%
  select(barcode, number, Sort_seq_phenotype, freq) %>%
  left_join(size_of_sort_seq_subset_table) %>% 
  left_join(size_of_seurat_cluster) %>%
  group_by(number, Sort_seq_phenotype, size_of_sort_seq_subset, size_of_seurat_cluster) %>%
  summarize(n_cell_Th_positive_within_cluster = sum(freq)) %>% 
  mutate(`% of Sort-Seq clonotypes mapped to cluster` = n_cell_Th_positive_within_cluster / size_of_sort_seq_subset) %>%
  mutate(x_axis = paste0(number, " (", size_of_seurat_cluster, " cells)"),
         y_axis = paste0(Sort_seq_phenotype, " (", size_of_sort_seq_subset, " cells)")) -> simple_heatmap_v2

ink_stain_plot_data <- simple_heatmap_v1 %>%
  full_join(simple_heatmap_v2, by = c("Sort_seq_phenotype", "number"), suffix = c("_v1", "_v2")) %>%
  mutate(x_axis = factor(x_axis, levels = paste0(.$number, " (", .$size_of_seurat_cluster, " cells)") %>% unique()))

y_order <- data.frame("Sort_seq_phenotype" = c("Tfh", "TREG", "Th2a", "Th2", "Th22", "Th17", "Th1-17", "Th1"),
                      "y_order" = 1:8)

color_max <- ink_stain_plot_data %>% 
  left_join(y_order) %>%
  select(`% of cluster stained by Sort-Seq`) %>%
  pull %>%
  max()

dot_plot_sortseq <- ink_stain_plot_data %>%
  left_join(y_order) %>%
  ggplot() +
  geom_point(aes(x = x_axis, y = fct_reorder(y_axis, y_order), color = `% of cluster stained by Sort-Seq`, size = `% of Sort-Seq clonotypes mapped to cluster`)) +
  scale_color_gradient(low = "white", high = "black", limits = c(0, color_max)) +
  scale_size_continuous(limits = c(0.000001, 1)) + #I introduced limits to remove dots corresponding to 0 from the plot
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.title = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("clusters") +
  ylab("Sort-Seq phenotypes") +
  coord_fixed()



# downsample data for cite-seq -------------------------------------

# subset 20000 cells for easier visualization
# for cite-seq heatmap we only take cells from public data
intgr_for_adt <- subset(intgr, subset = orig.ident %in% c("Ncl_EMTAB10026", "Sanger_EMTAB10026", "Cambridge_EMTAB10026"))
intgr_for_adt <- subset(intgr_for_adt, cells = sample(1:length(rownames(intgr_for_adt[[]])), 20000))

# dot plot citeseq --------------------------------------------------------

intgr_for_adt[[]] %>%
  rownames_to_column(var = "barcode") %>%
  left_join(annotation_table) %>%
  left_join(color_code) %>%
  mutate(number = as.factor(as.numeric(.$number))) -> plot_data

filenames <- list.files(folder_with_gating_results, pattern = "citeseq\\.txt", full.names = T)
filenames_short <- list.files(folder_with_gating_results, pattern = "citeseq\\.txt", full.names = F) %>% str_remove("\\_citeseq\\.txt") %>%
  str_replace("^t", "T") %>% str_replace("\\_", "-")

citeseq_gated_th <- map(filenames, read_tsv, col_names = "barcode") # reading files with results of in silico flow-cytometry-like gating
names(citeseq_gated_th) <- filenames_short

citeseq_gated_th <- bind_rows(citeseq_gated_th, .id = "phenotype") %>%
  mutate(freq = 1) %>%
  pivot_wider(names_from = phenotype, values_from = freq, values_fill = 0) %>%
  pivot_longer(names_to = "phenotype", values_to = "freq", cols = starts_with("T"))

data.frame("phenotype" = citeseq_gated_th$phenotype %>% unique()) %>%
  merge(plot_data) -> plot_data_long

plot_data_long %>%
  left_join(citeseq_gated_th, by = c("barcode", "phenotype")) %>%
  select(barcode, number, phenotype, freq) %>%
  as_tibble() %>%
  mutate(freq = case_when(is.na(.$freq) ~ 0,
                          TRUE ~ .$freq)) %>%
  split(f = .$phenotype) %>%
  map(function(x) {
    x %>%
      group_by(number) %>%
      summarise(n_cells_in_cluster = n(),
                Th_positive_cells = sum(freq),
                `% of cluster stained by CITE-Seq` = mean(freq))
  }) %>%
  bind_rows(.id = "phenotype") -> simple_heatmap_v1

plot_data_long %>%
  left_join(citeseq_gated_th, by = c("barcode", "phenotype")) %>%
  select(barcode, phenotype, freq) %>%
  mutate(freq = case_when(is.na(.$freq) ~ 0,
                          TRUE ~ .$freq)) %>%
  group_by(phenotype) %>%
  summarize(size_of_cite_seq_subset = sum(freq)) -> size_of_cite_seq_subset_table

plot_data_long %>%
  select(barcode, number) %>%
  distinct() %>%
  group_by(number) %>%
  summarize(size_of_seurat_cluster = n()) -> size_of_seurat_cluster


plot_data_long %>%
  left_join(citeseq_gated_th, by = c("barcode", "phenotype")) %>%
  select(barcode, number, phenotype, freq) %>%
  left_join(size_of_cite_seq_subset_table) %>% 
  left_join(size_of_seurat_cluster) %>%
  mutate(freq = case_when(is.na(.$freq) ~ 0,
                          TRUE ~ .$freq)) %>%
  group_by(number, phenotype, size_of_cite_seq_subset, size_of_seurat_cluster) %>%
  summarize(n_cell_Th_positive_within_cluster = sum(freq)) %>% 
  mutate(`% of CITE-Seq clonotypes mapped to cluster` = n_cell_Th_positive_within_cluster / size_of_cite_seq_subset) %>%
  mutate(x_axis = paste0(number, " (", size_of_seurat_cluster, " cells)"),
         y_axis = paste0(phenotype, " (", size_of_cite_seq_subset, " cells)")) -> simple_heatmap_v2

ink_stain_plot_data <- simple_heatmap_v1 %>%
  full_join(simple_heatmap_v2, by = c("number", "phenotype"), suffix = c("_v1", "_v2")) %>%
  mutate(x_axis = factor(x_axis, levels = paste0(.$number, " (", .$size_of_seurat_cluster, " cells)") %>% unique()))

y_order <- data.frame("phenotype" = c("Tfh", "Treg", "Th2a", "Th2", "Th22", "Th17", "Th1-17", "Th1"),
                      "y_order" = 1:8)

color_max <- ink_stain_plot_data %>% 
  left_join(y_order) %>%
  select(`% of cluster stained by CITE-Seq`) %>%
  pull %>%
  max()


dot_plot_citeseq <- ink_stain_plot_data %>% 
  left_join(y_order) %>% 
  ggplot() +
  geom_point(aes(x = x_axis, y = fct_reorder(y_axis, y_order), color = `% of cluster stained by CITE-Seq`, size = `% of CITE-Seq clonotypes mapped to cluster`)) +
  scale_color_gradient(low = "white", high = "black", limits = c(0, color_max)) +
  scale_size_continuous(limits = c(0.000001, 1)) + #I introduced limits to remove dots corresponding to 0 from the plot
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.title = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("clusters") +
  ylab("CITE-Seq phenotypes") +
  coord_fixed()

wrap_plots(list(dot_plot_sortseq, dot_plot_citeseq)) + plot_layout(nrow = 1, heights = c(1, 1))

ggsave(here("outs", "scRNAseq", "figures", "dot_plot.pdf"), height = 7, width = 13, bg = "white")


