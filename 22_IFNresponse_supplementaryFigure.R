library(tidyverse)
library(patchwork)
library(here)
library(Seurat)


intgr <- read_rds(path_to_intgr_seurat)
intgr <- subset(intgr, subset = orig.ident %in% c("Ncl_EMTAB10026", "Sanger_EMTAB10026", "Cambridge_EMTAB10026"))

umap_tx <- intgr@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(intgr[[]]) %>%
  rownames_to_column(var = "barcode") %>%
  left_join(annotation_table) %>%
  left_join(color_code)

umap_tx %>%
  split(f = .$Status_on_day_collection_summary) %>%
  imap(function(x, y) {
    ggplot() +
      geom_point(data = x, aes(x = UMAP_1, y = UMAP_2, color = colors), size = 0.01) +
      scale_color_identity(guide = "legend", labels = annotation_table$Subset, breaks = color_code$colors) +
      theme_minimal() +
      guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) +
      #geom_label(data = labels_coordinates, aes(x = x_position, y = y_position, label = Subset), fill = "white", alpha = 0.6) +
      coord_fixed() +
      theme_void() +
      theme(legend.position = "none") +
      ggtitle(paste(y))
  }) -> status

status <- status[c("Healthy", "Asymptomatic", "Mild", "Moderate", "Severe", "Critical", "LPS_90mins", "LPS_10hours")]

wrap_plots(status)

ggsave(here("outs", "scRNAseq", "figures", "IFNumaps.png"), dpi = 600, bg = "white", height = 6, width = 5)



umap_tx %>%
  group_by(sample_id, Status_on_day_collection_summary, Subset) %>%
  summarise(n_cells = n()) %>%
  ungroup() %>%
  group_by(sample_id) %>%
  mutate(proportion = n_cells / sum(n_cells)) -> data

x_order <- data.frame("Status_on_day_collection_summary" = c("Healthy", "Asymptomatic", "Mild", "Moderate", "Severe", "Critical", "LPS_90mins", "LPS_10hours"),
                      "x_axis" = c(1:8))

data %>%
  filter(Subset %in% c("EffMem_IFN_response", "Naive_IFN_response")) %>%
  inner_join(x_order) %>% 
  ggplot() +
  geom_jitter(aes(x = fct_reorder(Status_on_day_collection_summary, x_axis), y = proportion), size = 0.1, width = 0.1) +
  stat_summary(fun = median, aes(x = fct_reorder(Status_on_day_collection_summary, x_axis), y = proportion), color = "red", shape = 95, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid.major.y = element_line(size = 0.05, linetype = 'solid', colour = "grey90"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        strip.text = element_text(size = 8)) +
  ylab("Normalized cluster size") +
  xlab("Patient status on the day of sample collection") +
  facet_wrap(facets = "Subset")

ggsave(here("outs", "scRNAseq", "figures", "IFNscatterplots.pdf"), dpi = 600, bg = "white", height = 3, width = 5)



# apply stat test ---------------------------------------------------------

data %>%
  filter(Subset %in% c("EffMem_IFN_response")) -> data_aov

aov(data_aov$proportion ~ data_aov$Status_on_day_collection_summary) -> model_aov

summary(model_aov) %>% capture.output() %>% write_lines(here("outs", "anova_summary_Effmem_IFN_response.txt"))
TukeyHSD(model_aov) %>% broom::tidy() %>% filter(adj.p.value < 0.05) %>% write_tsv(here("outs", "TukeyHSD_Effmem_IFN_response.txt")) 
  


data %>%
  filter(Subset %in% c("Naive_IFN_response")) -> data_aov

aov(data_aov$proportion ~ data_aov$Status_on_day_collection_summary) -> model_aov

summary(model_aov) %>% capture.output() %>% write_lines(here("outs", "anova_summary_Naive_IFN_response.txt"))
TukeyHSD(model_aov) %>% broom::tidy() %>% filter(adj.p.value < 0.05) %>% write_tsv(here("outs", "TukeyHSD_Naive_IFN_response.txt")) 


