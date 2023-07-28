library(Seurat)
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)
library(scatterpie)
library(pals)

set.seed(42)

# data input - reference and query scRNAseq ----------------------------------

intgr <- read_rds(path_to_intgr_seurat)

umap_tx_intgr <- intgr@reductions$umap@cell.embeddings %>%
    as.data.frame() %>%
    rownames_to_column(var = "barcode") %>%
    left_join(intgr[[]] %>% rownames_to_column(var = "barcode"), by = "barcode")

rm(intgr)

query <- read_rds(path_to_query)

# data input - collect all the sets of responding clones -------------------

patient_hla <- read_lines(here("data", "patient_hla.txt"))
pattern_hla <-  patient_hla %>% str_replace("\\-", "\\\\-") %>% str_replace("\\*", "\\\\*") %>% str_replace("\\:", "\\\\:") %>% paste(collapse = "|")


read_tsv(here("outs", "identified_clones", "filtered_clones", "pogorelyy2022.cdr3aa-v-j-1mm.txt")) %>% 
  filter(str_detect(.$meta.best_HLA_assoc, pattern = pattern_hla)) %>%
  filter((PBMCtp01_1 + PBMCtp01_2)/2 >= 3 |
           (PBMCtp03_1 + PBMCtp03_2)/2 >= 3 |
           (PBMCtp08_1 + PBMCtp08_2)/2 >= 3) %>% #filter only clones having at least 3 UMIs in mean_count in vax / covid1 /covid2 timepoints
  select(cdr3aa, v, j) %>% distinct() -> pogorelyy2022

read_tsv(here("outs", "identified_clones", "filtered_clones", "TcellAssay.cdr3aa-cdr3nt-v-j.txt")) %>% 
  select(cdr3aa, v, j) %>% distinct -> TcellAssay

read_tsv(here("outs", "identified_clones", "filtered_clones", "vdjdb.cdr3aa-v-j-1mm.txt")) %>%
  filter(meta.species == "HomoSapiens") %>%
  filter(meta.antigen.species == "SARS-CoV-2") %>% 
  filter(str_detect(.$`meta.mhc.a`, pattern = pattern_hla) | str_detect(.$`meta.mhc.b`, pattern = pattern_hla)) %>% #filter only patient's mhc
  filter((PBMCtp01_1 + PBMCtp01_2)/2 >= 3 |
           (PBMCtp03_1 + PBMCtp03_2)/2 >= 3 |
           (PBMCtp08_1 + PBMCtp08_2)/2 >= 3) %>% #filter only clones having at least 3 UMIs in mean_count in vax / covid1 /covid2 timepoints
  filter(cdr3aa != "CASSETSGSTDTQYF") %>% # filter out irrelevant clone (selected based on very high frequency even before covid at tp00)
  select(cdr3aa, v, j) %>% distinct -> vdjdb


read_tsv(here("outs", "identified_clones", "filtered_clones", "edgeRhits.txt")) %>% 
  split(f = .$meta.comparison) %>%
  map(function(x) select(x, cdr3aa, v, j) %>% distinct) -> list_of_clones1

list_of_clones1 <- list_of_clones1[!str_detect(names(list_of_clones1), pattern = "control")]

list_of_clones2 <- list("pogorelyy2022" = pogorelyy2022, "TcellAssay" = TcellAssay, "vdjdb" = vdjdb)

list_of_clones <- c(list_of_clones1, list_of_clones2)


# Check phenotypes of SARS-CoV-2 specific cells -----------------------

umap_tx <- map(query, function(x) {
  x@reductions$ref.umap@cell.embeddings %>%
    as.data.frame() %>%
    rownames_to_column(var = "barcode") %>%
    left_join(x[[]] %>% rownames_to_column(var = "barcode"), by = "barcode")
}) %>%
  bind_rows()


# check TCRbeta clonotype overlaps between single cell data for D11 and the collection of SARS-CoV-2 specific TCR clonotypes from the same donor
overlaps <- map(list_of_clones, function(x) {
  umap_tx %>%
    filter(TRB_cdr3aa %in% x$cdr3aa) %>%
    filter(TRB_cdr3aa != "CASSETSGSTDTQYF") # filter out irrelevant clone (selected based on very high frequency even before covid at tp00)
}) 


clone_method <- overlaps %>%
  bind_rows(.id = "method") %>%
  select(TRB_cdr3aa, method) %>%
  distinct() %>%
  mutate(count = 1) # introduce count, to calculate, with how many methods was verified each clonotype (needed for scatterpie chart) 

data_scatter_pie <- clone_method %>%
  pivot_wider(id_cols = "TRB_cdr3aa", names_from = method, values_from = count, values_fill = 0)

umap_tx2 <- umap_tx %>% 
  inner_join(data_scatter_pie, by = "TRB_cdr3aa") %>% # inner join for filtering only cells specific to SARS-CoV-2
  mutate(radius = 0.22) # for scatterpie plotting

umap_tx %>%
  ggplot() +
  geom_point(data = umap_tx_intgr, aes(x = UMAP_1, y = UMAP_2), color = "grey96", size = 0.1, alpha = 1) +
  geom_point(data = umap_tx, aes(x = refUMAP_1, y = refUMAP_2), color = "grey96", size = 0.1, alpha = 1) +
  theme_void() +
  geom_scatterpie(aes(x = refUMAP_1, y = refUMAP_2, group = TRB_cdr3aa, r = radius * 1.1), cols = colnames(data_scatter_pie)[2:length(colnames(data_scatter_pie))], data = umap_tx2, color = NA) +
  coord_fixed() +
  scale_fill_manual(values = c( "#00F5A0", "#F600F6", "#4400FA", "darkgreen", "#F6BC00","#F55520")) +
  theme(legend.text = element_text(size = 10)) -> figA


overlaps %>%
  bind_rows(.id = "method") %>%
  ggplot() +
  geom_point(data = umap_tx_intgr, aes(x = UMAP_1, y = UMAP_2), color = "grey96", size = 0.1, alpha = 1) +
  geom_point(data = umap_tx, aes(x = refUMAP_1, y = refUMAP_2), color = "grey96", size = 0.1, alpha = 1) +
  theme_void() +
  geom_point(aes(x = refUMAP_1, y = refUMAP_2, color = TRB_cdr3aa), size = 2) +
  #geom_text_repel(aes(x = refUMAP_1, y = refUMAP_2, label = TRB_cdr3aa), color = "black", size = 5) +
  coord_fixed() +
  theme(legend.direction = "horizontal",
        legend.text = element_text(size = 7)) +
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) -> figB



overlaps %>%
  bind_rows(.id = "method") %>%
  ggplot() +
  geom_point(data = umap_tx_intgr, aes(x = UMAP_1, y = UMAP_2), color = "grey96", size = 0.1, alpha = 1) +
  geom_point(data = umap_tx, aes(x = refUMAP_1, y = refUMAP_2), color = "grey96", size = 0.1, alpha = 1) +
  theme_void() +
  geom_point(data = . %>% left_join(annotation_table, by = c("predicted.id" = "Subset"))  %>% left_join(color_code), aes(x = refUMAP_1, y = refUMAP_2, color = colors), size = 2) +
  scale_color_identity(guide = "legend", labels = annotation_table$Subset, breaks = color_code$colors) +
  #geom_text_repel(aes(x = refUMAP_1, y = refUMAP_2, label = predicted.id), color = "black", size = 4) +
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
  coord_fixed() +
  theme(legend.text = element_text(size = 10)) -> figC


wrap_plots(list(figA, figB, figC)) + plot_layout(nrow = 3, guides = "keep")

ggsave(here("outs", "scRNAseq", "figures", "clones_mapped_threePanels.png"), height = 10, width = 14, dpi = 600)


# draw piechart with summary ----------------------------------------------

overlaps %>%
  bind_rows(.id = "method") %>%
  left_join(annotation_table, by = c("predicted.id" = "Subset")) %>%
  left_join(color_code) %>%
  select(barcode, colors) %>%
  distinct() %>%
  ggplot(aes(x = factor(1), fill = colors)) +
  geom_bar(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_identity(guide = "legend", labels = annotation_table$Subset, breaks = color_code$colors) +
  guides(fill = guide_legend(override.aes = list(size=15, alpha = 1))) +
  theme_void()

ggsave(here("outs", "scRNAseq", "figures", "piechart.pdf"), height = 5, width = 7, dpi = 600)


# suplementary figure

# mapped clones (supplementary UMAPs) -------------------------------------

overlaps %>%
  bind_rows(.id = "method") %>%
  ggplot() +
  geom_point(data = umap_tx_intgr, aes(x = UMAP_1, y = UMAP_2), color = "grey96", size = 0.1, alpha = 1) +
  geom_point(data = umap_tx, aes(x = refUMAP_1, y = refUMAP_2), color = "grey96", size = 0.1, alpha = 1) +
  theme_void() +
  facet_wrap(facets = vars(method)) +
  geom_point(aes(x = refUMAP_1, y = refUMAP_2, fill = TRB_cdr3aa), size = 2.5, shape = 21, color = "white") +
  #scale_fill_identity(guide = "legend", labels = color_code$seurat_clusters, breaks = color_code$colors) +
  #geom_text_repel(aes(x = refUMAP_1, y = refUMAP_2, label = TRB_cdr3aa), color = "black", size = 4) +
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
  coord_fixed() +
  theme(legend.text = element_text(size = 10)) -> p1

overlaps %>%
  bind_rows(.id = "method") %>%
  ggplot() +
  geom_point(data = umap_tx_intgr, aes(x = UMAP_1, y = UMAP_2), color = "grey96", size = 0.1, alpha = 1) +
  geom_point(data = umap_tx, aes(x = refUMAP_1, y = refUMAP_2), color = "grey96", size = 0.1, alpha = 1) +
  theme_void() +
  facet_wrap(facets = vars(method)) +
  geom_point(data = . %>% left_join(annotation_table, by = c("predicted.id" = "Subset")) %>% left_join(color_code), aes(x = refUMAP_1, y = refUMAP_2, fill = colors), shape = 21, color = "white", size = 2.5) +
  scale_fill_identity(guide = "legend", labels = annotation_table$Subset, breaks = color_code$colors) +
  #geom_text_repel(aes(x = refUMAP_1, y = refUMAP_2, label = predicted.id), color = "black", size = 4) +
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
  coord_fixed() +
  theme(legend.text = element_text(size = 10)) -> p2

wrap_plots(list(p1, p2)) + plot_layout(ncol = 1)

ggsave(here("outs", "scRNAseq", "figures", "S6.png"), height = 12, width = 15, dpi = 600)





