library(Seurat)
library(tidyverse)
library(patchwork)
library(here)
library(pals)

set.seed(42)

#path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
intgr <- read_rds(path_to_intgr_seurat)

# functions for ggplot based feature plots --------------------------------

gg_feature_plot <- function(dataset_coordinates, seurat_object, assay, gene_name) {
  ifelse(assay == "integrated", 
         plot <- dataset_coordinates %>%
           inner_join(data.frame("gene" = seurat_object@assays$integrated@data[gene_name, ], "barcode" = colnames(seurat_object@assays$integrated@data)), by = "barcode") %>%
           ggplot(aes(x=UMAP_1, y=UMAP_2)) +
           geom_point(aes(color = gene), size = 0.01) +
           theme_void() +
           ggtitle(paste(noquote(assay), gene_name)) + 
           scale_color_gradient2(midpoint = 0, low = "seagreen2", high = "purple", mid = "grey90") +
           coord_fixed(),
         plot <- dataset_coordinates %>%
           inner_join(data.frame("gene" = seurat_object@assays$RNA@data[gene_name, ], "barcode" = colnames(seurat_object@assays$RNA@data)), by = "barcode") %>%
           ggplot(aes(x=UMAP_1, y=UMAP_2)) +
           geom_point(aes(color = gene), size = 0.01) +
           theme_void() +
           theme(legend.position = "none") +
           ggtitle(paste(noquote(assay), gene_name)) +
           scale_color_gradient(low = "grey90", high = "purple") + 
           coord_fixed())
  plot
}



gg_feature_plot_protein <- function(dataset_coordinates, seurat_object, assay, gene_name) {
  plot <- dataset_coordinates %>%
    filter(nCount_ADT > 0) %>%
    inner_join(data.frame("gene" = seurat_object@assays$ADT@data[gene_name, ], "barcode" = colnames(seurat_object@assays$ADT@data)), by = "barcode") %>%
    ggplot(aes(x=UMAP_1, y=UMAP_2)) +
    geom_point(aes(color = gene), size = 0.01) +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle(paste(noquote(assay), gene_name)) + 
    scale_color_gradient(low = "grey90", high = "purple") +
    coord_fixed()
  plot
}


# extract UMAP coordinates for cells --------------------------------------

umap_tx <- intgr@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(intgr[[]]) %>%
  rownames_to_column(var = "barcode") %>%
  slice_sample(prop = 1) #mix rows for more clear plotting



# lists of marker genes expression for supplements -----------------------

# RNA

vec_genes_paper <- c("FGFBP2", "GZMB", "GZMH", "CXCR3", "NKG7",
                     "GNLY", "CCL5", "GZMK", "KLRB1", "CCR6",
                     "LGALS1", "PTGDR2", "CCR4", "CCR10", "FOXP3",
                     "IL2RA", "TIGIT", "CTLA4", "CXCR5",
                     "LAG3", "PDCD1", "MKI67", "TYMS", "SELL", "CCR7", "CD69",
                     "IFI6", "IFI44L", "XAF1",
                     "RTKN2", "HLA-DRB1", "GBP5")

plots <- map(vec_genes_paper, function(x) gg_feature_plot(dataset_coordinates = umap_tx, seurat_object = intgr, assay = "RNA", gene_name = x)) 

ggsave(plot = (wrap_plots(plots) + plot_layout(ncol = 4)), filename = here("outs", "scRNAseq", "figures", "annotation_SI_RNA_umaps.png"), height = 12, width = 8, dpi = 600)


# protein

vec_proteins_paper <- c("AB-OX40L", "AB-CD3", "AB-CD45RA", "AB-CCR4", "AB-CD25",
                        "AB-CD45RO", "AB-PD1", "AB-TIGIT",
                        "AB-PTGDR2", "AB-TCRg-d", "AB-CXCR3", "AB-CCR6",
                        "AB-CD69", "AB-CD62L", "AB-CCR7", "AB-CTLA4",
                        "AB-LAG3", "AB-CD27", "AB-HLA-DR", "AB-CXCR5",
                        "AB-CD1d", "AB-ICOS", "AB-IL4R", "AB-IL7R")

plots <- map(vec_proteins_paper, function(x) gg_feature_plot_protein(dataset_coordinates = umap_tx, seurat_object = intgr, assay = "ADT", gene_name = x))
ggsave(plot = (wrap_plots(plots) + plot_layout(ncol = 4)), filename = here("outs", "scRNAseq", "figures", "annotation_SI_Protein_umaps.png"), height = 14, width = 8, dpi = 600)





# nice UMAP ---------------------------------------------------------------


umap_tx <- umap_tx %>%
  left_join(annotation_table) %>%
  left_join(color_code)

labels_coordinates <- umap_tx %>%
  group_by(Subset) %>%
  summarise(x_position = mean(UMAP_1),
            y_position = mean(UMAP_2))

ggplot() +
  geom_point(data = umap_tx, aes(x = UMAP_1, y = UMAP_2, color = colors), size = 0.2) +
  scale_color_identity(guide = "legend", labels = annotation_table$Subset, breaks = color_code$colors) +
  theme_minimal() +
  guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) +
  #geom_label(data = labels_coordinates, aes(x = x_position, y = y_position, label = Subset), fill = "white", alpha = 0.6) +
  coord_fixed() +
  theme_void() +
  theme(legend.text = element_text(size = 14), 
        legend.title = element_text(color = "white"))

ggsave(here("outs", "scRNAseq", "figures", "UMAP.png"), height = 6, width = 9, dpi = 600)


# supplementary UMAPs - batch effects -----------------------------------------------

color_cluster_annotation <- color_code %>% left_join(annotation_table)

umap_tx <- umap_tx %>%
  left_join(annotation_table) %>%
  left_join(color_code)

labels_coordinates <- umap_tx %>%
  group_by(Subset) %>%
  summarise(x_position = mean(UMAP_1),
            y_position = mean(UMAP_2))

ggplot() +
  geom_point(data = umap_tx, aes(x = UMAP_1, y = UMAP_2, color = colors), size = 0.2) +
  scale_color_identity(guide = "legend", labels = color_cluster_annotation$Subset, breaks = color_cluster_annotation$colors) +
  theme_void() +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  facet_wrap(facets = vars(orig.ident)) +
  coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) +
  #geom_label(data = labels_coordinates, aes(x = x_position, y = y_position, label = Subset), fill = "white", alpha = 0.6) +
  theme(legend.text = element_text(size = 14), 
        legend.title = element_text(color = "white"))

ggsave(here("outs", "scRNAseq", "figures", "UMAPbatches.png"), height = 9, width = 9, dpi = 600)





