
# packages ----------------------------------------------------------------
library(here)
library(patchwork)
library(Seurat)
library(tidyverse)
library(pals)

set.seed(42)


# data input -----------------------------------------------------------

intgr <- read_rds(path_to_intgr_seurat)


query <- map(path_to_query, function(x) read_rds(x)) # upload query files for reference transfer procedure
names(query) <- c("D11_COVID1", "D11_COVID2_rep1", "D11_COVID2_rep2")



# Seurat reference mapping ------------------------------------------------


anchors <- map(query, function(x) {
  FindTransferAnchors(
  reference = intgr,
  query = x,
  reference.reduction = "pca",
  k.anchor = 5)
})

write_rds(anchors, paste0(here("outs", "scRNAseq", "anchors"), date(), ".rds"))

query_scRNA2 <- map2(query, anchors, function(x, y) MapQuery(anchorset = y, reference = intgr, query = x,
                                                             refdata = "Subset", reference.reduction = "pca", reduction.model = "umap"))

write_rds(query_scRNA2, here("outs", "scRNAseq", "query_scRNA2.rds"))


umap_tx_intgr <- intgr@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = "barcode") %>%
  left_join(intgr[[]] %>% rownames_to_column(var = "barcode"), by = "barcode") %>%
  left_join(annotation_table) %>%
  left_join(color_code)

umap_tx_D11_COVID1 <- query_scRNA2$D11_COVID1@reductions$ref.umap@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = "barcode") %>%
  left_join(query_scRNA2$D11_COVID1[[]] %>% rownames_to_column(var = "barcode"), by = "barcode") %>%
  left_join(annotation_table, by = c("predicted.id" = "Subset")) %>%
  left_join(color_code)

umap_tx_D11_COVID2_rep1 <- query_scRNA2$D11_COVID2_rep1@reductions$ref.umap@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = "barcode") %>%
  left_join(query_scRNA2$D11_COVID2_rep1[[]] %>% rownames_to_column(var = "barcode"), by = "barcode") %>%
  left_join(annotation_table, by = c("predicted.id" = "Subset")) %>%
  left_join(color_code)

umap_tx_D11_COVID2_rep2 <- query_scRNA2$D11_COVID2_rep2@reductions$ref.umap@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = "barcode") %>%
  left_join(query_scRNA2$D11_COVID2_rep2[[]] %>% rownames_to_column(var = "barcode"), by = "barcode") %>%
  left_join(annotation_table, by = c("predicted.id" = "Subset")) %>%
  left_join(color_code)

ggplot() +
  geom_point(data = umap_tx_intgr, aes(x = UMAP_1, y = UMAP_2), size = 0.01, color = "grey96") +
  geom_point(data = umap_tx_D11_COVID1, aes(x = refUMAP_1, y = refUMAP_2, color = colors), size = 0.01, alpha = 1) +
  scale_color_identity(guide = "legend", labels = annotation_table$Subset, breaks = color_code$colors) +
  theme_void() +
  ggtitle("D11_COVID1") +
  coord_fixed() +
  theme(legend.position = "none") -> p1


ggplot() +
  geom_point(data = umap_tx_intgr, aes(x = UMAP_1, y = UMAP_2), size = 0.01, color = "grey96") +
  geom_point(data = umap_tx_D11_COVID2_rep1, aes(x = refUMAP_1, y = refUMAP_2, color = colors), size = 0.01, alpha = 1) +
  scale_color_identity(guide = "legend", labels = annotation_table$Subset, breaks = color_code$colors) +
  theme_void() +
  ggtitle("D11_COVID2_rep1") +
  coord_fixed() +
  theme(legend.position = "none")  -> p2

ggplot() +
  geom_point(data = umap_tx_intgr, aes(x = UMAP_1, y = UMAP_2), size = 0.01, color = "grey96") +
  geom_point(data = umap_tx_D11_COVID2_rep2, aes(x = refUMAP_1, y = refUMAP_2, color = colors), size = 0.01, alpha = 1) +
  scale_color_identity(guide = "legend", labels = annotation_table$Subset, breaks = color_code$colors) +
  theme_void() +
  ggtitle("D11_COVID2_rep2") +
  coord_fixed() +
  theme(legend.position = "none")  -> p3

wrap_plots(list(p1, p2, p3)) + plot_layout(ncol = 1)
ggsave(here("outs", "scRNAseq", "figures", "labelTransfer.png"), height = 6, width = 3, dpi = 600)







