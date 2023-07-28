
# packages ----------------------------------------------------------------
library(here)
library(tidyverse)
library(lubridate)


# Import and clean data ---------------------------------------------------


#Sys.setlocale("LC_TIME", "en_GB.UTF-8")

data <- read_tsv(here("data", "D11_samples.txt")) %>% 
  mutate(date = dmy(date_of_blood_sampling)) %>%
  mutate(assay_legend = str_remove(.$timepoint, "tp.*")) %>%
  mutate(assay_facets = case_when(assay_legend == "COVID" ~ "VAX_COVID",
                                  assay_legend == "VAX" ~ "VAX_COVID",
                                  TRUE ~ assay_legend)) %>% 
  mutate(assay_facets = factor(assay_facets, levels = c("VAX_COVID", "ELISPOT", "tcrPBMC", "SC", "CFSE")))

dates_vector <- unique(data$date)
dates_labels <- unique(data$days_after_VAX1)

# plot --------------------------------------------------------------------


ggplot(data, aes(x = date, y = 0)) +
  geom_vline(xintercept = dates_vector, linetype = "dashed", color = "grey70", size=0.5) +
  geom_hline(yintercept=0, color = "grey30", size=0.7) +
  geom_point(size=3, stroke = 1.2, aes(color = assay_legend, shape = assay_legend)) +
  geom_point(data = data %>%
               filter(assay_legend == "tcrPBMC" & replica == "rep2"),
             aes(x = date, y = 0, color = assay_legend, shape = assay_legend),
             size = 3, stroke = 1.2,
             position = position_nudge(x = 0, y = 0.5)) +
  scale_color_manual(values = c("#0077bb", "#009988", "#33bbee", "#EE7733", "#ee3377", "darkorchid3", "grey20", "#cc3311")) +
  scale_shape_manual(values = c(0,1,2,4,9,7,10,12)) +
  ylim(-0.3, 1) +
  #geom_text(size = 2, nudge_y = 0.5, angle = 45, alpha = 1) +
  facet_grid(rows = "assay_facets", switch = "y") +
  theme_classic() +
  theme(axis.text.y.left = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(),
        strip.text.y.left = element_text(angle = 0, vjust = 0.2),
        strip.background.y=element_rect(color = NA),
        legend.position = "bottom") +
  ylab("") +
  xlab("") +
  theme(panel.spacing.y=unit(0, "lines")) +
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 0.7, vjust = 0.5)) +
  scale_x_date(breaks = dates_vector, labels = dates_vector, date_labels = "%d %b %y",
               sec.axis = dup_axis(labels = dates_labels))

ggsave(here("outs", "plot_timeline.pdf"), width = 10, height = 5)




