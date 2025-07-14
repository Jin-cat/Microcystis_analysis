library(ggplot2)
library(readr)
library(dplyr)

setwd("C:\\Users\\Desktop\\Microcystis")

df <- read_delim("BGC_group_G.txt", delim = "\t")
df$`Similarity Confidence` <- factor(df$`Similarity Confidence`, levels = c("High", "Medium", "Low"))
df <- df %>% rename(Strain = ...1)

strain_order <- unique(df$Strain)
df$Strain <- factor(df$Strain, levels = strain_order)

df$Type <- factor(df$Type, levels = c("Betalactone", "Cyanobactin", "Microviridin", "NRPS", "PKS", "PKS-NRPS"))

cluster_order <- df %>%
  distinct(Type, Cluster) %>%
  arrange(Type, Cluster) %>%
  pull(Cluster)
df$Cluster <- factor(df$Cluster, levels = cluster_order)

ggplot(df, aes(x = Cluster, y = Strain, color = Type, size = `Similarity Confidence`)) +
  geom_point() +
  scale_color_manual(values = c(
    "Betalactone" = "#fe3939", "Cyanobactin" = "#fc8302",
    "Microviridin" = "#19d62f", "NRPS" = "#28bbe0",
    "PKS" = "#ab79f0", "PKS-NRPS" = "#c87f48")) +
  scale_size_manual(values = c("High" = 5, "Medium" = 3.5, "Low" = 2.5)) +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_y_discrete(limits = rev(levels(df$Strain))) +
  theme_bw() +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major = element_line(color = "gray90")) +
  labs(x = NULL, y = "Strain", color = "BGC Type", size = "Confidence")
