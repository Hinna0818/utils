setwd("/Users/hinna/work/utils/R/plot")
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggthemes)

rm(list = ls())
tpm <- openxlsx::read.xlsx("/Users/hinna/work/utils/data/multi_gene_rna_tpm.xlsx", rowNames = TRUE)
coldata <- data.frame(
  condition = sub("[-_.]\\d+$", "", colnames(tpm)),
  row.names = colnames(tpm)
)
coldata$condition <- factor(coldata$condition,
                            levels = c("MII", "zygote", "2cell", "4cell", "morula", "blastocyst"))


goi <- "Cdx2"
goi_exp <- tpm[goi, ]

plot_data <- data.frame(
  sample = colnames(tpm),
  tpm = as.numeric(goi_exp),
  group = coldata$condition
)

summary_data <- plot_data %>%
  group_by(group) %>%
  summarise(
    mean_tpm = mean(tpm),
    sd_tpm = sd(tpm)
  )
save(plot_data, summary_data, file = "/Users/hinna/work/utils/data/Zscan4c_plotdata.rda")


## 可视化
ggplot(data = summary_data, aes(x = group, y = mean_tpm, fill = group)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_errorbar(
    aes(ymin = mean_tpm, ymax = mean_tpm + sd_tpm),
    width = 0.2,
    color = "black",
    linewidth = 0.8
  ) +
  geom_jitter(
    data = plot_data,
    aes(y = tpm),
    width = 0.15,
    size = 2.0,
    alpha = 0.6,
    shape = 21,
    fill = "white",
    stroke = 1
  ) +
  labs(
    title = paste(goi,"Expression"),
    x = NULL,
    y = "Expression (TPM)"
  ) +
  theme_few(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", linewidth = 1.5, linetype = "solid", fill = NA)
  ) +
  scale_fill_npg() +
  theme(aspect.ratio = 1)

source("./ggbar.R")
p1 <- ggbar(plot_data, summary_data, gene_name = "Cdx2")
p1
