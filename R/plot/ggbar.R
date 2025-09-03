## ggbar with error bar
ggbar <- function(
  plot_data, summary_data
){
  p <- ggplot(data = summary_data, aes(x = .data[[group]], 
    y = .data[[mean_tpm]], 
    fill = .data[[group]])) +
    geom_bar(stat = "identity", width = 0.5, color = "black") +
    geom_errorbar(
      aes(ymin = .data[[mean_tpm]], 
          ymax = .data[[mean_tpm]] + .data[[sd_tpm]]),
      width = 0.2,
      color = "black",
      linewidth = 0.8
    ) +
    geom_jitter(
      data = plot_data,
      aes(y = .data[[tpm]]),
      width = 0.15,
      size = 2.0,
      alpha = 0.6,
      shape = 21,
      fill = "white",
      stroke = 1
    ) +
    labs(
      title = paste(gene, "Expression"),
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
      panel.border = element_rect(color = "black", size = 1.5, linetype = "solid", fill = NA)
    ) +
    scale_fill_npg() +
    theme(aspect.ratio = 1)

  return(p)
}