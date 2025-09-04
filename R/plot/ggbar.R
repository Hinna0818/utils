#' @title Bar Plot of Gene Expression with Error Bars and Jitter
#' 
#' @param plot_data A data.frame containing sample-level data. Must include columns: 
#'   \code{group} (factor), \code{tpm} (numeric).
#' @param summary_data A data.frame containing group-level summary statistics. Must include columns: 
#'   \code{group}, \code{mean_tpm}, and \code{sd_tpm}.
#' @param gene_name A character string specifying the name of the gene being plotted. Used as plot title.
#' 
#' @import ggplot2
#' @import dplyr
#' @import ggsci
#' @import ggthemes
#' @return A \code{ggplot2} object.
#' 
ggbar <- function(
  plot_data, 
  summary_data,
  gene_name
){
  
  library(dplyr)
  library(ggplot2)
  library(ggsci)
  library(ggthemes)

  p <- ggplot(data = summary_data, aes(x = .data[["group"]], 
    y = .data[["mean_tpm"]], 
    fill = .data[["group"]])) +
    geom_bar(stat = "identity", width = 0.5, color = "black") +
    geom_errorbar(
      aes(ymin = .data[["mean_tpm"]], 
          ymax = .data[["mean_tpm"]] + .data[["sd_tpm"]]),
      width = 0.2,
      color = "black",
      linewidth = 0.8
    ) +
    geom_jitter(
      data = plot_data,
      aes(y = .data[["tpm"]]),
      width = 0.15,
      size = 2.0,
      alpha = 0.6,
      shape = 21,
      fill = "white",
      stroke = 1
    ) +
    labs(
      title = paste(gene_name, "Expression"),
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

  return(p)
}