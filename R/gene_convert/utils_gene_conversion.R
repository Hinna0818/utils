#' Convert Human Gene Symbols to Mouse Gene Symbols via Ortholog Mapping
#'
#' This function maps a vector of human gene symbols to their corresponding mouse gene symbols
#' using NCBI ortholog mapping and annotation databases (`org.Hs.eg.db` and `org.Mm.eg.db`).
#'
#' @param human_symbols A character vector of human gene symbols.
#' @param ortholog_path Path to the ortholog mapping file,
#'        containing columns such as `#tax_id`, `GeneID`, `Other_tax_id`, `Other_GeneID`.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{HumanSymbol}{Original human gene symbol}
#'   \item{HumanEntrez}{NCBI Entrez ID for human gene}
#'   \item{MouseEntrez}{Mapped mouse Entrez ID}
#'   \item{MouseSymbol}{Mapped mouse gene symbol}
#' }
#'
#' @examples
#' \dontrun{
#'   human_genes <- c("IL10", "TNF", "CCL2")
#'   result <- convert_human_to_mouse(human_genes)
#'   print(result)
#' }
#'
#' @export

convert_human_to_mouse <- function(human_symbols, 
                                   ortholog_path = NULL) {

  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(readr)
  
  ortholog <- read_tsv(ortholog_path, show_col_types = FALSE)
  ortholog_clean <- ortholog %>%
    filter(`#tax_id` == 9606, Other_tax_id == 10090) %>%
    select(HumanEntrez = GeneID, MouseEntrez = Other_GeneID) %>%
    mutate(across(everything(), as.character))
  
  human_entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = human_symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  
  human_df <- data.frame(
    HumanSymbol = names(human_entrez),
    HumanEntrez = as.character(human_entrez),
    stringsAsFactors = FALSE
  )
  
  merged <- inner_join(human_df, ortholog_clean, by = "HumanEntrez")
  
  merged$MouseSymbol <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = merged$MouseEntrez,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  
  return(merged)
}

#' Convert Mouse Gene Symbols to Human Gene Symbols via Ortholog Mapping
#'
#' This function maps a vector of mouse gene symbols to their corresponding human gene symbols
#' using NCBI ortholog data and annotation databases (`org.Mm.eg.db` and `org.Hs.eg.db`).
#'
#' @param mouse_symbols A character vector of mouse gene symbols.
#' @param ortholog_path Path to the ortholog mapping file (e.g., from NCBI's `gene_orthologs.gz`), 
#'        containing columns such as `#tax_id`, `GeneID`, `Other_tax_id`, `Other_GeneID`.
#'
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{MouseSymbol}{Original mouse gene symbol}
#'   \item{MouseEntrez}{NCBI Entrez ID for mouse gene}
#'   \item{HumanEntrez}{Mapped human Entrez ID}
#'   \item{HumanSymbol}{Mapped human gene symbol}
#' }
#'
#' @examples
#' \dontrun{
#'   mouse_genes <- c("Il10", "Tnf", "Ccl2")
#'   result <- convert_mouse_to_human(mouse_genes)
#'   print(result)
#' }
#'
#' @export

convert_mouse_to_human <- function(mouse_symbols, 
                                   ortholog_path = NULL) {
  
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(readr)
  
  ortholog <- read_tsv(ortholog_path, show_col_types = FALSE)
  ortholog_clean <- ortholog %>%
    filter(`#tax_id` == 9606, Other_tax_id == 10090) %>%
    dplyr::select(HumanEntrez = GeneID, MouseEntrez = Other_GeneID) %>%
    mutate(across(everything(), as.character))  
  
  
  mouse_entrez <- mapIds(
    org.Mm.eg.db,
    keys = mouse_symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  
  mouse_df <- data.frame(
    MouseSymbol = names(mouse_entrez),
    MouseEntrez = as.character(mouse_entrez),
    stringsAsFactors = FALSE
  )
  
  merged <- mouse_df %>%
    inner_join(ortholog_clean, by = "MouseEntrez")
  
  merged$HumanSymbol <- mapIds(
    org.Hs.eg.db,
    keys = merged$HumanEntrez,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  
  return(merged)
}








