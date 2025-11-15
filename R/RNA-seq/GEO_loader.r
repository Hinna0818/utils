#' download RNA-seq expression data from GEO datasets
#' @param GSE_id Character; An valid query GSE id number.
#' @param get_sample Logical; Whether to download sample information. Default is TRUE.
#' @param transform Logical; Whether to transform raw counts into `log10(counts)`. Default is TRUE.
#' @import org.Hs.eg.db
#' @importFrom data.table fread
#' @importFrom AnnotationDbi select
#' @importFrom stats na.omit
#' @importFrom GEOquery getGEO
#' @importFrom Biobase pData
#' @importFrom edgeR cpm
#' 
#' @return a list containing RNA-seq raw counts, group information(if needed), and log10(counts) (if needed).
getGEOcts <- function(GSE_id, get_sample = TRUE, transform = TRUE){
  if (is.null(GSE_id)){
    stop("Please enter a valid GSE id.")
  }
  
  # load counts table from GEO
  urld <-"https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&"
  gse <- GSE_id
  path <- paste0(urld, "acc=", gse, "&format=file&file=", gse, "_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  file <- paste0(gse,"_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  download.file(url = path, destfile = file)
  
  tbl <- as.matrix(data.table::fread(file, header = T, colClasses = "integer"), rownames = 1)
  ensembl_matrix <- as.data.frame(tbl)
  
  e2s <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = rownames(ensembl_matrix),
                               columns ="SYMBOL",
                               keytype ="ENTREZID")
  ids <- stats::na.omit(e2s)
  ids <- ids[!duplicated(ids$SYMBOL), ]
  ids <- ids[!duplicated(ids$ENTREZID), ]
  symbol_matrix <- ensembl_matrix[match(ids$ENTREZID, rownames(ensembl_matrix)), ]
  rownames(symbol_matrix) <- ids$SYMBOL
  
  if (get_sample){
    gset <- GEOquery::getGEO(gse, destdir ='.', getGPL = F)
    pd <- Biobase::pData(gset[[1]])
    com <- intersect(rownames(pd), colnames(symbol_matrix))
    symbol_matrix <- symbol_matrix[, com]
    pd <- pd[com,]
    #group_list <- pd$`group:ch1` # depends on your data
  }
  
  if (transform){
    dat <- log10(edgeR::cpm(symbol_matrix) + 1)
  }
  
  res <- list(
    symbol_matrix = symbol_matrix,
    group_info = if (transform) pd else NULL,  # or group_list = group_list instead
    dat = if (transform) dat else NULL
  )
  
  return(res)
}
