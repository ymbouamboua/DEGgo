# ========================================================= #
# EXTRACT GENES FOUND IN GO TERMS
# ========================================================= #
#' Extract genes found in GO terms
#'
#' @param go_results GO enrichment result data.frame.
#' @param genes Character vector of genes to search.
#' @param gene_sep Separator used in the geneID column.
#'
#' @return A data.frame with GO terms containing selected genes.
#' @export
deggo_extract_go_genes <- function(go_results, genes, gene_sep = "/") {
  
  if (is.null(go_results) || !is.data.frame(go_results) || nrow(go_results) == 0) {
    return(data.frame())
  }
  
  required <- c("ID", "Description", "geneID")
  missing <- setdiff(required, colnames(go_results))
  if (length(missing) > 0) {
    stop("Missing columns in go_results: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  
  out <- lapply(seq_len(nrow(go_results)), function(i) {
    
    gene_vec <- strsplit(
      as.character(go_results$geneID[i]),
      gene_sep,
      fixed = TRUE
    )[[1]]
    
    hits <- intersect(gene_vec, genes)
    if (length(hits) == 0) return(NULL)
    
    data.frame(
      GO_ID = go_results$ID[i],
      Description = go_results$Description[i],
      Gene = paste(hits, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, out)
  if (is.null(out)) return(data.frame())
  
  rownames(out) <- NULL
  out
}



# ========================================================= #
# EXTRACT GO GENES PAIRWISE
# ========================================================= #
#' Extract selected genes from pairwise GO results
#'
#' @param results DEGgo result object.
#' @param genes Character vector of genes.
#'
#' @return A data.frame.
#' @export
deggo_extract_go_genes_pairwise <- function(results, genes) {
  
  if (is.null(results$go_results)) {
    stop("results must contain results$go_results.", call. = FALSE)
  }
  
  out <- lapply(names(results$go_results), function(comp) {
    
    go_df <- results$go_results[[comp]]$go_results
    x <- deggo_extract_go_genes(go_df, genes)
    
    if (nrow(x) == 0) return(NULL)
    
    x$comparison <- comp
    x
  })
  
  out <- do.call(rbind, out)
  if (is.null(out)) return(data.frame())
  
  rownames(out) <- NULL
  out
}



# ========================================================= #
# EXTRACT GO KEYWORDS
# ========================================================= #
#' Extract GO terms matching keywords
#'
#' @param results DEGgo result object.
#' @param keywords Character vector of keywords.
#'
#' @return A data.frame.
#' @export
deggo_extract_go_keywords <- function(results, keywords) {
  
  if (is.null(results$go_results)) {
    stop("results must contain results$go_results.", call. = FALSE)
  }
  
  out <- lapply(names(results$go_results), function(comp) {
    
    go_df <- results$go_results[[comp]]$go_results
    
    if (is.null(go_df) || !is.data.frame(go_df) || nrow(go_df) == 0) {
      return(NULL)
    }
    
    idx <- Reduce(
      "|",
      lapply(
        keywords,
        function(k) grepl(k, go_df$Description, ignore.case = TRUE)
      )
    )
    
    if (!any(idx)) return(NULL)
    
    x <- go_df[idx, , drop = FALSE]
    x$comparison <- comp
    x
  })
  
  out <- do.call(rbind, out)
  if (is.null(out)) return(data.frame())
  
  rownames(out) <- NULL
  out
}



# ========================================================= #
# EXTRACT DEG GENES
# ========================================================= #
#' Extract selected DEG genes from DEGgo results
#'
#' @param results DEGgo result object.
#' @param genes Character vector of gene symbols.
#'
#' @return A data.frame.
#' @export
deggo_extract_deg_genes <- function(results, genes) {
  
  if (is.null(results$results)) {
    stop("results must contain results$results.", call. = FALSE)
  }
  
  out <- lapply(names(results$results), function(comp) {
    
    x <- results$results[[comp]]
    
    if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) {
      return(NULL)
    }
    
    required <- c("SYMBOL", "log2FoldChange", "pvalue", "padj")
    missing <- setdiff(required, colnames(x))
    if (length(missing) > 0) {
      stop("Missing columns in results$results[[", comp, "]]: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    
    keep <- x$SYMBOL %in% genes
    if (!any(keep)) return(NULL)
    
    cols <- intersect(
      c("SYMBOL", "log2FoldChange", "pvalue", "padj", "significance"),
      colnames(x)
    )
    
    y <- x[keep, cols, drop = FALSE]
    y$comparison <- comp
    y
  })
  
  out <- do.call(rbind, out)
  if (is.null(out)) return(data.frame())
  
  rownames(out) <- NULL
  out
}



# ========================================================= #
# EXTRACT RHYTHMIC GENES
# ========================================================= #
#' Extract selected genes from DEGgo rhythmicity results
#'
#' @param results A `deggo_rhythm_results` object returned by
#'   [run_deggo_rhythmicity()].
#' @param genes Character vector of gene identifiers.
#'
#' @return A data.frame.
#' @export
deggo_extract_rhythmic_genes <- function(results, genes) {

  if (is.null(results$summary)) {
    stop("results must contain results$summary.", call. = FALSE)
  }

  summary <- results$summary

  if (!is.data.frame(summary) || nrow(summary) == 0) {
    return(data.frame())
  }

  keep <- summary$gene %in% genes
  if (!any(keep)) return(data.frame())

  out <- summary[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}
