# =========================================================
# GET ORGDB
# =========================================================

.get_orgdb <- function(
    organism = c("human", "mouse", "rat"),
    orgdb = NULL
) {
  
  if (!is.null(orgdb)) {
    return(orgdb)
  }
  
  organism <- match.arg(organism)
  
  if (organism == "human") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("Install org.Hs.eg.db first.", call. = FALSE)
    }
    return(org.Hs.eg.db::org.Hs.eg.db)
  }
  
  if (organism == "mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Install org.Mm.eg.db first.", call. = FALSE)
    }
    return(org.Mm.eg.db::org.Mm.eg.db)
  }
  
  if (organism == "rat") {
    if (!requireNamespace("org.Rn.eg.db", quietly = TRUE)) {
      stop("Install org.Rn.eg.db first.", call. = FALSE)
    }
    return(org.Rn.eg.db::org.Rn.eg.db)
  }
}


# =========================================================
# MAP ENSEMBL TO SYMBOL
# =========================================================

map_ensembl_to_feature <- function(
    ensembl_vec,
    orgdb,
    mart = NULL
) {
  
  ens_clean <- sub("\\..*$", "", ensembl_vec)
  mapping <- NULL
  
  if (!is.null(mart)) {
    
    log_msg("Mapping Ensembl -> SYMBOL via biomaRt...")
    
    mapping <- tryCatch({
      
      biomaRt::getBM(
        attributes = c("ensembl_gene_id", "external_gene_name"),
        filters = "ensembl_gene_id",
        values = unique(ens_clean),
        mart = mart
      ) |>
        tibble::as_tibble() |>
        dplyr::rename(
          ENSEMBL = ensembl_gene_id,
          SYMBOL = external_gene_name
        ) |>
        dplyr::mutate(
          SYMBOL = ifelse(SYMBOL == "", NA, SYMBOL)
        )
      
    }, error = function(e) {
      log_msg("biomaRt mapping failed. Falling back to OrgDb.")
      NULL
    })
  }
  
  if (is.null(mapping)) {
    
    log_msg("Mapping Ensembl -> SYMBOL using OrgDb...")
    
    mapping <- AnnotationDbi::select(
      x = orgdb,
      keys = unique(ens_clean),
      keytype = "ENSEMBL",
      columns = c("SYMBOL", "ENTREZID")
    ) |>
      tibble::as_tibble() |>
      dplyr::mutate(
        SYMBOL = ifelse(SYMBOL == "", NA, SYMBOL),
        FEATURE = ifelse(!is.na(SYMBOL), SYMBOL, ENSEMBL)
      ) |>
      dplyr::select(ENSEMBL, FEATURE)
    
  } else {
    
    mapping <- mapping |>
      dplyr::mutate(FEATURE = SYMBOL) |>
      dplyr::select(ENSEMBL, FEATURE)
  }
  
  id_map <- mapping |>
    dplyr::distinct(ENSEMBL, .keep_all = TRUE)
  
  feat <- id_map$FEATURE[
    match(ens_clean, id_map$ENSEMBL)
  ]
  
  feat[is.na(feat)] <- ens_clean[is.na(feat)]
  
  feat
}


# =========================================================
# MAP ENTREZ IDS
# =========================================================

map_entrez_ids <- function(
    res_df,
    orgdb
) {
  
  log_msg("Mapping ENTREZ IDs...")
  
  ens_clean <- sub("\\..*$", "", res_df$ENSEMBL)
  
  entrez_map <- AnnotationDbi::select(
    x = orgdb,
    keys = unique(ens_clean),
    columns = "ENTREZID",
    keytype = "ENSEMBL"
  )
  
  entrez_map <- entrez_map |>
    dplyr::distinct(ENSEMBL, .keep_all = TRUE)
  
  res_df$ENTREZID <- entrez_map$ENTREZID[
    match(ens_clean, entrez_map$ENSEMBL)
  ]
  
  res_df
}

# =========================================================
# ANNOTATE DE RESULTS
# =========================================================
#' Annotate differential expression results
#'
#' @param res_df Differential expression result table.
#' @param orgdb OrgDb annotation database.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param logfc_cutoff Absolute log2 fold-change cutoff.
#'
#' @return Annotated and processed differential expression table.
#'
#' @keywords internal
#' @noRd
annotate_de_results <- function(
    res_df,
    orgdb,
    padj_cutoff = 0.05,
    logfc_cutoff = 0.25
) {
  
  res_df$SYMBOL <- map_ensembl_to_feature(
    ensembl_vec = res_df$ENSEMBL,
    orgdb = orgdb
  )
  
  res_df <- map_entrez_ids(
    res_df = res_df,
    orgdb = orgdb
  )
  
  processed <- process_deg_results(
    res_df = res_df,
    padj_cutoff = padj_cutoff,
    logfc_cutoff = logfc_cutoff
  )
  
  processed$res_df
}