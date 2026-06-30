# ========================================================= #
# EXPORT RESULTS
# ========================================================= #

export_deg_results <- function(
    res_df,
    sig_deg,
    output_dir = "DEGgo_out"
) {
  
  log <- .deggo_msg(verbose = TRUE, prefix = "DEGgo")
  log("Exporting DEG results...", type = "info")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cols <- intersect(
    c(
      "feature",
      "symbol",
      "SYMBOL",
      "ENSEMBL",
      "gene_id",
      "gene_name",
      "entrez_id",
      "ENTREZID",
      "DEG_status",
      "significance",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj",
      "negLog10Padj",
      "contrast",
      "comparison"
    ),
    colnames(res_df)
  )
  
  res_df <- res_df[, c(cols, setdiff(colnames(res_df), cols)), drop = FALSE]
  
  if (!is.null(sig_deg) && is.data.frame(sig_deg) && nrow(sig_deg)) {
    sig_deg <- sig_deg[
      order(sig_deg$padj),
      ,
      drop = FALSE
    ]
    
    sig_deg <- sig_deg[
      ,
      c(cols, setdiff(colnames(sig_deg), cols)),
      drop = FALSE
    ]
  }
  
  utils::write.table(
    res_df,
    file.path(output_dir, "All_DEGs.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  utils::write.table(
    sig_deg,
    file.path(output_dir, "Significant_DEGs.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}


# ============================================================ #
# DEGgo export pipeline
# ============================================================ #

.deggo_front_cols <- function(df) {
  if (is.null(df) || !is.data.frame(df)) {
    return(df)
  }
  
  front <- intersect(
    c(
      "comparison",
      "contrast",
      "feature",
      "symbol",
      "SYMBOL",
      "ENSEMBL",
      "gene_id",
      "gene_name",
      "entrez_id",
      "ENTREZID",
      "DEG_status",
      "significance"
    ),
    colnames(df)
  )
  
  df[, c(front, setdiff(colnames(df), front)), drop = FALSE]
}


.deggo_extract_go_df <- function(go_obj) {
  if (is.null(go_obj)) {
    return(NULL)
  }
  
  if (is.data.frame(go_obj)) {
    return(go_obj)
  }
  
  if (is.list(go_obj)) {
    if (!is.null(go_obj$go_results) && is.data.frame(go_obj$go_results)) {
      return(go_obj$go_results)
    }
    
    if (!is.null(go_obj$results) && is.data.frame(go_obj$results)) {
      return(go_obj$results)
    }
  }
  
  NULL
}


.deggo_export_overall <- function(
    tables,
    output_file,
    comparison_col = "comparison",
    extractor = NULL
) {
  
  if (is.null(tables) || !length(tables)) {
    return(invisible(NULL))
  }
  
  overall <- do.call(
    rbind,
    lapply(names(tables), function(nm) {
      
      x <- tables[[nm]]
      
      if (!is.null(extractor)) {
        x <- extractor(x)
      }
      
      if (is.null(x) || !is.data.frame(x) || !nrow(x)) {
        return(NULL)
      }
      
      x[[comparison_col]] <- nm
      x <- .deggo_front_cols(x)
      
      x
    })
  )
  
  if (!is.null(overall) && nrow(overall)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    
    utils::write.table(
      overall,
      output_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  invisible(overall)
}


#' Export DEG tables
#'
#' @keywords internal
#' @noRd
.deggo_export_results <- function(
    de_results,
    dirs
) {
  
  if (is.null(de_results$results)) {
    return(de_results)
  }
  
  dir.create(dirs$results, recursive = TRUE, showWarnings = FALSE)
  
  for (nm in names(de_results$results)) {
    df <- .deggo_front_cols(de_results$results[[nm]])
    
    utils::write.table(
      df,
      file.path(dirs$results, paste0(nm, ".tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  invisible(de_results)
}


#' Export significant DEG tables
#'
#' @keywords internal
#' @noRd
.deggo_export_significant <- function(
    de_results,
    dirs
) {
  
  if (is.null(de_results$sig_deg)) {
    return(de_results)
  }
  
  dir.create(dirs$significant, recursive = TRUE, showWarnings = FALSE)
  
  for (nm in names(de_results$sig_deg)) {
    df <- .deggo_front_cols(de_results$sig_deg[[nm]])
    
    utils::write.table(
      df,
      file.path(dirs$significant, paste0(nm, "_significant.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  invisible(de_results)
}


#' Export cleaned significant DEG tables
#'
#' @keywords internal
#' @noRd
.deggo_export_clean_significant <- function(
    de_results,
    dirs
) {
  
  if (is.null(de_results$sig_deg_clean)) {
    return(de_results)
  }
  
  dir.create(dirs$significant, recursive = TRUE, showWarnings = FALSE)
  
  for (nm in names(de_results$sig_deg_clean)) {
    df <- .deggo_front_cols(de_results$sig_deg_clean[[nm]])
    
    utils::write.table(
      df,
      file.path(dirs$significant, paste0(nm, "_significant_clean.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  invisible(de_results)
}


#' Export overall DEG tables
#'
#' @keywords internal
#' @noRd
.deggo_export_overall_tables <- function(
    de_results,
    dirs,
    analysis_mode
) {
  
  all_de <- .deggo_export_overall(
    tables = de_results$results,
    output_file = file.path(
      dirs$results,
      paste0(analysis_mode, "_overall_DE.tsv")
    )
  )
  
  all_sig <- .deggo_export_overall(
    tables = de_results$sig_deg,
    output_file = file.path(
      dirs$significant,
      paste0(analysis_mode, "_overall_significant.tsv")
    )
  )
  
  all_sig_clean <- .deggo_export_overall(
    tables = de_results$sig_deg_clean,
    output_file = file.path(
      dirs$significant,
      paste0(analysis_mode, "_overall_significant_clean.tsv")
    )
  )
  
  all_go <- .deggo_export_overall(
    tables = de_results$go_results,
    output_file = file.path(
      dirs$go,
      paste0(analysis_mode, "_overall_GO.tsv")
    ),
    extractor = .deggo_extract_go_df
  )
  
  de_results$overall_de <- all_de
  de_results$overall_significant <- all_sig
  de_results$overall_significant_clean <- all_sig_clean
  de_results$overall_go <- all_go
  
  de_results
}


#' Export DEGgo summary
#'
#' @keywords internal
#' @noRd
.deggo_export_summary <- function(
    de_results,
    output_dir,
    analysis_mode
) {
  
  if (is.null(de_results$summary)) {
    return(de_results)
  }
  
  utils::write.table(
    de_results$summary,
    file.path(output_dir, paste0(analysis_mode, "_summary.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  invisible(de_results)
}


#' Export reproducibility bundle
#'
#' @keywords internal
#' @noRd
.deggo_export_reproducibility <- function(
    de_results,
    repro_dir,
    save_reproducibility = TRUE
) {
  
  if (!isTRUE(save_reproducibility)) {
    return(de_results)
  }
  
  .write_session_info(dirname(repro_dir))
  
  .save_repro(
    res = de_results,
    repro_dir = repro_dir,
    save_reproducibility = TRUE
  )
  
  invisible(de_results)
}


#' Run all DEGgo exports
#'
#' @keywords internal
#' @noRd
.deggo_export_all <- function(
    de_results,
    dirs,
    output_dir,
    analysis_mode,
    repro_dir,
    save_reproducibility = TRUE
) {
  
  de_results <- .deggo_export_results(de_results, dirs)
  
  de_results <- .deggo_export_significant(de_results, dirs)
  
  de_results <- .deggo_export_clean_significant(de_results, dirs)
  
  de_results <- .deggo_export_overall_tables(
    de_results = de_results,
    dirs = dirs,
    analysis_mode = analysis_mode
  )
  
  de_results <- .deggo_export_summary(
    de_results,
    output_dir,
    analysis_mode
  )
  
  de_results <- .deggo_export_reproducibility(
    de_results,
    repro_dir,
    save_reproducibility
  )
  
  invisible(de_results)
}