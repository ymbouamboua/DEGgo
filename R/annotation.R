# =========================================================
# GET ORGDB
# =========================================================
#' Retrieve an organism-specific OrgDb annotation database
#'
#' Returns a Bioconductor OrgDb annotation database corresponding
#' to the selected organism. If an OrgDb object is supplied,
#' it is returned unchanged.
#'
#' Supported organisms are human, mouse, and rat.
#'
#' @param organism Character string specifying the organism.
#'   One of \code{"human"}, \code{"mouse"}, or \code{"rat"}.
#' @param orgdb Optional OrgDb annotation database. If provided,
#'   this object is returned directly.
#'
#' @return An OrgDb annotation database object.
#'
#' @details
#' The following annotation packages are used:
#' \itemize{
#'   \item \code{org.Hs.eg.db} for human
#'   \item \code{org.Mm.eg.db} for mouse
#'   \item \code{org.Rn.eg.db} for rat
#' }
#'
#' If the required package is not installed, the function
#' stops with an informative error message.
#'
#' This helper is used internally by DEGgo annotation,
#' enrichment, and reporting functions.
#'
#' @seealso
#' \code{\link{map_ensembl_to_feature}},
#' \code{\link{map_entrez_ids}}
#'
#' @keywords internal
#' @noRd
# .get_orgdb <- function(
#     organism = c("human", "mouse", "rat"),
#     orgdb = NULL
# ) {
#
#   if (!is.null(orgdb)) {
#     return(orgdb)
#   }
#
#   organism <- match.arg(organism)
#
#   pkg <- switch(
#     organism,
#     human = "org.Hs.eg.db",
#     mouse = "org.Mm.eg.db",
#     rat   = "org.Rn.eg.db"
#   )
#
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     stop("Install ", pkg, " first.", call. = FALSE)
#   }
#
#   getExportedValue(pkg, pkg)
# }
.get_orgdb <- function(organism, orgdb = NULL) {

  if (!is.null(orgdb)) {
    return(orgdb)
  }

  organism <- match.arg(
    organism,
    choices = c("human", "mouse", "rat", "custom")
  )

  switch(
    organism,

    human = {
      if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        stop(
          "Package 'org.Hs.eg.db' is required.",
          call. = FALSE
        )
      }
      org.Hs.eg.db::org.Hs.eg.db
    },

    mouse = {
      if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
        stop(
          "Package 'org.Mm.eg.db' is required.",
          call. = FALSE
        )
      }
      org.Mm.eg.db::org.Mm.eg.db
    },

    rat = {
      if (!requireNamespace("org.Rn.eg.db", quietly = TRUE)) {
        stop(
          "Package 'org.Rn.eg.db' is required.",
          call. = FALSE
        )
      }
      org.Rn.eg.db::org.Rn.eg.db
    },

    custom = {
      stop(
        paste0(
          "For organism = 'custom', a valid OrgDb object must be supplied via `orgdb`.\n\n",

          "Example:\n",
          "library(org.Custom.eg.db)\n\n",

          "results <- run_deggo(\n",
          "  counts = counts,\n",
          "  metadata = metadata,\n",
          "  organism = 'custom',\n",
          "  orgdb = org.Custom.eg.db::org.Custom.eg.db\n",
          ")\n\n",

          "Custom OrgDb packages can be generated with:\n",
          "  AnnotationForge::makeOrgPackageFromNCBI()\n",
          "or\n",
          "  AnnotationForge::makeOrgPackage()\n"
        ),
        call. = FALSE
      )
    }
  )
}



#' Map Ensembl identifiers to feature names
#'
#' Converts Ensembl gene identifiers into human-readable feature names,
#' preferentially using gene symbols and falling back to Ensembl IDs
#' when annotation is unavailable.
#'
#' Mapping can be performed using either a supplied
#' \code{biomaRt::Mart} object or a Bioconductor OrgDb annotation
#' database.
#'
#' Version suffixes in Ensembl identifiers
#' (e.g. \code{ENSG000001234.5}) are automatically removed prior to
#' annotation.
#'
#' @param ensembl_vec Character vector of Ensembl gene identifiers.
#' @param orgdb OrgDb annotation database (e.g.
#'   \code{org.Hs.eg.db}, \code{org.Mm.eg.db}).
#' @param mart Optional \code{biomaRt} Mart object. If provided,
#'   annotation is attempted using biomaRt before falling back to
#'   OrgDb.
#'
#' @return A character vector containing mapped feature names.
#'   Gene symbols are returned when available; otherwise the cleaned
#'   Ensembl identifier is retained.
#'
#' @details
#' Annotation proceeds in two steps:
#'
#' \enumerate{
#'   \item If \code{mart} is supplied, gene symbols are retrieved using
#'   \code{biomaRt::getBM()}.
#'   \item If biomaRt is unavailable or fails, annotation is performed
#'   using \code{AnnotationDbi::select()} with the supplied
#'   OrgDb database.
#' }
#'
#' Duplicate mappings are resolved by keeping the first match.
#' Missing symbols are replaced by their corresponding Ensembl IDs.
#'
#' @seealso
#' \code{\link{map_entrez_ids}},
#' \code{\link{AnnotationDbi::select}}
#'
#' @keywords internal
#' @noRd
map_ensembl_to_feature <- function(
    ensembl_vec,
    orgdb,
    mart = NULL
) {

  log <- .msg(verbose = TRUE)

  ens_clean <- sub("\\..*$", "", ensembl_vec)
  mapping <- NULL

  if (!is.null(mart)) {

    log("Mapping Ensembl -> SYMBOL via biomaRt...")

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
      log("biomaRt mapping failed. Falling back to OrgDb.")
      NULL
    })
  }

  if (is.null(mapping)) {

    log("Mapping Ensembl -> SYMBOL using OrgDb...")

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


#' Map Ensembl identifiers to Entrez identifiers
#'
#' Annotates a differential expression result table by mapping
#' Ensembl gene identifiers to Entrez Gene identifiers using
#' a supplied Bioconductor OrgDb annotation database.
#'
#' Version suffixes in Ensembl identifiers (e.g. ENSG000001.5)
#' are automatically removed before annotation.
#'
#' @param res_df Differential expression result table containing
#'   an \code{ENSEMBL} column.
#' @param orgdb OrgDb annotation database (e.g.
#'   \code{org.Hs.eg.db}, \code{org.Mm.eg.db}).
#'
#' @return A data frame identical to \code{res_df} with an
#'   additional \code{ENTREZID} column.
#'
#' @seealso
#' \code{\link{map_ensembl_to_feature}},
#' \code{\link{run_go_enrichment}}
#'
#' @details
#' Gene annotation is performed using
#' \code{AnnotationDbi::select()} with
#' \code{keytype = "ENSEMBL"} and
#' \code{columns = "ENTREZID"}.
#'
#' When multiple mappings exist, the first match is retained.
#'
#' @keywords internal
#' @noRd
map_entrez_ids <- function(
    res_df,
    orgdb
) {

  log <- .msg(verbose = TRUE)

  log("Mapping ENTREZ IDs...")

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


#' @keywords internal
#' @noRd
.annotate_dds <- function(dds, orgdb) {

  genes <- rownames(dds)

  SummarizedExperiment::rowData(dds)$SYMBOL <-
    AnnotationDbi::mapIds(
      orgdb,
      keys = genes,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )

  SummarizedExperiment::rowData(dds)$ENTREZID <-
    AnnotationDbi::mapIds(
      orgdb,
      keys = genes,
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )

  SummarizedExperiment::rowData(dds)$GENENAME <-
    AnnotationDbi::mapIds(
      orgdb,
      keys = genes,
      column = "GENENAME",
      keytype = "ENSEMBL",
      multiVals = "first"
    )

  dds
}


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
