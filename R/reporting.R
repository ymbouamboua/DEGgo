
generate_report <- function(
    sig_deg,
    go_results,
    output_dir = "DEGgo_out",
    report_template = NULL
) {

  log <- .msg(verbose = TRUE)

  log("Generating HTML report...")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # -------------------------------------------------------
  # Find report template
  # -------------------------------------------------------

  if (is.null(report_template)) {

    candidates <- c(
      file.path(getwd(), "report.Rmd"),
      file.path(".", "report.Rmd"),
      file.path(output_dir, "report.Rmd"),
      system.file("report.Rmd", package = "DEGgo"),
      system.file("inst", "report.Rmd", package = "DEGgo")
    )

    report_template <- candidates[file.exists(candidates)][1]
  }

  if (is.na(report_template) || !file.exists(report_template)) {
    log(
      "report.Rmd not found. Skipping HTML report.",
      type = "warn"
    )
    return(invisible(NULL))
  }

  report_template <- normalizePath(
    report_template,
    winslash = "/",
    mustWork = TRUE
  )

  # -------------------------------------------------------
  # Plot paths
  # -------------------------------------------------------

  volcano_path <- file.path(output_dir, "Volcano_Plot.png")
  pca_path     <- file.path(output_dir, "PCA_Plot.png")
  heatmap_path <- file.path(output_dir, "Heatmap.png")
  go_path      <- file.path(output_dir, "GO_Dotplot.png")

  has_volcano <- file.exists(volcano_path)
  has_pca     <- file.exists(pca_path)
  has_heatmap <- file.exists(heatmap_path)
  has_go      <- file.exists(go_path)

  if (!has_volcano) log("Volcano plot missing. Report will continue.", type = "warn")
  if (!has_pca)     log("PCA plot missing. Report will continue.", type = "warn")
  if (!has_heatmap) log("Heatmap missing. Report will continue.", type = "warn")
  if (!has_go)      log("GO dotplot missing. Report will continue.", type = "warn")

  # -------------------------------------------------------
  # GO table
  # -------------------------------------------------------

  go_df <- tryCatch({
    if (is.null(go_results)) {
      data.frame()
    } else if (!is.null(go_results$go_results)) {
      as.data.frame(go_results$go_results)
    } else {
      as.data.frame(go_results)
    }
  }, error = function(e) {
    data.frame()
  })

  # -------------------------------------------------------
  # Render report
  # -------------------------------------------------------

  rmarkdown::render(
    input = report_template,
    output_file = "DEGgo_Report.html",
    output_dir = output_dir,
    params = list(
      sig_deg = sig_deg,
      go_results = go_df,

      volcano_path = if (has_volcano) volcano_path else NA_character_,
      pca_path = if (has_pca) pca_path else NA_character_,
      heatmap_path = if (has_heatmap) heatmap_path else NA_character_,
      go_barplot_path = if (has_go) go_path else NA_character_,

      has_volcano = has_volcano,
      has_pca = has_pca,
      has_heatmap = has_heatmap,
      has_go = has_go
    ),
    envir = new.env(parent = globalenv())
  )
}

