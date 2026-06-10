
generate_report <- function(
    sig_deg,
    go_results,
    output_dir = "DEGgo_out",
    report_template = NULL
) {
  
  log_msg("Generating HTML report...")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # -------------------------------------------------------
  # Find report template
  # -------------------------------------------------------
  
  if (is.null(report_template)) {
    
    candidates <- c(
      file.path(getwd(), "report_template.Rmd"),
      file.path(".", "report_template.Rmd"),
      file.path(output_dir, "report_template.Rmd"),
      system.file("report_template.Rmd", package = "cellDEverse"),
      system.file("inst", "report_template.Rmd", package = "cellDEverse")
    )
    
    report_template <- candidates[file.exists(candidates)][1]
  }
  
  if (is.na(report_template) || !file.exists(report_template)) {
    log_msg(
      "report_template.Rmd not found. Skipping HTML report.",
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
  go_path      <- file.path(output_dir, "GO_Barplot.png")
  
  has_volcano <- file.exists(volcano_path)
  has_pca     <- file.exists(pca_path)
  has_heatmap <- file.exists(heatmap_path)
  has_go      <- file.exists(go_path)
  
  if (!has_volcano) log_msg("Volcano plot missing. Report will continue.", type = "warn")
  if (!has_pca)     log_msg("PCA plot missing. Report will continue.", type = "warn")
  if (!has_heatmap) log_msg("Heatmap missing. Report will continue.", type = "warn")
  if (!has_go)      log_msg("GO barplot missing. Report will continue.", type = "warn")
  
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

