
interpretation_block <- function(section) {
  blocks <- list(
    qc = list(
      title = "QC interpretation",
      text = paste(
        "Samples should display comparable library sizes and numbers of detected genes.",
        "No strong outlier should dominate sample correlations or clustering.",
        "If obvious outliers are detected, they should be investigated before",
        "interpreting differential expression results."
      )
    ),
    pca = list(
      title = "PCA interpretation",
      text = paste(
        "Principal components should ideally separate samples according to the",
        "main biological variables of interest. Tight clustering of biological",
        "replicates supports the consistency of the dataset."
      )
    ),
    de = list(
      title = "Differential expression interpretation",
      text = paste(
        "Significant genes represent transcripts whose expression differs",
        "between the tested conditions according to the selected thresholds.",
        "Genes with large absolute log2 fold changes and low adjusted p-values",
        "are strong candidates for downstream validation."
      )
    ),
    volcano = list(
      title = "Volcano plot interpretation",
      text = paste(
        "Genes located in the upper-right region are significantly upregulated,",
        "whereas genes located in the upper-left region are significantly downregulated.",
        "Highlighted genes correspond to markers or genes of biological interest."
      )
    ),
    heatmap = list(
      title = "Heatmap interpretation",
      text = paste(
        "Samples sharing similar expression profiles should cluster together.",
        "Distinct expression patterns between experimental groups support",
        "the biological relevance of the detected differential expression."
      )
    ),
    go = list(
      title = "Gene Ontology interpretation",
      text = paste(
        "Enriched GO terms summarize biological processes overrepresented",
        "among differentially expressed genes. Related GO terms should be",
        "interpreted as broader biological themes rather than independent findings."
      )
    ),
    reproducibility = list(
      title = "Reproducibility",
      text = paste(
        "All analysis parameters, software versions and session information",
        "have been recorded to facilitate reproducibility and future re-analysis."
      )
    )
  )
  
  x <- blocks[[section]]
  
  if (!is.null(x)) {
    success_box(x$title, x$text)
  }
}


expected_box <- function(items, title = "Expected pattern") {
  html_box(
    title,
    paste0(
      "<ul>",
      paste0("<li>", items, "</li>", collapse = ""),
      "</ul>"
    ),
    type = "expected"
  )
}


section_note <- function(title, text = NULL) {
  html_box(title, text, type = "section")
}

info_box <- function(title, text) {
  html_box(title, text, type = "info")
}

success_box <- function(title, text) {
  html_box(title, text, type = "success")
}

warning_box <- function(title, text) {
  html_box(title, text, type = "warning")
}

