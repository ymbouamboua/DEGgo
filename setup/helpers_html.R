

html_box <- function(title, text = NULL, type = "info") {
  if (is.null(text)) {
    text <- title
    title <- switch(
      type,
      section = "About this section",
      expected = "Expected pattern",
      success = "Interpretation",
      warning = "Warning",
      "Information"
    )
  }
  
  cat(
    '<div class="', type, '-box">',
    '<div class="', type, '-title">', html_escape(title), '</div>',
    '<div class="', type, '-text">', text, '</div>',
    '</div>'
  )
}


html_card <- function(title, value, subtitle = NULL) {
  cat(
    '<div class="metric-card">',
    '<div class="metric-title">', html_escape(title), '</div>',
    '<div class="metric-value">', html_escape(value), '</div>',
    if (!is.null(subtitle)) {
      paste0('<div class="metric-subtitle">', html_escape(subtitle), '</div>')
    } else {
      ''
    },
    '</div>'
  )
}



html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}


