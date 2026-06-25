
plot_gallery <- function(files, badge = NULL, empty = "No figures found.") {
  files <- files[file.exists(files)]
  
  if (!length(files)) {
    cat("<p><em>", empty, "</em></p>")
    return(invisible(NULL))
  }
  
  cat('<div class="figure-grid">')
  
  for (f in sort(files)) {
    title <- clean_title(f)
    img_src <- html_img_src(f)
    badge_i <- badge %||% figure_badge(f)
    
    cat('<div class="figure-card">')
    cat('<div class="figure-header">')
    cat('<span class="figure-badge">', html_escape(badge_i), '</span>')
    cat('<h3>', html_escape(title), '</h3>')
    cat('</div>')
    
    cat('<a href="', img_src, '" target="_blank">')
    cat('<img class="figure-img" src="', img_src, '" alt="', html_escape(title), '" title="Click to open full-size figure">')
    cat('</a>')
    
    cat('<div class="figure-footer">')
    cat('<a href="', img_src, '" target="_blank">Open full figure</a>')
    cat('</div>')
    
    cat('</div>')
  }
  
  cat('</div>')
  invisible(NULL)
}


clean_title <- function(x) {
  x <- tools::file_path_sans_ext(basename(x))
  x <- gsub("_", " ", x)
  x <- gsub("Volcano Plot", "", x, ignore.case = TRUE)
  x <- gsub("GO terms", "", x, ignore.case = TRUE)
  x <- gsub("Heatmap", "", x, ignore.case = TRUE)
  x <- gsub("Sample Correlation", "Sample correlation", x, ignore.case = TRUE)
  x <- gsub("Hierarchical Clustering", "Hierarchical clustering", x, ignore.case = TRUE)
  trimws(x)
}


figure_badge <- function(f) {
  b <- basename(f)
  
  if (grepl("Correlation|Clustering|QC|Library|Detected", b, ignore.case = TRUE)) return("QC")
  if (grepl("PCA", b, ignore.case = TRUE)) return("PCA")
  if (grepl("Volcano", b, ignore.case = TRUE)) return("Volcano")
  if (grepl("Heatmap", b, ignore.case = TRUE)) return("Heatmap")
  if (grepl("GO", b, ignore.case = TRUE)) return("GO")
  
  "Figure"
}


