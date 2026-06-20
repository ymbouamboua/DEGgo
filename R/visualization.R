# ========================================================= #
# PLOT THEME
# ========================================================= #
#' Flexible ggplot2 theme used across DEGgo visualizations
#'
#' @param style Theme style.
#' @param txtsize Base text size.
#' @param xy.val Show axis tick labels.
#' @param x.ang X-axis text angle.
#' @param hjust,vjust Horizontal and vertical justification.
#' @param xlab,ylab Show x/y axis tick labels.
#' @param xy.lab Show all axis tick labels.
#' @param facet.face Facet label font face.
#' @param ttl.face Plot title font face.
#' @param txt.face Text font face.
#' @param ttl.pos Plot title position.
#' @param x.ttl,y.ttl Show x/y axis titles.
#' @param ticks,line,border Logical overrides.
#' @param grid.major,grid.minor Logical overrides for grid lines.
#' @param panel.fill Panel background fill.
#' @param facet.bg Show facet background.
#' @param mode Theme mode.
#' @param leg.pos Legend position.
#' @param leg.dir Legend direction.
#' @param leg.size Legend text size.
#' @param leg.ttl Legend title size.
#' @param leg.ttl.size Legend title text size.
#' @param leg.just Legend justification.
#' @param leg.ttl.text Optional legend title text.
#' @param ... Additional arguments passed to \code{ggplot2::theme()}.
#'
#' @return A ggplot2 theme object.
#' @keywords internal
#' @noRd
.deggo_theme <- function(
    style = c("classic", "minimal", "bw", "test", "void", "dirty", "gray"),
    txtsize = 12,
    xy.val = TRUE,
    x.ang = 0,
    hjust = NULL,
    vjust = NULL,
    xlab = TRUE,
    ylab = TRUE,
    xy.lab = TRUE,
    facet.face = "bold",
    ttl.face = "bold",
    txt.face = c("plain", "italic", "bold"),
    ttl.pos = c("center", "left", "right"),
    x.ttl = TRUE,
    y.ttl = TRUE,
    ticks = NULL,
    line = NULL,
    border = NULL,
    grid.major = NULL,
    grid.minor = NULL,
    panel.fill = "white",
    facet.bg = TRUE,
    mode = c("light", "dark"),
    leg.pos = "right",
    leg.dir = "vertical",
    leg.size = 10,
    leg.ttl = 10,
    leg.ttl.size = 10,
    leg.just = "center",
    leg.ttl.text = NULL,
    ...
) {

  style <- match.arg(style)
  ttl.pos <- match.arg(ttl.pos)
  txt.face <- match.arg(txt.face)
  mode <- match.arg(mode)

  lw <- 0.3

  if (is.null(line)) {
    line <- identical(style, "classic")
  }

  if (mode == "light") {
    col.txt <- "#1A1A1A"
    col.grid <- "#D9D9D9"
    col.panel <- panel.fill
    col.strip <- "#EFEFEF"
  } else {
    col.txt <- "#DDDDDD"
    col.grid <- "#444444"
    col.panel <- "#1E1E1E"
    col.strip <- "#383838"
  }

  if (is.null(hjust) || is.null(vjust)) {
    if (x.ang == 0) {
      hjust <- 0.5
      vjust <- 0.5
    } else if (x.ang == 45) {
      hjust <- 1
      vjust <- 1
    } else if (x.ang == 90) {
      hjust <- 1
      vjust <- 0.5
    } else if (x.ang == 270) {
      hjust <- 0
      vjust <- 0.5
    } else {
      hjust <- 1
      vjust <- 1
    }
  }

  ttl.hjust <- switch(
    ttl.pos,
    left = 0,
    center = 0.5,
    right = 1
  )

  base <- ggplot2::theme(
    text = ggplot2::element_text(
      color = col.txt,
      size = txtsize,
      family = "Helvetica",
      face = txt.face
    ),
    axis.text.x = ggplot2::element_text(color = col.txt, size = txtsize),
    axis.text.y = ggplot2::element_text(color = col.txt, size = txtsize),
    axis.title = ggplot2::element_text(size = txtsize),
    plot.title = ggplot2::element_text(
      hjust = ttl.hjust,
      face = ttl.face,
      size = txtsize + 2,
      color = col.txt
    ),
    strip.text = ggplot2::element_text(face = facet.face, color = col.txt),
    strip.background = ggplot2::element_rect(fill = col.strip, color = NA),
    panel.background = ggplot2::element_rect(fill = col.panel, color = NA),
    legend.title = ggplot2::element_text(size = leg.ttl.size, face = "bold"),
    legend.text = ggplot2::element_text(size = leg.size),
    legend.position = leg.pos,
    legend.direction = leg.dir,
    legend.justification = leg.just,
    legend.key.height = grid::unit(0.4, "cm"),
    legend.key.width = grid::unit(0.4, "cm"),
    legend.background = ggplot2::element_blank(),
    legend.box.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.spacing.y = grid::unit(0.05, "cm"),
    legend.margin = ggplot2::margin(1, 1, 1, 1),
    ...
  )

  preset <- switch(
    style,
    minimal = ggplot2::theme_minimal(base_size = txtsize),
    classic = ggplot2::theme_classic(base_size = txtsize),
    bw = ggplot2::theme_bw(base_size = txtsize),
    test = ggplot2::theme_test(base_size = txtsize),
    void = ggplot2::theme_void(base_size = txtsize),
    dirty = ggplot2::theme_minimal(base_size = txtsize) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      ),
    gray = ggplot2::theme_gray(base_size = txtsize) +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "#EDEDED", color = NA),
        panel.grid.major = ggplot2::element_line(color = "#CCCCCC", linewidth = lw),
        panel.grid.minor = ggplot2::element_line(color = "#DDDDDD", linewidth = lw / 2)
      )
  )

  th <- preset + base

  if (style %in% c("bw", "test", "gray")) {
    th <- th + ggplot2::theme(
      panel.border = ggplot2::element_rect(
        linewidth = lw,
        color = col.txt,
        fill = NA
      ),
      axis.line = ggplot2::element_blank()
    )
  }

  if (line && identical(style, "classic")) {
    th <- th + ggplot2::theme(
      axis.line.x = ggplot2::element_line(color = col.txt, linewidth = lw),
      axis.line.y = ggplot2::element_line(color = col.txt, linewidth = lw)
    )
  } else {
    th <- th + ggplot2::theme(axis.line = ggplot2::element_blank())
  }

  if (xy.val && xlab) {
    th <- th + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = x.ang,
        hjust = hjust,
        vjust = vjust
      )
    )
  }

  if (!xy.lab) {
    th <- th + ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  }

  if (!xlab) {
    th <- th + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }

  if (!ylab) {
    th <- th + ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  }

  if (!x.ttl) th <- th + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  if (!y.ttl) th <- th + ggplot2::theme(axis.title.y = ggplot2::element_blank())

  if (!is.null(ticks)) {
    th <- th + if (ticks) {
      ggplot2::theme(axis.ticks = ggplot2::element_line(color = col.txt, linewidth = lw))
    } else {
      ggplot2::theme(axis.ticks = ggplot2::element_blank())
    }
  }

  if (!is.null(border)) {
    th <- th + if (border) {
      ggplot2::theme(panel.border = ggplot2::element_rect(color = col.grid, fill = NA, linewidth = lw))
    } else {
      ggplot2::theme(panel.border = ggplot2::element_blank())
    }
  }

  if (!is.null(grid.major)) {
    th <- th + if (grid.major) {
      ggplot2::theme(panel.grid.major = ggplot2::element_line(color = col.grid, linewidth = lw))
    } else {
      ggplot2::theme(panel.grid.major = ggplot2::element_blank())
    }
  }

  if (!is.null(grid.minor)) {
    th <- th + if (grid.minor) {
      ggplot2::theme(panel.grid.minor = ggplot2::element_line(color = col.grid, linewidth = lw / 2))
    } else {
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    }
  }

  if (!facet.bg) {
    th <- th + ggplot2::theme(strip.background = ggplot2::element_blank())
  }

  if (!is.null(leg.ttl.text)) {
    th <- th + ggplot2::labs(color = leg.ttl.text, fill = leg.ttl.text)
  }

  if (identical(style, "void")) {
    th <- th + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )
  }

  th
}



# ========================================================= #
# DEGgo colors
# ========================================================= #

#' DEGgo default color palettes
#'
#' @return A named list of color vectors.
#' @keywords internal
#' @noRd
.deggo_colors <- function() {
  list(
    discrete = c(
      "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
      "#66A61E", "#E6AB02", "#A6761D", "#666666",
      "#1F78B4", "#B2DF8A", "#FB9A99", "#CAB2D6"
    ),
    regulation = c(
      Up = "#740001",
      Down = "#6497B1",
      NS = "grey70"
    ),
    qc = c(
      PASS = "#00A087",
      FAIL = "#E64B35",
      WARN = "#E6AB02"
    ),
    sex = c(
      Female = "#CC79A7",
      Male = "#0072B2",
      F = "#CC79A7",
      M = "#0072B2"
    ),
    treatment = c(
      Control = "#4DBBD5",
      Treated = "#D55E00",
      CTRL = "#4DBBD5",
      PBS = "#4DBBD5"
    )
  )
}



#' Make annotation colors
#'
#' @param annotation_col Annotation data frame.
#'
#' @return A named list of colors compatible with \code{pheatmap}.
#' @keywords internal
#' @noRd
.deggo_annotation_colors <- function(annotation_col) {
  if (is.null(annotation_col) || !ncol(annotation_col)) return(NULL)

  pal <- .deggo_colors()
  base <- pal[setdiff(names(pal), "discrete")]
  discrete <- pal$discrete

  out <- list()

  for (nm in colnames(annotation_col)) {
    levs <- levels(factor(annotation_col[[nm]]))

    known_cols <- NULL

    if (nm %in% names(base)) {
      known_cols <- base[[nm]]
    }

    if (!is.null(known_cols)) {
      matched <- known_cols[intersect(levs, names(known_cols))]
      missing <- setdiff(levs, names(matched))

      if (length(missing) > 0) {
        extra <- stats::setNames(
          rep(discrete, length.out = length(missing)),
          missing
        )
        matched <- c(matched, extra)
      }

      out[[nm]] <- matched[levs]

    } else {

      out[[nm]] <- stats::setNames(
        rep(discrete, length.out = length(levs)),
        levs
      )
    }
  }

  out
}


# ========================================================= #
# PLOT GENE EXPRESSION
# ========================================================= #
#' Plot normalized gene expression
#'
#' Generates publication-ready violin, boxplot, or barplot visualizations
#' from normalized expression values produced by
#' \code{extract_normalized_expression()}.
#'
#' Statistical comparisons can be added using
#' \code{ggpubr::stat_compare_means()}.
#'
#' @param expr_df Long-format expression table.
#' @param gene Gene symbol to plot.
#' @param x Metadata variable displayed on the x-axis.
#' @param color Metadata variable used for fill colors.
#' @param facet Optional metadata variable used for faceting.
#' @param geom Plot type: \code{"violin"}, \code{"boxplot"},
#'   or \code{"barplot"}.
#' @param comparisons List of comparisons passed to
#'   \code{ggpubr::stat_compare_means()}.
#' @param stats Logical. Add statistical comparisons.
#' @param stat_method Statistical test method.
#' @param stat_label Label style:
#'   \code{"p.signif"} or \code{"p.format"}.
#' @param output_dir Output directory.
#' @param filename Output filename.
#' @param width Plot width.
#' @param height Plot height.
#' @param dpi Plot resolution.
#' @param txtsize Base font size.
#' @param ncol Number of plot columns.
#' @param x_ang Numeric. Angle of x-axis labels.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_gene_expression(
#'   expr_df,
#'   gene = "Adipoq",
#'   x = "treatment",
#'   color = "treatment",
#'   facet = "tissue",
#'   geom = "violin"
#' )
#' }
#'
#' @export
#'
plot_gene_expression <- function(
    expr_df,
    gene,
    x = "treatment",
    color = x,
    facet = NULL,
    geom = c("violin", "boxplot", "barplot"),
    comparisons = list(c("PBS", "PAMH")),
    stats = TRUE,
    stat_method = "wilcox.test",
    stat_label = c("p.signif", "p.format"),
    output_dir = NULL,
    filename = NULL,
    width = 6,
    height = 5,
    dpi = 300,
    txtsize = 12,
    x_ang = 45,
    ncol = NULL
) {

  geom <- match.arg(geom)
  stat_label <- match.arg(stat_label)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }

  if (isTRUE(stats) && !requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is required for stats.", call. = FALSE)
  }

  df <- as.data.frame(expr_df)

  required <- c("gene", "expression", x)
  if (!is.null(color)) required <- c(required, color)
  if (!is.null(facet)) required <- c(required, facet)

  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  df <- df[df$gene == gene, , drop = FALSE]

  if (nrow(df) == 0) {
    stop("Gene not found in expr_df: ", gene, call. = FALSE)
  }

  df[[x]] <- factor(df[[x]])

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[x]], y = .data[["expression"]], fill = .data[[color]])
  )

  if (geom == "violin") {
    p <- p +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.8, linewidth = 0.1) +
      ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85) +
      ggplot2::geom_jitter(width = 0.08, size = 0.5, alpha = 0.8)
  }

  if (geom == "boxplot") {
    p <- p +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8, linewidth = 0.1) +
      ggplot2::geom_jitter(width = 0.12, size = 0.5, alpha = 0.8)
  }

  if (geom == "barplot") {
    p <- p +
      ggplot2::stat_summary(fun = mean, geom = "bar", alpha = 0.8, width = 0.7) +
      ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 0.2, linewidth = 0.1) +
      ggplot2::geom_jitter(width = 0.12, size = 0.5, alpha = 0.8)
  }

  if (!is.null(facet)) {
    p <- p +
      ggplot2::facet_wrap(stats::as.formula(paste("~", facet)), scales = "free_y", ncol = ncol)
  }

  if (isTRUE(stats)) {
    p <- p +
      ggpubr::stat_compare_means(
        comparisons = comparisons,
        method = stat_method,
        label = stat_label
      )
  }

  y_label <- switch(
    unique(df$assay)[1],
    raw = "Raw counts",
    normalized = "Normalized counts",
    log2_normalized = "log2(Normalized counts + 1)",
    vst = "VST expression",
    fpkm = "FPKM",
    tpm = "TPM",
    "Expression"
  )

  p <- p +
    ggplot2::labs(title = gene, x = x, y = y_label) +
    .deggo_theme(style = "classic", txtsize = txtsize, x.ang = x_ang)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    if (is.null(filename)) {
      filename <- paste0(gene, "_", geom)
    }

    ggplot2::ggsave(
      file.path(output_dir, paste0(filename, ".png")),
      p, width = width, height = height, dpi = dpi, bg = "white"
    )

    ggplot2::ggsave(
      file.path(output_dir, paste0(filename, ".pdf")),
      p, width = width, height = height, bg = "white"
    )
  }

  p
}




# ========================================================= #
# VOLCANO PLOT
# ========================================================= #
#' Generate volcano plot for differential expression results
#'
#' Creates a customizable volcano plot for DEGgo differential expression
#' results. Genes are automatically classified as Up, Down,
#' or Ns based on log2 fold-change and p-value thresholds.
#'
#' @param res_df Differential expression result data frame.
#' @param gene_col Column containing gene labels. If unavailable, the function
#'   automatically tries \code{"SYMBOL"}, \code{"gene"}, then \code{"ENSEMBL"}.
#' @param lfc_col Column containing log2 fold-change values.
#' @param pval_col Column containing raw p-values.
#' @param padj_col Column containing adjusted p-values.
#' @param use_padj Logical. If \code{TRUE}, use \code{padj_col}; otherwise use
#'   \code{pval_col}.
#' @param logfc_cutoff Absolute log2 fold-change threshold.
#' @param padj_cutoff P-value or adjusted p-value threshold.
#' @param top_n_labels Number of top significant genes to label.
#' @param genes_highlight Optional character vector of genes to highlight.
#' @param colors Named vector of colors for \code{Up},
#'   \code{Down}, and \code{Ns}.
#' @param point_size Point size.
#' @param alpha Point transparency.
#' @param title Plot title.
#' @param legend Logical. Show legend.
#' @param output_dir Output directory. If not \code{NULL}, the plot is saved.
#' @param filename Output file name without extension.
#' @param save_pdf Logical. Also save a PDF version.
#' @param width Plot width.
#' @param height Plot height.
#' @param txtsize Base text size.
#' @param style Plot style passed to \code{.deggo_theme()}.
#' @param dpi PNG resolution.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_volcano <- function(
    res_df,
    gene_col = NULL,
    lfc_col = "log2FoldChange",
    pval_col = "pvalue",
    padj_col = "padj",
    use_padj = TRUE,
    logfc_cutoff = 0.25,
    padj_cutoff = 0.05,
    top_n_labels = 40,
    genes_highlight = NULL,
    colors = c(
      "Up" = "#740001",
      "Down" = "#6497b1",
      "Ns" = "gray70"
    ),
    point_size = 1.5,
    alpha = 0.65,
    title = "Volcano Plot",
    legend = TRUE,
    output_dir = "DEGgo_out",
    filename = "Volcano_Plot",
    save_pdf = TRUE,
    width = 10,
    height = 8,
    dpi = 300,
    style = "classic",
    txtsize = 12
) {

  log <- .msg(verbose = TRUE, prefix = "DEGgo")
  log("Generating volcano plot...", type = "info")

  res_df <- as.data.frame(res_df, stringsAsFactors = FALSE)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required.", call. = FALSE)
  }

  res_df <- as.data.frame(res_df, stringsAsFactors = FALSE)

  if (is.null(gene_col)) {
    gene_col <- if ("SYMBOL" %in% colnames(res_df)) {
      "SYMBOL"
    } else if ("gene" %in% colnames(res_df)) {
      "gene"
    } else if ("ENSEMBL" %in% colnames(res_df)) {
      "ENSEMBL"
    } else {
      stop(
        "No gene label column found. Provide 'gene_col'.",
        call. = FALSE
      )
    }
  }

  p_col <- if (isTRUE(use_padj)) padj_col else pval_col

  required_cols <- c(gene_col, lfc_col, p_col)
  missing_cols <- setdiff(required_cols, colnames(res_df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  res_df[[lfc_col]] <- suppressWarnings(
    as.numeric(res_df[[lfc_col]])
  )

  res_df[[p_col]] <- suppressWarnings(
    as.numeric(res_df[[p_col]])
  )

  res_df <- res_df[
    !is.na(res_df[[lfc_col]]) &
      !is.na(res_df[[p_col]]) &
      res_df[[p_col]] > 0,
    ,
    drop = FALSE
  ]

  if (nrow(res_df) == 0) {
    stop("No valid rows available for volcano plot.", call. = FALSE)
  }

  res_df$volcano_status <- "Ns"

  res_df$volcano_status[
    res_df[[p_col]] < padj_cutoff &
      res_df[[lfc_col]] > logfc_cutoff
  ] <- "Up"

  res_df$volcano_status[
    res_df[[p_col]] < padj_cutoff &
      res_df[[lfc_col]] < -logfc_cutoff
  ] <- "Down"

  res_df$volcano_status <- factor(
    res_df$volcano_status,
    levels = c(
      "Up",
      "Down",
      "Ns"
    )
  )

  res_df$neg_log10_p <- -log10(res_df[[p_col]])

  counts <- table(res_df$volcano_status)

  legend_labels <- c(
    "Up" = paste0(
      "Up (",
      counts[["Up"]],
      ")"
    ),
    "Down" = paste0(
      "Down (",
      counts[["Down"]],
      ")"
    ),
    "Ns" = paste0(
      "Ns (",
      counts[["Ns"]],
      ")"
    )
  )

  max_abs <- max(
    abs(res_df[[lfc_col]]),
    logfc_cutoff,
    na.rm = TRUE
  ) + 0.5

  top_genes <- res_df[
    res_df$volcano_status != "Ns" &
      !is.na(res_df[[gene_col]]) &
      res_df[[gene_col]] != "",
    ,
    drop = FALSE
  ]

  if (nrow(top_genes) > 0 && top_n_labels > 0) {
    top_genes <- top_genes[
      order(top_genes[[p_col]]),
      ,
      drop = FALSE
    ]

    top_genes <- utils::head(
      top_genes,
      top_n_labels
    )
  } else {
    top_genes <- top_genes[0, , drop = FALSE]
  }

  p <- ggplot2::ggplot(
    res_df,
    ggplot2::aes(
      x = .data[[lfc_col]],
      y = .data[["neg_log10_p"]]
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        color = .data[["volcano_status"]]
      ),
      size = point_size,
      alpha = alpha
    ) +
    ggplot2::scale_color_manual(
      values = colors,
      labels = legend_labels,
      name = "Regulation",
      drop = FALSE
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      xintercept = c(
        -logfc_cutoff,
        logfc_cutoff
      ),
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.4
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-max_abs, max_abs)
    ) +
    ggplot2::labs(
      title = title,
      x = expression(log[2]~Fold~Change),
      y = if (isTRUE(use_padj)) {
        expression(-log[10]~adjusted~italic(P))
      } else {
        expression(-log[10]~italic(P))
      }
    ) + .deggo_theme(style = "classic",txtsize = txtsize)

  if (nrow(top_genes) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = top_genes,
        ggplot2::aes(
          label = .data[[gene_col]]
        ),
        size = 3,
        box.padding = 0.4,
        point.padding = 0.3,
        max.overlaps = 30,
        segment.color = "gray50"
      )
  }

  if (!is.null(genes_highlight)) {
    highlight_df <- res_df[
      res_df[[gene_col]] %in% genes_highlight,
      ,
      drop = FALSE
    ]

    if (nrow(highlight_df) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = highlight_df,
          shape = 21,
          size = point_size + 1.2,
          color = "black",
          stroke = 0.8
        ) +
        ggrepel::geom_text_repel(
          data = highlight_df,
          ggplot2::aes(
            label = .data[[gene_col]]
          ),
          fontface = "bold",
          size = 3.5,
          box.padding = 0.5,
          max.overlaps = Inf
        )
    }
  }

  if (!isTRUE(legend)) {
    p <- p + ggplot2::theme(
      legend.position = "none"
    )
  }

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(
        output_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )
    }

    ggplot2::ggsave(
      filename = file.path(
        output_dir,
        paste0(filename, ".png")
      ),
      plot = p,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )

    if (isTRUE(save_pdf)) {
      ggplot2::ggsave(
        filename = file.path(
          output_dir,
          paste0(filename, ".pdf")
        ),
        plot = p,
        width = width,
        height = height,
        bg = "white"
      )
    }
  }

  p
}


# ========================================================= #
# PCA PLOT
# ========================================================= #
#' Generate DEGgo PCA plot
#'
#' @param dds DESeq2 dataset.
#' @param metadata Sample metadata.
#' @param output_dir Output directory.
#' @param filename Output filename without extension.
#' @param intgroup Metadata columns used for PCA coloring/grouping.
#' @param title Plot title.
#' @param color_by Metadata column used for point color.
#' @param shape_by Metadata column used for point shape.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot resolution.
#'
#' @return A list containing the VST object, PCA data, and plot.
#' @export
plot_pca <- function(
    dds,
    metadata,
    output_dir = "DEGgo_out",
    filename = "PCA_Plot",
    intgroup = "condition",
    color_by = NULL,
    shape_by = NULL,
    title = "PCA Plot",
    width = 8,
    height = 6,
    dpi = 300
) {

  log <- .msg(verbose = TRUE, prefix = "DEGgo")
  log("Generating PCA plot...", type = "info")

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required.", call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!"sample" %in% colnames(metadata)) {
    metadata$sample <- rownames(metadata)
  }

  groups <- unique(c(intgroup, color_by, shape_by, "sample"))
  groups <- groups[!is.na(groups)]
  groups <- intersect(groups, colnames(metadata))

  if (!length(groups)) {
    groups <- "sample"
  }

  SummarizedExperiment::colData(dds)[, groups] <- metadata[colnames(dds), groups, drop = FALSE]

  vsd <- DESeq2::varianceStabilizingTransformation(
    dds,
    blind = FALSE
  )

  pca_data <- DESeq2::plotPCA(
    vsd,
    intgroup = groups,
    returnData = TRUE
  )

  percent_var <- round(100 * attr(pca_data, "percentVar"))

  if (is.null(color_by) || !color_by %in% colnames(pca_data)) {
    color_by <- groups[1]
  }

  if (!is.null(shape_by) && !shape_by %in% colnames(pca_data)) {
    shape_by <- NULL
  }

  pca_data[[color_by]] <- factor(pca_data[[color_by]])

  n_col <- nlevels(pca_data[[color_by]])

  pca_colors <- c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666",
    "#1F78B4", "#B2DF8A", "#FB9A99", "#CAB2D6"
  )

  if (n_col > length(pca_colors)) {
    pca_colors <- grDevices::hcl.colors(n_col, palette = "Dark 3")
  } else {
    pca_colors <- pca_colors[seq_len(n_col)]
  }

  names(pca_colors) <- levels(pca_data[[color_by]])

  pca_plot <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = .data[[color_by]]
    )
  ) +
    ggplot2::geom_point(size = 5, alpha = 0.9)

  if (!is.null(shape_by)) {
    pca_data[[shape_by]] <- factor(pca_data[[shape_by]])

    pca_plot <- ggplot2::ggplot(
      pca_data,
      ggplot2::aes(
        x = PC1,
        y = PC2,
        color = .data[[color_by]],
        shape = .data[[shape_by]]
      )
    ) +
      ggplot2::geom_point(size = 5, alpha = 0.9)
  }

  pca_plot <- pca_plot +
    ggplot2::scale_color_manual(values = pca_colors) +
    .deggo_theme(style = "classic",txtsize = 10) +
    ggplot2::labs(
      title = title,
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance"),
      color = color_by,
      shape = shape_by
    )

  if (identical(color_by, "sample")) {
    pca_plot <- pca_plot +
      ggplot2::geom_text(
        ggplot2::aes(label = sample),
        size = 2.8,
        vjust = -0.8,
        show.legend = FALSE
      )
  }

  ggplot2::ggsave(
    file.path(output_dir, paste0(filename, ".png")),
    pca_plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  ggplot2::ggsave(
    file.path(output_dir, paste0(filename, ".pdf")),
    pca_plot,
    width = width,
    height = height,
    bg = "white"
  )

  list(
    vsd = vsd,
    pca_data = pca_data,
    pca_plot = pca_plot
  )
}



# ========================================================= #
# .create_placeholder_png
# ========================================================= #
#' Create placeholder PNG
#'
#' @param file Output file.
#' @param text Text displayed.
#' @param width Image width.
#' @param height Image height.
#'
#' @keywords internal
#' @noRd
.create_placeholder_png <- function(
    file,
    text = "Plot not available",
    width = 1200,
    height = 900
  ) {
    grDevices::png(file, width = width, height = height, res = 150)
    graphics::plot.new()
    graphics::text(0.5, 0.5, labels = text, cex = 1.2)
    grDevices::dev.off()
    invisible(file)
  }

# ========================================================= #
# HEATMAP
# ========================================================= #

#' Generate DEGgo heatmap
#'
#' Generates a clustered heatmap of top differentially expressed genes.
#'
#' @param vsd Variance-stabilized dataset.
#' @param res_df Differential expression results data frame.
#' @param metadata Sample metadata data frame.
#' @param top_n_heatmap Number of genes displayed in the heatmap.
#' @param padj_cutoff Adjusted p-value threshold.
#' @param main Heatmap title.
#' @param output_dir Directory for exported plots.
#' @param filename Output file name without extension.
#' @param fallback Logical. If TRUE, use top ranked genes when no significant
#'   genes pass \code{padj_cutoff}.
#' @param contrast Contrast name to plot.
#' @param sample_subset Optional sample vector to retain.
#' @param metadata_filter Optional named list used to filter metadata.
#' @param annotation_cols Metadata columns shown as heatmap annotations.
#' @param annotation_colors Named list of annotation colors.
#' @param order_by Metadata columns used to order samples.
#' @param scale_rows Logical; scale rows before plotting.
#' @param cluster_rows Logical; cluster genes.
#' @param cluster_cols Logical; cluster samples.
#' @param fontsize_row Row label font size.
#' @param fontsize_col Column label font size.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#'
#' @return Heatmap expression matrix, invisibly.
#'
#' @export
plot_heatmap <- function(
    vsd,
    res_df,
    metadata,
    contrast = NULL,
    sample_subset = NULL,
    metadata_filter = NULL,
    top_n_heatmap = 50,
    padj_cutoff = 0.05,
    main = "Top Differentially Expressed Genes",
    output_dir = "DEGgo_out",
    filename = "Heatmap",
    fallback = TRUE,
    annotation_cols = c("condition", "treatment", "sex", "tissue"),
    annotation_colors = NULL,
    order_by = NULL,
    scale_rows = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 8,
    fontsize_col = 7,
    width = 8,
    height = 10
) {

  log <- .msg(verbose = TRUE, prefix = "DEGgo")
  log("Generating heatmap...", type = "info")

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required.", call. = FALSE)
  }

  mat <- SummarizedExperiment::assay(vsd)

  if (!is.null(sample_subset)) {

    sample_subset <- intersect(sample_subset, colnames(mat))
    mat <- mat[ ,sample_subset,drop = FALSE]
    metadata <- metadata[sample_subset,,drop = FALSE]
  }

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!is.null(metadata_filter)) {

    keep <- rep(TRUE, nrow(metadata))

    for (nm in names(metadata_filter)) {
      keep <- keep &
        metadata[[nm]] %in% metadata_filter[[nm]]
    }

    metadata <- metadata[keep, , drop = FALSE]
    mat <- mat[,rownames(metadata), drop = FALSE]
  }

  if (!"sample" %in% colnames(metadata)) {
    metadata$sample <- rownames(metadata)
  }

  rownames(metadata) <- metadata$sample

  if (!is.null(contrast)) {
    contrast_col <- contrast[1]
    contrast_levels <- contrast[2:3]

    if (!contrast_col %in% colnames(metadata)) {
      stop("Contrast column not found in metadata: ", contrast_col, call. = FALSE)
    }

    keep_samples <- rownames(metadata)[
      metadata[[contrast_col]] %in% contrast_levels
    ]

    keep_samples <- intersect(keep_samples, colnames(mat))

    if (!length(keep_samples)) {
      stop("No samples found for contrast: ", paste(contrast, collapse = " "), call. = FALSE)
    }

    mat <- mat[, keep_samples, drop = FALSE]
    metadata <- metadata[keep_samples, , drop = FALSE]
  } else {
    metadata <- metadata[colnames(mat), , drop = FALSE]
  }

  top <- res_df |>
    subset(!is.na(padj) & !is.na(log2FoldChange) & padj < padj_cutoff) |>
    {\(x) x[order(x$padj, -abs(x$log2FoldChange)), , drop = FALSE]}() |>
    utils::head(top_n_heatmap)

  if (!nrow(top) && isTRUE(fallback)) {
    top <- res_df |>
      subset(!is.na(padj) & !is.na(log2FoldChange)) |>
      {\(x) x[order(x$padj, -abs(x$log2FoldChange)), , drop = FALSE]}() |>
      utils::head(top_n_heatmap)
  }

  if (!nrow(top)) return(NULL)

  genes_use <- intersect(top$ENSEMBL, rownames(mat))
  if (!length(genes_use)) return(NULL)

  top <- top[match(genes_use, top$ENSEMBL), , drop = FALSE]
  mat_use <- mat[genes_use, , drop = FALSE]

  gene_labels <- ifelse(
    is.na(top$SYMBOL) | top$SYMBOL == "",
    top$ENSEMBL,
    top$SYMBOL
  )

  rownames(mat_use) <- make.unique(gene_labels)

  metadata_use <- metadata[colnames(mat_use), , drop = FALSE]

  if (!is.null(order_by)) {
    missing <- setdiff(order_by, colnames(metadata_use))
    if (length(missing)) {
      stop("order_by column(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    ord <- do.call(order, metadata_use[, order_by, drop = FALSE])
    metadata_use <- metadata_use[ord, , drop = FALSE]
    mat_use <- mat_use[, rownames(metadata_use), drop = FALSE]
  }

  annotation_col <- metadata_use[
    ,
    intersect(annotation_cols, colnames(metadata_use)),
    drop = FALSE
  ]

  if (ncol(annotation_col)) {
    annotation_col[] <- lapply(annotation_col, factor)
  } else {
    annotation_col <- NULL
  }

  if (isTRUE(scale_rows)) {
    mat_use <- t(scale(t(mat_use)))
    mat_use <- mat_use[stats::complete.cases(mat_use), , drop = FALSE]
  }

  if (nrow(mat_use) < 2) return(NULL)

  if (is.null(annotation_colors) || is.character(annotation_colors)) {
    annotation_colors <- .deggo_annotation_colors(annotation_col)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pheatmap::pheatmap(
    mat_use,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    treeheight_row = 0,
    treeheight_col = 0,
    angle_col = 90,
    main = main,
    border_color = "#9C9391",
    color = grDevices::colorRampPalette(
      c("#6497b1", "#F7F7F7", "#740001")
    )(100),
    breaks = seq(-2, 2, length.out = 101),
    #color = grDevices::colorRampPalette(c("purple", "black", "yellow"))(100),
    filename = file.path(output_dir, paste0(filename, ".png")),
    width = width,
    height = height
  )

  invisible(mat_use)
}


# ========================================================= #
# GENE EXPRESSION HEATMAP
# ========================================================= #
#' Plot expression heatmap for selected genes
#'
#' Generates a clustered heatmap for user-defined genes.
#'
#' @param counts Count matrix or count table.
#' @param metadata Sample metadata.
#' @param genes Genes to display.
#' @param gene_col Gene ID column.
#' @param feature_col Gene symbol column.
#' @param sample_col Sample column.
#' @param assay_transform Expression transformation.
#' @param annotation_cols Metadata columns shown above heatmap.
#' @param annotation_colors Named list of annotation colors.
#' @param order_by Metadata columns used to order samples.
#' @param output_dir Output directory.
#' @param filename Output filename.
#' @param main Heatmap title.
#' @param scale_rows Scale genes.
#' @param cluster_rows Cluster genes.
#' @param cluster_cols Cluster samples.
#' @param width Plot width.
#' @param height Plot height.
#' @param fontsize_row Row names size.
#' @param fontsize_col Column names size.
#' @param color Heatmap color palette.
#' @param breaks Numeric vector of color breaks.
#' @param border_color Border color for heatmap cells.
#'
#' @return Expression matrix used for plotting.
#' @export
#'
plot_gene_heatmap <- function(
    counts,
    metadata,
    genes,
    gene_col = c("gene_id", "GeneID", "gene", "Gene", "ENSEMBL", "ensembl", "ensembl_id"),
    feature_col = c("gene_name", "SYMBOL", "symbol", "gene_symbol", "external_gene_name"),
    sample_col = c("sample", "Sample", "SAMPLE"),
    assay_transform = c("log2", "log2cpm"),
    annotation_cols = c("condition", "treatment", "sex", "tissue"),
    annotation_colors = NULL,
    order_by = NULL,
    output_dir = "DEGgo_out",
    filename = "Gene_Expression_Heatmap",
    main = "Selected gene expression heatmap",
    color = grDevices::colorRampPalette(c("#6497b1", "#F7F7F7", "#740001"))(100),
    breaks = seq(-2, 2, length.out = 101),
    border_color = "#9C9391",
    scale_rows = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 8,
    fontsize_col = 7,
    width = 10,
    height = 7
) {

  assay_transform <- match.arg(assay_transform)

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required.", call. = FALSE)
  }

  counts <- as.data.frame(counts, check.names = FALSE)
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  gene_col <- gene_col[gene_col %in% colnames(counts)][1]
  feature_col <- feature_col[feature_col %in% colnames(counts)][1]
  sample_col <- sample_col[sample_col %in% colnames(metadata)][1]

  if (is.na(gene_col)) stop("No gene ID column found.", call. = FALSE)
  if (is.na(sample_col)) stop("No sample column found.", call. = FALSE)

  metadata[[sample_col]] <- as.character(metadata[[sample_col]])
  rownames(metadata) <- metadata[[sample_col]]

  sample_cols <- intersect(colnames(counts), metadata[[sample_col]])
  if (!length(sample_cols)) stop("No matching sample columns found.", call. = FALSE)

  mat <- as.matrix(counts[, sample_cols, drop = FALSE])
  rownames(mat) <- as.character(counts[[gene_col]])
  storage.mode(mat) <- "numeric"

  feature <- if (!is.na(feature_col)) as.character(counts[[feature_col]]) else rownames(mat)

  keep <- rownames(mat) %in% genes | toupper(feature) %in% toupper(genes)
  if (!any(keep)) stop("None of the requested genes were found.", call. = FALSE)

  mat_use <- mat[keep, , drop = FALSE]
  rownames(mat_use) <- make.unique(feature[keep])

  if (assay_transform == "log2") {
    mat_use <- log2(mat_use + 1)
  } else {
    mat_use <- log2(t(t(mat_use) / colSums(mat) * 1e6) + 1)
  }

  metadata_use <- metadata[colnames(mat_use), , drop = FALSE]

  if (!is.null(order_by)) {
    missing <- setdiff(order_by, colnames(metadata_use))
    if (length(missing)) {
      stop("order_by column(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    ord <- do.call(order, metadata_use[, order_by, drop = FALSE])
    metadata_use <- metadata_use[ord, , drop = FALSE]
    mat_use <- mat_use[, rownames(metadata_use), drop = FALSE]
  }

  annotation_col <- metadata_use[, intersect(annotation_cols, colnames(metadata_use)), drop = FALSE]
  if (!ncol(annotation_col)) annotation_col <- NULL

  if (isTRUE(scale_rows)) {
    mat_use <- t(scale(t(mat_use)))
    mat_use <- mat_use[stats::complete.cases(mat_use), , drop = FALSE]
  }

  if (!nrow(mat_use)) stop("No valid genes after scaling.", call. = FALSE)

  annotation_col <- metadata_use[, intersect(annotation_cols, colnames(metadata_use)), drop = FALSE]
  annotation_col[] <- lapply(annotation_col, factor)

  if (is.null(annotation_colors) || is.character(annotation_colors)) {
    annotation_colors <- .deggo_annotation_colors(annotation_col)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pheatmap::pheatmap(
    mat_use,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    treeheight_row = 0,
    treeheight_col = 0,
    angle_col = 90,
    main = main,
    color = color,
    breaks = breaks,
    border_color = border_color,
    filename = file.path(output_dir, paste0(filename, ".png")),
    width = width,
    height = height
  )
}



# ========================================================= #
# BULK RNA-SEQ QC
# ========================================================= #
#' Explore and QC bulk RNA-seq count data
#'
#' Generates pre-differential-expression quality control plots and metrics for
#' bulk RNA-seq count data, including library size, detected genes, expression
#' distribution, density curves, sample correlation, hierarchical clustering,
#' PCA plots, top variable gene heatmap, optional marker heatmap, and automatic
#' sample QC flags.
#'
#' @param counts Count table or matrix. Genes should be rows and samples columns.
#'   If \code{counts} is a data frame, one column must contain gene identifiers.
#' @param metadata Sample metadata data frame.
#' @param gene_col Candidate gene identifier columns in \code{counts}.
#' @param feature_col Candidate gene symbol/name columns in \code{counts}.
#' @param sample_col Candidate sample identifier columns in \code{metadata}.
#' @param annotation_cols Metadata columns shown as annotations in heatmaps.
#' @param color_by Metadata column used to color library size and density plots.
#'   If \code{NULL}, DEGgo tries tissue, condition, treatment, then sex.
#' @param markers Optional marker genes for marker heatmap.
#' @param top_variable_genes Number of top variable genes for heatmap.
#' @param output_dir Output directory.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot resolution.
#'
#' @examples
#' \dontrun{
#' qc <- explore_bulk_rnaseq(counts, metadata)
#'
#' cleaned <- remove_flagged_samples(
#'   counts = counts,
#'   metadata = metadata,
#'   qc_table = qc$qc
#' )
#'
#' counts <- cleaned$counts
#' metadata <- cleaned$metadata
#' }
#' @return Invisibly returns a list containing QC table, plots, correlation
#' matrix, PCA plots, top variable genes, marker heatmap object, and output path.
#'
#' @export
explore_bulk_rnaseq <- function(
    counts,
    metadata,
    gene_col = c("gene_id", "GeneID", "gene", "Gene", "ENSEMBL", "ensembl", "ensembl_id"),
    feature_col = c("gene_name", "SYMBOL", "symbol", "gene_symbol", "external_gene_name"),
    sample_col = c("sample", "Sample", "SAMPLE"),
    annotation_cols = c("condition", "treatment", "sex", "tissue"),
    color_by = NULL,
    markers = c("Adipoq", "Lep", "Ucp1", "Amh", "Amhr2", "Ar", "Akr1c3"),
    top_variable_genes = 500,
    output_dir = "DEGgo_QC",
    width = 8,
    height = 6,
    dpi = 300
) {

  for (pkg in c("ggplot2", "pheatmap")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required.", call. = FALSE)
    }
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  counts <- as.data.frame(counts, check.names = FALSE)
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  gene_col <- gene_col[gene_col %in% colnames(counts)][1]
  feature_col <- feature_col[feature_col %in% colnames(counts)][1]
  sample_col <- sample_col[sample_col %in% colnames(metadata)][1]

  if (is.na(gene_col)) stop("No gene ID column found.", call. = FALSE)
  if (is.na(sample_col)) stop("No sample column found.", call. = FALSE)

  rownames(metadata) <- as.character(metadata[[sample_col]])

  sample_cols <- intersect(colnames(counts), rownames(metadata))
  if (!length(sample_cols)) stop("No matching samples.", call. = FALSE)

  mat <- as.matrix(counts[, sample_cols, drop = FALSE])
  rownames(mat) <- as.character(counts[[gene_col]])
  storage.mode(mat) <- "numeric"

  metadata <- metadata[sample_cols, , drop = FALSE]
  logmat <- log2(mat + 1)

  if (is.null(color_by)) {
    color_by <- intersect(c("tissue", "condition", "treatment", "sex"), colnames(metadata))[1]
  }
  color_by <- if (!is.na(color_by) && color_by %in% colnames(metadata)) color_by else NULL

  qc <- data.frame(
    sample = colnames(mat),
    library_size = colSums(mat),
    detected_genes = colSums(mat > 0),
    percent_zero = 100 * colMeans(mat == 0),
    median_expression = apply(logmat, 2, stats::median),
    stringsAsFactors = FALSE
  )

  qc <- merge(qc, metadata, by.x = "sample", by.y = sample_col, all.x = TRUE)

  fill_aes <- if (!is.null(color_by)) ggplot2::aes(fill = .data[[color_by]]) else ggplot2::aes()
  color_aes <- if (!is.null(color_by)) ggplot2::aes(color = .data[[color_by]], group = sample) else ggplot2::aes(color = sample)

  p_lib <- ggplot2::ggplot(qc, ggplot2::aes(stats::reorder(sample, library_size), library_size / 1e6)) +
    fill_aes + ggplot2::geom_col() + ggplot2::coord_flip() +
    .deggo_theme(style = "classic", txtsize = 8) +
    ggplot2::labs(x = NULL, y = "Library size (million reads)", title = "Library size")

  p_det <- ggplot2::ggplot(qc, ggplot2::aes(stats::reorder(sample, detected_genes), detected_genes)) +
    fill_aes + ggplot2::geom_col() + ggplot2::coord_flip() +
    .deggo_theme(style = "classic", txtsize = 8) +
    ggplot2::labs(x = NULL, y = "Detected genes", title = "Detected genes")

  ggplot2::ggsave(file.path(output_dir, "Library_Size.png"), p_lib, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(file.path(output_dir, "Detected_Genes.png"), p_det, width = width, height = height, dpi = dpi)

  df_long <- data.frame(
    expression = as.vector(logmat),
    sample = rep(colnames(logmat), each = nrow(logmat)),
    stringsAsFactors = FALSE
  )

  df_long <- merge(df_long, metadata, by.x = "sample", by.y = sample_col, all.x = TRUE)

  p_box <- ggplot2::ggplot(df_long, ggplot2::aes(sample, expression)) +
    ggplot2::geom_boxplot(outlier.size = 0.15) +
    ggplot2::coord_flip() +
    .deggo_theme(style = "classic", txtsize = 8) +
    ggplot2::labs(x = NULL, y = "log2(count + 1)", title = "Expression distribution")

  p_density <- ggplot2::ggplot(df_long, ggplot2::aes(expression)) +
    color_aes + ggplot2::geom_density(linewidth = 0.4, alpha = 0.7) +
    .deggo_theme(style = "classic", txtsize = 8) +
    ggplot2::labs(x = "log2(count + 1)", y = "Density", title = "Expression density")

  ggplot2::ggsave(file.path(output_dir, "Expression_Boxplot.png"), p_box, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(file.path(output_dir, "Expression_Density.png"), p_density, width = width, height = height, dpi = dpi)

  annotation_col <- metadata[, intersect(annotation_cols, colnames(metadata)), drop = FALSE]
  if (ncol(annotation_col)) annotation_col[] <- lapply(annotation_col, factor) else annotation_col <- NULL
  ann_colors <- .deggo_annotation_colors(annotation_col)

  cor_mat <- stats::cor(logmat, method = "spearman")

  pheatmap::pheatmap(
    cor_mat,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    filename = file.path(output_dir, "Sample_Correlation_Heatmap.png"),
    main = "Sample Correlation",
    color = grDevices::colorRampPalette(c("#6497b1", "#F7F7F7", "#740001"))(100),
    fontsize = 8,
    width = 12,
    height = 10
  )

  grDevices::png(file.path(output_dir, "Hierarchical_Clustering.png"), width = 12, height = 6, units = "in", res = dpi)
  plot(stats::hclust(stats::dist(t(logmat))), main = "Hierarchical clustering", xlab = "", sub = "")
  grDevices::dev.off()

  vars <- apply(logmat, 1, stats::var, na.rm = TRUE)
  vars <- vars[is.finite(vars) & vars > 0]

  pca_plots <- NULL
  pca <- NULL
  pca_df <- NULL

  if (length(vars) >= 2) {

    pca_genes <- names(sort(vars, decreasing = TRUE))[seq_len(min(5000, length(vars)))]

    pca <- stats::prcomp(
      t(logmat[pca_genes, , drop = FALSE]),
      scale. = TRUE
    )

    meta_use <- metadata[rownames(pca$x), , drop = FALSE]
    meta_use <- meta_use[, !colnames(meta_use) %in% c("sample", "PC1", "PC2"), drop = FALSE]

    pca_df <- data.frame(
      sample = rownames(pca$x),
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      meta_use,
      check.names = FALSE
    )

    pca_df <- pca_df[, !duplicated(colnames(pca_df)), drop = FALSE]
    percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

    pca_plot <- function(color_by, shape_by = NULL, filename, title) {
      if (!color_by %in% colnames(pca_df)) return(NULL)

      pca_df[[color_by]] <- factor(pca_df[[color_by]])
      levs <- levels(pca_df[[color_by]])

      base_cols <- .deggo_colors()
      cols <- if (color_by %in% names(base_cols)) {
        tmp <- base_cols[[color_by]]
        missing <- setdiff(levs, names(tmp))
        if (length(missing)) tmp <- c(tmp, stats::setNames(grDevices::hcl.colors(length(missing), "Dark 3"), missing))
        tmp[levs]
      } else {
        stats::setNames(grDevices::hcl.colors(length(levs), "Dark 3"), levs)
      }

      aes_use <- ggplot2::aes(PC1, PC2, color = .data[[color_by]])

      if (!is.null(shape_by) && shape_by %in% colnames(pca_df)) {
        pca_df[[shape_by]] <- factor(pca_df[[shape_by]])
        aes_use <- ggplot2::aes(PC1, PC2, color = .data[[color_by]], shape = .data[[shape_by]])
      }

      p <- ggplot2::ggplot(pca_df, aes_use) +
        ggplot2::geom_point(size = 4, alpha = 0.9) +
        ggplot2::scale_color_manual(values = cols) +
        .deggo_theme(style = "classic", txtsize = 8) +
        ggplot2::labs(
          title = title,
          x = paste0("PC1: ", percent_var[1], "% variance"),
          y = paste0("PC2: ", percent_var[2], "% variance"),
          color = color_by,
          shape = shape_by
        )

      if (identical(color_by, "sample")) {
        p <- p + ggplot2::geom_text(
          ggplot2::aes(label = sample),
          size = 2.5,
          vjust = -0.8,
          show.legend = FALSE
        )
      }

      ggplot2::ggsave(file.path(output_dir, paste0(filename, ".png")), p, width = width, height = height, dpi = dpi)
      ggplot2::ggsave(file.path(output_dir, paste0(filename, ".pdf")), p, width = width, height = height)

      p
    }

    pca_plots <- list(
      sample = pca_plot("sample", NULL, "PCA_by_sample", "PCA by sample"),
      tissue = pca_plot("tissue", NULL, "PCA_by_tissue", "PCA by tissue"),
      treatment = pca_plot("treatment", NULL, "PCA_by_treatment", "PCA by treatment"),
      sex = pca_plot("sex", NULL, "PCA_by_sex", "PCA by sex"),
      tissue_treatment = pca_plot("tissue", "treatment", "PCA_tissue_treatment", "PCA tissue + treatment")
    )
  }

  top <- names(sort(vars, decreasing = TRUE))[seq_len(min(top_variable_genes, length(vars)))]

  if (length(top) >= 2) {
    top_mat <- logmat[top, , drop = FALSE]
    top_mat <- t(scale(t(top_mat)))
    top_mat <- top_mat[stats::complete.cases(top_mat), , drop = FALSE]

    if (nrow(top_mat) >= 2) {
      pheatmap::pheatmap(
        top_mat,
        annotation_col = annotation_col,
        annotation_colors = ann_colors,
        show_rownames = FALSE,
        show_colnames = FALSE,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        #color = grDevices::colorRampPalette(c("purple", "black", "yellow"))(100),
        color = grDevices::colorRampPalette(
          c("#6497b1", "#F7F7F7", "#740001")
        )(100),
        breaks = seq(-2, 2, length.out = 101),
        filename = file.path(output_dir, "TopVariableGenes_Heatmap.png"),
        width = 9,
        height = 10
      )
    }
  }

  marker_heatmap <- NULL

  if (!is.null(markers) && length(markers) && exists("plot_gene_heatmap")) {
    marker_heatmap <- tryCatch(
      plot_gene_heatmap(
        counts = counts,
        metadata = metadata,
        genes = markers,
        gene_col = gene_col,
        feature_col = feature_col,
        sample_col = sample_col,
        annotation_cols = annotation_cols,
        annotation_colors = ann_colors,
        order_by = intersect(c("tissue", "sex", "treatment", "condition"), colnames(metadata)),
        color = grDevices::colorRampPalette(
          c("#6497b1", "#F7F7F7", "#740001")
        )(100),
        breaks = seq(-2, 2, length.out = 101),
        output_dir = output_dir,
        filename = "Marker_Heatmap",
        main = "Marker gene expression"
      ),
      error = function(e) NULL
    )
  }

  qc$flag_library <- .flag_low_mad(qc$library_size)
  qc$flag_genes <- .flag_low_mad(qc$detected_genes)

  sample_cor <- apply(cor_mat, 1, function(x) mean(x[-which.max(x)], na.rm = TRUE))
  qc$mean_correlation <- sample_cor[match(qc$sample, names(sample_cor))]
  qc$flag_correlation <- .flag_low_mad(qc$mean_correlation)

  if (!is.null(pca)) {
    pc_df <- data.frame(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
    center <- c(mean(pc_df$PC1), mean(pc_df$PC2))
    pc_df$pca_distance <- sqrt((pc_df$PC1 - center[1])^2 + (pc_df$PC2 - center[2])^2)
    qc$pca_distance <- pc_df$pca_distance[match(qc$sample, pc_df$sample)]
    qc$flag_pca <- .flag_high_mad(qc$pca_distance)
  } else {
    qc$pca_distance <- NA_real_
    qc$flag_pca <- FALSE
  }

  flag_cols <- grep("^flag_", colnames(qc), value = TRUE)
  qc$n_flags <- rowSums(qc[, flag_cols, drop = FALSE], na.rm = TRUE)
  qc$qc_status <- ifelse(qc$n_flags > 0, "FAIL", "PASS")
  qc$recommend_remove <- qc$n_flags >= 2

  utils::write.table(qc, file.path(output_dir, "QC_Flags.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(qc[qc$recommend_remove, , drop = FALSE], file.path(output_dir, "Samples_Recommended_For_Removal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  p_qc <- ggplot2::ggplot(qc, ggplot2::aes(stats::reorder(sample, n_flags), n_flags, fill = qc_status)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = .deggo_colors()$qc) +
    .deggo_theme(style = "classic", txtsize = 8) +
    ggplot2::labs(x = NULL, y = "Number of QC flags", fill = "QC status", title = "Sample QC flags")

  ggplot2::ggsave(file.path(output_dir, "QC_Flags.png"), p_qc, width = width, height = height, dpi = dpi)

  utils::write.table(qc, file.path(output_dir, "QC_Summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(list(
    qc = qc,
    library_size = p_lib,
    detected_genes = p_det,
    boxplot = p_box,
    density = p_density,
    correlation = cor_mat,
    pca = pca_plots,
    top_variable_genes = top,
    marker_heatmap = marker_heatmap,
    qc_flags_plot = p_qc,
    output_dir = output_dir
  ))
}



# ========================================================= #
# MARKER SCORE CHECK
# ========================================================= #
#' Validate sample identity using marker gene signatures
#'
#' Computes signature scores from predefined marker gene sets,
#' predicts the most likely sample identity, and flags potential
#' sample swaps or annotation mismatches.
#'
#' For each marker set, the function calculates the mean expression
#' of all detected marker genes across samples. The predicted group
#' corresponds to the marker set with the highest score.
#'
#' @param counts A data frame or matrix containing gene expression counts.
#'   Rows correspond to genes and columns correspond to samples.
#' @param metadata A data frame containing sample metadata.
#' @param marker_sets Named list of marker gene vectors. Each list element
#'   defines a biological group or tissue signature.
#' @param sample_col Column name in \code{metadata} containing sample IDs.
#'   Default is \code{"sample"}.
#' @param group_col Column name in \code{metadata} containing the expected
#'   sample annotation (e.g. tissue, condition). Default is \code{"tissue"}.
#' @param feature_col Column name in \code{counts} containing gene symbols.
#'   Default is \code{"gene_name"}.
#' @param log_transform Logical; if \code{TRUE}, computes scores from
#'   \code{log2(count + 1)} transformed values. Default is \code{TRUE}.
#' @param plot Logical; if \code{TRUE}, returns a barplot of marker scores.
#'   Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{scores}{Data frame with marker scores, predicted groups,
#'   and swap flags.}
#'   \item{long}{Long-format version of the score table.}
#'   \item{swaps}{Subset of samples where predicted and expected
#'   groups differ.}
#'   \item{plot}{A ggplot object showing marker scores per sample
#'   if \code{plot = TRUE}; otherwise \code{NULL}.}
#' }
#'
#' @details
#' The function performs a simple signature-based classification:
#' \enumerate{
#'   \item Matches marker genes against the feature column.
#'   \item Computes the mean expression of matched genes for each sample.
#'   \item Predicts the group with the highest signature score.
#'   \item Flags samples where predicted and expected annotations differ.
#' }
#'
#' Gene matching is case-insensitive.
#'
#' @examples
#' \dontrun{
#' marker_sets <- list(
#'   BAT = c("Ucp1", "Cidea", "Ppargc1a"),
#'   WAT = c("Adipoq", "Lep", "Fabp4"),
#'   TESTIS = c("Amh", "Sox9", "Ddx4"),
#'   OVARY = c("Foxl2", "Fshr", "Bmp15")
#' )
#'
#' res <- marker_score_check(
#'   counts = counts,
#'   metadata = metadata,
#'   marker_sets = marker_sets,
#'   sample_col = "sample",
#'   group_col = "tissue"
#' )
#'
#' head(res$scores)
#' res$swaps
#' print(res$plot)
#'}
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs
#'
#' @export
marker_score_check <- function(
    counts,
    metadata,
    marker_sets,
    sample_col = "sample",
    group_col = "tissue",
    feature_col = "gene_name",
    log_transform = TRUE,
    plot = TRUE
) {

  counts <- as.data.frame(counts, check.names = FALSE)
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  sample_cols <- metadata[[sample_col]]
  sample_cols <- sample_cols[sample_cols %in% colnames(counts)]

  if (length(sample_cols) == 0) {
    stop("No metadata samples found in counts columns.", call. = FALSE)
  }

  gene_names <- toupper(counts[[feature_col]])

  compute_one_score <- function(genes) {

    keep <- gene_names %in% toupper(genes)

    if (!any(keep)) {
      return(rep(NA_real_, length(sample_cols)))
    }

    mat <- as.matrix(counts[keep, sample_cols, drop = FALSE])
    storage.mode(mat) <- "numeric"

    if (log_transform) {
      mat <- log2(mat + 1)
    }

    colMeans(mat, na.rm = TRUE)
  }

  scores <- data.frame(
    sample = sample_cols,
    stringsAsFactors = FALSE
  )

  for (nm in names(marker_sets)) {
    scores[[paste0(nm, "_score")]] <- compute_one_score(marker_sets[[nm]])
  }

  md <- metadata[, c(sample_col, group_col), drop = FALSE]
  colnames(md) <- c("sample", "true_group")

  scores <- merge(scores, md, by = "sample", all.x = TRUE)

  score_cols <- paste0(names(marker_sets), "_score")
  score_mat <- as.matrix(scores[, score_cols, drop = FALSE])

  scores$predicted_group <- names(marker_sets)[
    max.col(score_mat, ties.method = "first")
  ]

  scores$possible_swap <- (
    scores$true_group != scores$predicted_group
  )

  scores_long <- reshape2::melt(
    scores,
    id.vars = c(
      "sample",
      "true_group",
      "predicted_group",
      "possible_swap"
    ),
    measure.vars = score_cols,
    variable.name = "marker_set",
    value.name = "score"
  )

  scores_long$marker_set <- sub(
    "_score$",
    "",
    scores_long$marker_set
  )

  p <- NULL

  if (plot) {

    p <- ggplot2::ggplot(
      scores_long,
      ggplot2::aes(
        x = sample,
        y = score,
        fill = marker_set
      )
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = NULL,
        y = "Mean marker score",
        fill = "Marker set"
      ) +
      .deggo_theme()
  }

  list(
    scores = scores,
    long = scores_long,
    swaps = scores[scores$possible_swap, , drop = FALSE],
    plot = p
  )
}




# ========================================================= #
# PLOT GO TERMS
# ========================================================= #
#' Plot GO terms by regulation status
#'
#' @param go_df GO enrichment data frame from \code{run_go_enrichment()}.
#' @param go_terms_of_interest Optional GO terms to display.
#' @param comparison Optional comparison name used as plot title.
#' @param top_n Number of GO terms to show per regulation.
#' @param wrap_width Width for wrapping GO term labels.
#' @param color_values Named colors for Up and Down.
#' @param size_range Point size range.
#' @param style Theme style passed to \code{.deggo_theme()}.
#' @param x_angle X-axis text angle.
#' @param font_size Base font size.
#'
#' @return A ggplot object.
#' @export
plot_go_terms <- function(
    go_df,
    go_terms_of_interest = NULL,
    comparison = NULL,
    top_n = 10,
    wrap_width = 50,
    color_values = c(
      "Up" = "#740001",
      "Down" = "#6497b1"
    ),
    size_range = c(2, 6),
    style = "bw",
    x_angle = 0,
    font_size = 12
) {

  for (pkg in c("ggplot2", "dplyr", "stringr", "forcats")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required.", call. = FALSE)
    }
  }

  if (is.null(go_df) || !is.data.frame(go_df) || nrow(go_df) == 0) {
    stop("go_df is empty.", call. = FALSE)
  }

  required <- c("Description", "FoldEnrichment", "p.adjust", "Count", "Regulation")
  missing <- setdiff(required, colnames(go_df))

  if (length(missing) > 0) {
    stop(
      "Missing column(s) in go_df: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  df <- go_df |>
    dplyr::mutate(
      Description = stringr::str_to_upper(.data[["Description"]]),
      Regulation = factor(.data[["Regulation"]], levels = c("Up", "Down")),
      FoldEnrichment = as.numeric(.data[["FoldEnrichment"]]),
      p.adjust = as.numeric(.data[["p.adjust"]]),
      Count = as.integer(.data[["Count"]])
    ) |>
    dplyr::filter(
      !is.na(.data[["Regulation"]]),
      is.finite(.data[["FoldEnrichment"]]),
      is.finite(.data[["p.adjust"]]),
      .data[["p.adjust"]] > 0
    )

  if (nrow(df) == 0) {
    stop("No valid GO terms after cleaning.", call. = FALSE)
  }

  if (!is.null(go_terms_of_interest)) {
    go_terms_of_interest <- stringr::str_to_upper(go_terms_of_interest)

    df <- df |>
      dplyr::filter(.data[["Description"]] %in% go_terms_of_interest)
  }

  df <- df |>
    dplyr::group_by(.data[["Regulation"]], .data[["Description"]]) |>
    dplyr::arrange(
      .data[["p.adjust"]],
      dplyr::desc(.data[["FoldEnrichment"]]),
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()

  if (is.null(go_terms_of_interest)) {
    df <- df |>
      dplyr::group_by(.data[["Regulation"]]) |>
      dplyr::arrange(
        .data[["p.adjust"]],
        dplyr::desc(.data[["FoldEnrichment"]]),
        .by_group = TRUE
      ) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::ungroup()
  }

  df <- df |>
    dplyr::mutate(
      log10FDR = -log10(.data[["p.adjust"]]),
      Description_wrapped = stringr::str_wrap(
        .data[["Description"]],
        width = wrap_width
      ),
      Regulation = droplevels(.data[["Regulation"]]),
      Description_wrapped = forcats::fct_reorder(
        .data[["Description_wrapped"]],
        .data[["log10FDR"]]
      )
    )

  present_regs <- unique(as.character(df$Regulation))
  color_values <- color_values[names(color_values) %in% present_regs]

  if (nrow(df) == 0) {
    stop("No GO terms left after filtering.", call. = FALSE)
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["log10FDR"]],
      y = .data[["Description_wrapped"]]
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = .data[["log10FDR"]],
        yend = .data[["Description_wrapped"]],
        color = .data[["Regulation"]]
      ),
      linewidth = 0.8,
      alpha = 1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = .data[["FoldEnrichment"]],
        color = .data[["Regulation"]]
      ),
      alpha = 1
    ) +
    ggplot2::scale_color_manual(
      values = color_values,
      drop = TRUE
    ) +
    ggplot2::scale_size(range = size_range) +
    ggplot2::labs(
      x = expression(-log[10](FDR)),
      y = NULL,
      size = "Fold enrichment",
      color = "Regulation",
      title = comparison %||% "GO enrichment"
    ) +
    .deggo_theme(style = style, x.ang = x_angle, txtsize = font_size)
}


# ========================================================= #
# PLOT ALL GO TERMS
# ========================================================= #
#' Plot GO enrichment terms for all pairwise comparisons
#'
#' Generates GO enrichment plots for every comparison returned by
#' \code{run_deggo()} and optionally exports figures to disk.
#'
#' @param results DEGgo results object returned by \code{run_deggo()}.
#' @param top_n Number of GO terms to display per regulation class
#'   (Up and Down). Default is \code{10}.
#' @param font_size Base font size used in plots. Default is \code{8}.
#' @param output_dir Optional output directory for exported figures.
#'   If \code{NULL}, plots are returned but not saved.
#' @param width Plot width in inches. Default is \code{8}.
#' @param height Plot height in inches. Default is \code{6}.
#' @param dpi Plot resolution. Default is \code{300}.
#'
#' @return
#' A named list of ggplot objects, one per comparison.
#'
#' @details
#' This function iterates over all pairwise comparisons stored in
#' \code{results$go_results} and calls \code{plot_go_terms()} on each
#' GO enrichment table.
#'
#' If \code{output_dir} is supplied, PNG and PDF versions of all plots
#' are automatically exported.
#'
#' @examples
#' \dontrun{
#'
#' go_plots <- plot_all_go_terms(
#'   results = results,
#'   top_n = 10,
#'   font_size = 8,
#'   output_dir = "GO_plots"
#' )
#'
#' names(go_plots)
#'
#' go_plots$WAT_Female_PAMH_vs_PBS
#'
#' }
#'
#' @export
plot_all_go_terms <- function(
    results,
    top_n = 10,
    font_size = 8,
    output_dir = NULL,
    width = 8,
    height = 6,
    dpi = 300
) {

  plots <- lapply(names(results$go_results), function(comp) {

    go_df <- results$go_results[[comp]]$go_results

    if (is.null(go_df) || nrow(go_df) == 0) {
      return(NULL)
    }

    p <- plot_go_terms(
      go_df = go_df,
      comparison = comp,
      top_n = top_n,
      font_size = font_size
    )

    if (!is.null(output_dir)) {

      dir.create(
        output_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )

      ggplot2::ggsave(
        filename = file.path(
          output_dir,
          paste0(comp, "_GO_terms.png")
        ),
        plot = p,
        width = width,
        height = height,
        dpi = dpi
      )

      ggplot2::ggsave(
        filename = file.path(
          output_dir,
          paste0(comp, "_GO_terms.pdf")
        ),
        plot = p,
        width = width,
        height = height
      )
    }

    p
  })

  names(plots) <- names(results$go_results)

  plots
}

