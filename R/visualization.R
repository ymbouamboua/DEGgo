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
    line <- style == "classic"
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

  ttl.pos <- switch(
    ttl.pos,
    left = 0,
    center = 0.5,
    right = 1
  )

  base <- ggplot2::theme(
    text = ggplot2::element_text(
      color = col.txt,
      size = txtsize,
      family = "Helvetica"
    ),
    axis.text.x = ggplot2::element_text(
      color = col.txt,
      size = txtsize
    ),
    axis.text.y = ggplot2::element_text(
      color = col.txt,
      size = txtsize
    ),
    axis.title = ggplot2::element_text(
      size = txtsize
    ),
    plot.title = ggplot2::element_text(
      hjust = ttl.pos,
      face = ttl.face,
      size = txtsize + 2,
      color = col.txt
    ),
    strip.text = ggplot2::element_text(
      face = facet.face,
      color = col.txt
    ),
    legend.title = ggplot2::element_text(
      size = leg.ttl.size + 2,
      face = "bold"
    ),
    legend.text = ggplot2::element_text(
      size = leg.size
    ),
    legend.position = leg.pos,
    legend.direction = leg.dir,
    legend.justification = leg.just,
    legend.key.height = grid::unit(0.4, "cm"),
    legend.key.width = grid::unit(0.4, "cm"),
    legend.background = ggplot2::element_blank(),
    legend.box.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.box = "vertical",
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
        panel.background = ggplot2::element_rect(
          fill = "#EDEDED",
          color = NA
        ),
        panel.grid.major = ggplot2::element_line(
          color = "#CCCCCC",
          linewidth = lw
        ),
        panel.grid.minor = ggplot2::element_line(
          color = "#DDDDDD",
          linewidth = lw / 2
        )
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

  if (line && style == "classic") {
    th <- th + ggplot2::theme(
      axis.line.x = ggplot2::element_line(
        color = col.txt,
        linewidth = lw
      ),
      axis.line.y = ggplot2::element_line(
        color = col.txt,
        linewidth = lw
      )
    )
  } else {
    th <- th + ggplot2::theme(
      axis.line = ggplot2::element_blank()
    )
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

  if (!x.ttl) {
    th <- th + ggplot2::theme(
      axis.title.x = ggplot2::element_blank()
    )
  }

  if (!y.ttl) {
    th <- th + ggplot2::theme(
      axis.title.y = ggplot2::element_blank()
    )
  }

  if (!is.null(ticks)) {
    th <- th + if (isTRUE(ticks)) {
      ggplot2::theme(
        axis.ticks = ggplot2::element_line(
          color = col.txt,
          linewidth = lw
        )
      )
    } else {
      ggplot2::theme(axis.ticks = ggplot2::element_blank())
    }
  }

  if (!is.null(border)) {
    th <- th + if (isTRUE(border)) {
      ggplot2::theme(
        panel.border = ggplot2::element_rect(
          color = col.grid,
          fill = NA,
          linewidth = lw
        )
      )
    } else {
      ggplot2::theme(panel.border = ggplot2::element_blank())
    }
  }

  if (!is.null(grid.major)) {
    th <- th + if (isTRUE(grid.major)) {
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(
          color = col.grid,
          linewidth = lw
        )
      )
    } else {
      ggplot2::theme(panel.grid.major = ggplot2::element_blank())
    }
  }

  if (!is.null(grid.minor)) {
    th <- th + if (isTRUE(grid.minor)) {
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_line(
          color = col.grid,
          linewidth = lw / 2
        )
      )
    } else {
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    }
  }

  if (!facet.bg) {
    th <- th + ggplot2::theme(
      strip.background = ggplot2::element_blank()
    )
  }

  if (!is.null(leg.ttl.text)) {
    th <- th + ggplot2::labs(color = leg.ttl.text)
  }

  if (style == "void") {
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

  if (style == "dirty") {
    th <- th + ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )
  }

  th
}


# ========================================================= #
# PUBLICATION DEFAULTS
# ========================================================= #

.deggo_plot_defaults <- function(
    context = c("single", "multi", "heatmap", "report")
) {
  context <- match.arg(context)

  switch(
    context,
    single = list(width = 6.5, height = 5.0, dpi = 300, txtsize = 8),
    multi  = list(width = 8.0, height = 6.0, dpi = 300, txtsize = 8),
    heatmap = list(width = 7.0, height = 8.0, dpi = 300, txtsize = 7),
    report = list(width = 9.0, height = 6.5, dpi = 300, txtsize = 8)
  )
}



# ========================================================= #
# DEGgo colors
# ========================================================= #

#' Retrieve a DEGgo color palette
#'
#' Returns a qualitative, sequential, or diverging color palette suitable
#' for DEGgo visualizations and annotation variables.
#'
#' @param palette Character string specifying the palette name. Available
#'   palettes are `"default"`, `"nature"`, `"jama"`, `"nejm"`, `"lancet"`,
#'   `"viridis"`, and `"okabe"`.
#' @param n Optional positive integer giving the number of colors to return.
#'   When `NULL`, the complete palette is returned.
#' @param type Character string specifying the palette type. One of
#'   `"discrete"`, `"sequential"`, or `"diverging"`.
#' @param direction Integer equal to `1` or `-1`. Use `-1` to reverse the
#'   palette.
#'
#' @return A character vector of hexadecimal colors.
#'
#' @examples
#' deggo_palette("default")
#' deggo_palette("okabe", n = 4)
#' deggo_palette("viridis", n = 8, type = "sequential")
#' deggo_palette("default", n = 11, type = "diverging")
#'
#' @export
deggo_palette <- function(
    palette = c(
      "default",
      "nature",
      "jama",
      "nejm",
      "lancet",
      "viridis",
      "okabe"
    ),
    n = NULL,
    type = c("discrete", "sequential", "diverging"),
    direction = 1
) {

  palette <- match.arg(palette)
  type <- match.arg(type)

  if (
    length(direction) != 1L ||
    is.na(direction) ||
    !direction %in% c(-1, 1)
  ) {
    stop(
      "`direction` must be either 1 or -1.",
      call. = FALSE
    )
  }

  if (!is.null(n)) {

    if (
      length(n) != 1L ||
      is.na(n) ||
      !is.numeric(n) ||
      n < 1 ||
      n != as.integer(n)
    ) {
      stop(
        "`n` must be a positive integer or NULL.",
        call. = FALSE
      )
    }

    n <- as.integer(n)
  }

  palettes <- .deggo_palette_registry()

  cols <- palettes[[palette]][[type]]

  if (is.null(cols) || !length(cols)) {
    stop(
      sprintf(
        "Palette '%s' does not provide a '%s' color scale.",
        palette,
        type
      ),
      call. = FALSE
    )
  }

  cols <- unname(cols)

  if (direction == -1) {
    cols <- rev(cols)
  }

  if (is.null(n)) {
    return(cols)
  }

  .deggo_expand_palette(
    colors = cols,
    n = n,
    type = type
  )
}



#' @keywords internal
#' @noRd
.deggo_palette_registry <- function() {

  default_diverging <- c(
    "#2166AC",
    "#4393C3",
    "#92C5DE",
    "#D1E5F0",
    "#F7F7F7",
    "#FDDBC7",
    "#F4A582",
    "#D6604D",
    "#B2182B"
  )

  list(

    default = list(

      discrete = c(
        "#0072B2",
        "#D55E00",
        "#009E73",
        "#CC79A7",
        "#E69F00",
        "#56B4E9",
        "#882255",
        "#44AA99",
        "#AA4499",
        "#999999",
        "#000000",
        "#F0E442"
      ),

      sequential = c(
        "#F7FBFF",
        "#DEEBF7",
        "#C6DBEF",
        "#9ECAE1",
        "#6BAED6",
        "#4292C6",
        "#2171B5",
        "#08519C",
        "#08306B"
      ),

      diverging = default_diverging
    ),

    nature = list(

      discrete = c(
        "#3C5488",
        "#E64B35",
        "#00A087",
        "#4DBBD5",
        "#F39B7F",
        "#8491B4",
        "#91D1C2",
        "#DC0000",
        "#7E6148",
        "#B09C85"
      ),

      sequential = c(
        "#F1F5F9",
        "#D9E4EE",
        "#B7CDDE",
        "#8AAEC8",
        "#5D8DB3",
        "#3C6D9D",
        "#274F7C",
        "#17365D"
      ),

      diverging = c(
        "#3C5488",
        "#6686B7",
        "#A5B9D5",
        "#D9E1EC",
        "#F7F7F7",
        "#F4D4CD",
        "#EA9A8C",
        "#E06451",
        "#B52A1E"
      )
    ),

    jama = list(

      discrete = c(
        "#374E55",
        "#DF8F44",
        "#00A1D5",
        "#B24745",
        "#79AF97",
        "#6A6599",
        "#80796B"
      ),

      sequential = c(
        "#F4F7F7",
        "#DCE5E7",
        "#BACDD1",
        "#91AFB6",
        "#688F98",
        "#4B7079",
        "#374E55"
      ),

      diverging = c(
        "#374E55",
        "#63787E",
        "#93A3A7",
        "#C8D0D2",
        "#F7F7F7",
        "#F2D4BA",
        "#E9AF7A",
        "#DF8F44",
        "#B8651F"
      )
    ),

    nejm = list(

      discrete = c(
        "#BC3C29",
        "#0072B5",
        "#E18727",
        "#20854E",
        "#7876B1",
        "#6F99AD",
        "#FFDC91",
        "#EE4C97"
      ),

      sequential = c(
        "#F2F7FA",
        "#D4E7F0",
        "#ACD2E3",
        "#7EB7D2",
        "#4A98BE",
        "#207CAD",
        "#0072B5"
      ),

      diverging = c(
        "#0072B5",
        "#4C99C2",
        "#91BED5",
        "#C9DDE7",
        "#F7F7F7",
        "#EFCBC5",
        "#DD9186",
        "#C95C4C",
        "#BC3C29"
      )
    ),

    lancet = list(

      discrete = c(
        "#00468B",
        "#ED0000",
        "#42B540",
        "#0099B4",
        "#925E9F",
        "#FDAF91",
        "#AD002A",
        "#ADB6B6",
        "#1B1919"
      ),

      sequential = c(
        "#F2F7FA",
        "#D5E4EF",
        "#AAC9DF",
        "#73A9CC",
        "#3D88B5",
        "#146A9E",
        "#00468B"
      ),

      diverging = c(
        "#00468B",
        "#477EAC",
        "#89ACCA",
        "#C7D7E3",
        "#F7F7F7",
        "#F3C7C7",
        "#E98D8D",
        "#DD4C4C",
        "#C60000"
      )
    ),

    viridis = list(

      discrete = c(
        "#440154",
        "#482878",
        "#3E4989",
        "#31688E",
        "#26828E",
        "#1F9E89",
        "#35B779",
        "#6DCD59",
        "#B4DE2C",
        "#FDE725"
      ),

      sequential = c(
        "#440154",
        "#482878",
        "#3E4989",
        "#31688E",
        "#26828E",
        "#1F9E89",
        "#35B779",
        "#6DCD59",
        "#B4DE2C",
        "#FDE725"
      ),

      diverging = c(
        "#440154",
        "#3B528B",
        "#21918C",
        "#5EC962",
        "#FDE725",
        "#F4B942",
        "#E76F51",
        "#9D4EDD",
        "#5A189A"
      )
    ),

    okabe = list(

      discrete = c(
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#000000"
      ),

      sequential = c(
        "#F4F9FC",
        "#D7EBF5",
        "#ACD6EA",
        "#7CBCE0",
        "#56B4E9",
        "#278FC4",
        "#0072B2"
      ),

      diverging = c(
        "#0072B2",
        "#4E9BCD",
        "#91BFE0",
        "#C9DDED",
        "#F7F7F7",
        "#F4D6BA",
        "#E9A66F",
        "#DC7532",
        "#D55E00"
      )
    )
  )
}



#' @keywords internal
#' @noRd
.deggo_expand_palette <- function(
    colors,
    n,
    type = c("discrete", "sequential", "diverging")
) {

  type <- match.arg(type)

  colors <- unname(colors)

  if (!length(colors)) {
    stop(
      "`colors` must contain at least one color.",
      call. = FALSE
    )
  }

  if (n <= length(colors)) {
    return(colors[seq_len(n)])
  }

  if (type == "discrete") {

    if (n <= 12L) {
      return(rep(colors, length.out = n))
    }

    return(
      grDevices::hcl.colors(
        n = n,
        palette = "Dark 3"
      )
    )
  }

  grDevices::colorRampPalette(colors)(n)
}



#' @keywords internal
#' @noRd
.deggo_colors <- function(
    palette = "default"
) {

  list(

    discrete = deggo_palette(
      palette = palette,
      type = "discrete"
    ),

    sequential = deggo_palette(
      palette = palette,
      type = "sequential"
    ),

    diverging = deggo_palette(
      palette = palette,
      type = "diverging"
    ),

    heatmap = c(
      low  = "#6497B1",
      mid  = "#F7F7F7",
      high = "#740001"
    ),

    regulation = c(
      Up   = "#740001",
      Down = "#6497B1",
      NS   = "grey70"
    ),

    qc = c(
      PASS = "#009E73",
      WARN = "#E69F00",
      FAIL = "#D55E00"
    ),

    sex = c(
      Female = "#CC79A7",
      Male   = "#0072B2",
      "F"    = "#CC79A7",
      "M"    = "#0072B2"
    )
  )
}


#' Resolve Heatmap Annotation Columns
#'
#' @param metadata Sample metadata.
#' @param annotation_cols Character vector, `"auto"`, or `NULL`.
#' @param analysis_mode Analysis mode.
#' @param design_formula Optional model formula.
#' @param contrast Optional single-analysis contrast.
#' @param pairwise_group_cols Pairwise grouping columns.
#' @param current_contrast Current pairwise contrast.
#' @param sample_col Sample identifier column.
#'
#' @return Character vector of valid annotation columns.
#'
#' @keywords internal
#' @noRd
.deggo_heatmap_annotation_cols <- function(
    metadata,
    annotation_cols = "auto",
    analysis_mode = c("single", "pairwise"),
    design_formula = NULL,
    contrast = NULL,
    pairwise_group_cols = NULL,
    current_contrast = NULL,
    sample_col = "sample"
) {

  analysis_mode <- match.arg(analysis_mode)

  metadata <- as.data.frame(
    metadata,
    stringsAsFactors = FALSE
  )

  if (
    is.null(annotation_cols) ||
    identical(annotation_cols, FALSE)
  ) {
    return(character(0))
  }

  is_auto <- (
    length(annotation_cols) == 1L &&
      is.character(annotation_cols) &&
      identical(annotation_cols, "auto")
  )

  if (is_auto) {

    annotation_cols <- character(0)

    if (identical(analysis_mode, "single")) {

      formula_vars <- if (!is.null(design_formula)) {
        all.vars(design_formula)
      } else {
        character(0)
      }

      contrast_var <- if (
        !is.null(contrast) &&
        length(contrast) >= 1L
      ) {
        contrast[1]
      } else {
        character(0)
      }

      annotation_cols <- c(
        formula_vars,
        contrast_var
      )
    }

    if (identical(analysis_mode, "pairwise")) {

      current_contrast_var <- if (
        !is.null(current_contrast) &&
        length(current_contrast) >= 1L
      ) {
        current_contrast[1]
      } else {
        character(0)
      }

      annotation_cols <- c(
        pairwise_group_cols,
        current_contrast_var
      )
    }
  }

  annotation_cols <- unique(
    as.character(annotation_cols)
  )

  annotation_cols <- annotation_cols[
    !is.na(annotation_cols) &
      nzchar(annotation_cols)
  ]

  annotation_cols <- setdiff(
    annotation_cols,
    sample_col
  )

  annotation_cols <- intersect(
    annotation_cols,
    colnames(metadata)
  )

  if (!length(annotation_cols)) {
    return(character(0))
  }

  annotation_cols[
    vapply(
      annotation_cols,
      function(v) {
        values <- metadata[[v]]
        values <- values[
          !is.na(values) &
            nzchar(as.character(values))
        ]

        length(unique(values)) >= 1L
      },
      logical(1)
    )
  ]
}


#' Generate colors for heatmap annotations
#'
#' Internal helper used to assign deterministic colors to categorical
#' annotation variables. Known semantic variables, such as sex or gender,
#' use predefined colors when available. Other variables use the general
#' discrete palette.
#'
#' @param annotation_col A data frame containing categorical annotation
#'   variables.
#' @param palettes A named list of color palettes, usually returned by
#'   `.deggo_colors()`.
#'
#' @return A named list of named color vectors, or `NULL` when no valid
#'   annotation variable is available.
#' @keywords internal
#' @noRd
.deggo_annotation_colors <- function(
    annotation_col,
    palettes = NULL,
    palette = "default"
) {

  if (
    is.null(annotation_col) ||
    !is.data.frame(annotation_col) ||
    ncol(annotation_col) == 0L
  ) {
    return(NULL)
  }

  if (is.null(palettes)) {
    palettes <- .deggo_colors(palette = palette)
  }

  discrete <- palettes$discrete

  if (is.null(discrete) || !length(discrete)) {
    discrete <- grDevices::hcl.colors(
      n = 12,
      palette = "Dark 3"
    )
  }

  semantic_palettes <- list()

  if (!is.null(palettes$sex)) {
    semantic_palettes$sex <- palettes$sex
    semantic_palettes$gender <- palettes$sex
  }

  out <- vector(
    mode = "list",
    length = ncol(annotation_col)
  )

  names(out) <- colnames(annotation_col)

  for (nm in colnames(annotation_col)) {

    x <- annotation_col[[nm]]

    if (is.factor(x)) {

      observed <- as.character(x)

      levels_use <- levels(x)

      levels_use <- levels_use[
        levels_use %in%
          observed[!is.na(observed) & nzchar(observed)]
      ]

    } else {

      values <- as.character(x)

      levels_use <- unique(
        values[!is.na(values) & nzchar(values)]
      )
    }

    if (!length(levels_use)) {
      out[[nm]] <- NULL
      next
    }

    nm_key <- tolower(trimws(nm))
    known <- semantic_palettes[[nm_key]]

    if (
      !is.null(known) &&
      length(known) &&
      !is.null(names(known))
    ) {

      known <- known[
        !is.na(names(known)) &
          nzchar(names(known))
      ]

      known <- known[
        !duplicated(names(known))
      ]

      matched <- intersect(
        levels_use,
        names(known)
      )

      missing <- setdiff(
        levels_use,
        names(known)
      )

      cols <- known[matched]

      if (length(missing)) {

        available <- discrete[
          !discrete %in% unname(cols)
        ]

        if (!length(available)) {
          available <- discrete
        }

        extra <- stats::setNames(
          rep(
            available,
            length.out = length(missing)
          ),
          missing
        )

        cols <- c(cols, extra)
      }

      out[[nm]] <- cols[levels_use]

    } else {

      out[[nm]] <- stats::setNames(
        deggo_palette(
          palette = palette,
          n = length(levels_use),
          type = "discrete"
        ),
        levels_use
      )
    }
  }

  out <- out[
    !vapply(out, is.null, logical(1))
  ]

  if (!length(out)) {
    return(NULL)
  }

  out
}



#' @keywords internal
#' @noRd
.deggo_label_col <- function(df) {
  if ("feature" %in% colnames(df)) {
    "feature"
  } else if ("symbol" %in% colnames(df)) {
    "symbol"
  } else if ("SYMBOL" %in% colnames(df)) {
    "SYMBOL"
  } else if ("gene_name" %in% colnames(df)) {
    "gene_name"
  } else if ("ENSEMBL" %in% colnames(df)) {
    "ENSEMBL"
  } else if ("gene_id" %in% colnames(df)) {
    "gene_id"
  } else {
    stop("No gene label column found.", call. = FALSE)
  }
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
#' @param style Theme style passed to \code{.deggo_theme()}.
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
    style = "classic",
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
    .deggo_theme(style = style, txtsize = txtsize, x.ang = x_ang)

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
#' @param palette Optional color palette used for regulation classes. It may
#'   be a named DEGgo palette or a named character vector defining colors for
#'   `"Up"`, `"Down"`, and `"NS"`.
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
#' @param label_size Size of standard gene labels.
#' @param highlight_label_size Size of highlighted gene labels.
#' @param label_bg Fill color for label backgrounds.
#' @param label_bg_r Radius of label background boxes.
#' @param label_force Force parameter passed to ggrepel.
#' @param label_force_pull Pull force parameter passed to ggrepel.
#' @param label_max_overlaps Maximum number of overlapping labels.
#' @param seed Random seed for reproducible label placement.
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
    palette = "default",
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
    style = c("classic", "minimal", "bw", "test", "void", "dirty", "gray"),
    txtsize = 12,
    label_size = 3,
    highlight_label_size = 3.5,
    label_bg = "grey95",
    label_bg_r = 0.12,
    label_force = 1.2,
    label_force_pull = 0.4,
    label_max_overlaps = Inf,
    seed = 42
) {

  style <- match.arg(style)

  log <- .deggo_msg(verbose = TRUE, prefix = "DEGgo")
  log("Generating volcano plot...", type = "info")

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required.", call. = FALSE)
  }

  res_df <- as.data.frame(res_df, stringsAsFactors = FALSE)

  if (is.null(gene_col)) {
    gene_col <- .deggo_label_col(res_df)
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

  res_df[[lfc_col]] <- suppressWarnings(as.numeric(res_df[[lfc_col]]))
  res_df[[p_col]] <- suppressWarnings(as.numeric(res_df[[p_col]]))

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
    levels = c("Up", "Down", "Ns")
  )

  res_df$neg_log10_p <- -log10(res_df[[p_col]])

  counts <- table(res_df$volcano_status)

  legend_labels <- c(
    "Up" = paste0("Up (", counts[["Up"]], ")"),
    "Down" = paste0("Down (", counts[["Down"]], ")"),
    "Ns" = paste0("Ns (", counts[["Ns"]], ")")
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
    top_genes$label_score <-
      abs(top_genes[[lfc_col]]) *
      top_genes$neg_log10_p

    top_genes <- top_genes[
      order(top_genes$label_score, decreasing = TRUE),
      ,
      drop = FALSE
    ]

    top_genes <- utils::head(top_genes, top_n_labels)
  } else {
    top_genes <- top_genes[0, , drop = FALSE]
  }

  if (!is.null(genes_highlight) && nrow(top_genes) > 0) {
    top_genes <- top_genes[
      !(top_genes[[gene_col]] %in% genes_highlight),
      ,
      drop = FALSE
    ]
  }

  p <- ggplot2::ggplot(
    res_df,
    ggplot2::aes(
      x = .data[[lfc_col]],
      y = .data[["neg_log10_p"]]
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[["volcano_status"]]),
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
      xintercept = c(-logfc_cutoff, logfc_cutoff),
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.4
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-max_abs, max_abs),
      clip = "off"
    ) +
    ggplot2::labs(
      title = title,
      x = expression(log[2]~Fold~Change),
      y = if (isTRUE(use_padj)) {
        expression(-log[10]~adjusted~italic(P))
      } else {
        expression(-log[10]~italic(P))
      }
    ) +
    .deggo_theme(
      style = style,
      txtsize = txtsize
    )

  if (nrow(top_genes) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = top_genes,
        ggplot2::aes(label = .data[[gene_col]]),
        color = "black",
        size = label_size,
        fontface = "plain",
        bg.color = label_bg,
        bg.r = label_bg_r,
        box.padding = 0.45,
        point.padding = 0.25,
        min.segment.length = 0,
        segment.color = "grey55",
        max.overlaps = label_max_overlaps,
        force = label_force,
        force_pull = label_force_pull,
        seed = seed
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
          ggplot2::aes(
            fill = .data[["volcano_status"]]
          ),
          shape = 21,
          size = point_size + 0.5,
          color = "black",
          stroke = 0.2
        ) +
        ggplot2::scale_fill_manual(
          values = colors,
          guide = "none",
          drop = FALSE
        ) +
        ggrepel::geom_text_repel(
          data = highlight_df,
          ggplot2::aes(label = .data[[gene_col]]),
          color = "black",
          size = highlight_label_size,
          fontface = "bold",
          bg.color = label_bg,
          bg.r = label_bg_r,
          box.padding = 0.65,
          point.padding = 0.35,
          min.segment.length = 0,
          segment.color = "grey35",
          max.overlaps = Inf,
          force = label_force + 0.8,
          force_pull = label_force_pull,
          seed = seed
        )
    }
  }

  if (!isTRUE(legend)) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    ggplot2::ggsave(
      filename = file.path(output_dir, paste0(filename, ".png")),
      plot = p,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )

    if (isTRUE(save_pdf)) {
      ggplot2::ggsave(
        filename = file.path(output_dir, paste0(filename, ".pdf")),
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
#' @param style Plot style passed to \code{.deggo_theme()}.
#' @param txtsize Base text size.
#' @param palette Optional color palette. Can be a named palette supported by
#'   DEGgo or a character vector of colors.
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
    dpi = 300,
    style = "classic",
    txtsize = 12,
    palette = "default"
) {

  log <- .deggo_msg(verbose = TRUE, prefix = "DEGgo")
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
    .deggo_theme(style = style, txtsize = txtsize) +
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

#' Compute Automatic Heatmap Dimensions
#'
#' Internal helper used by DEGgo to compute publication-friendly heatmap
#' dimensions from the number of genes and samples.
#'
#' @param n_genes Integer. Number of genes displayed in the heatmap.
#' @param n_samples Integer. Number of samples displayed in the heatmap.
#' @param base_width Numeric. Base plot width in inches.
#' @param base_height Numeric. Base plot height in inches.
#' @param gene_unit Numeric. Height added per gene.
#' @param sample_unit Numeric. Width added per sample.
#' @param min_width,max_width Numeric. Minimum and maximum plot width.
#' @param min_height,max_height Numeric. Minimum and maximum plot height.
#'
#' @return A list with `width` and `height`.
#'
#' @keywords internal
#' @noRd
.deggo_heatmap_size <- function(
    n_genes,
    n_samples,
    base_width = 5.5,
    base_height = 4.5,
    gene_unit = 0.18,
    sample_unit = 0.22,
    min_width = 6,
    max_width = 14,
    min_height = 5,
    max_height = 16
) {
  width <- base_width + n_samples * sample_unit
  height <- base_height + n_genes * gene_unit

  width <- max(min_width, min(width, max_width))
  height <- max(min_height, min(height, max_height))

  list(width = width, height = height)
}


# ========================================================= #
# plot_heatmap
# ========================================================= #

#' @keywords internal
#' @noRd
.prepare_heatmap_annotations <- function(
    metadata,
    annotation_cols,
    annotation_colors = NULL
){

  annotation <- metadata[, intersect(annotation_cols,colnames(metadata)), drop=FALSE]

  if(ncol(annotation)){
    annotation[] <- lapply(annotation,factor)
  }else{
    annotation <- NULL
  }

  if(is.null(annotation_colors))
    annotation_colors <- .deggo_annotation_colors(annotation)
  list(
    annotation=annotation,
    colors=annotation_colors
  )
}




#' Prepare DEGgo heatmap matrix
#'
#' Internal helper used by `plot_heatmap()` to extract the VST matrix,
#' subset samples, select top DEG genes and optionally scale rows.
#'
#' @param vsd Variance-stabilized object, usually returned by DESeq2::vst().
#' @param res_df Differential expression result table.
#' @param metadata Sample metadata data frame.
#' @param contrast Optional contrast vector used to subset samples.
#' @param sample_subset Optional character vector of samples to keep.
#' @param metadata_filter Optional named list used to filter metadata rows.
#' @param top_n_heatmap Number of top genes to display.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param fallback Logical. If TRUE, use top ranked genes when no significant
#'   genes pass `padj_cutoff`.
#' @param order_by Optional metadata columns used to order samples.
#' @param scale_rows Logical. If TRUE, z-score scale genes by row.
#' @param log Optional DEGgo logger.
#'
#' @return A list containing `matrix`, `metadata` and `top`, or NULL.
#'
#' @keywords internal
#' @noRd
.prepare_heatmap_matrix <- function(
    vsd,
    res_df,
    metadata,
    contrast = NULL,
    sample_subset = NULL,
    metadata_filter = NULL,
    top_n_heatmap = 20,
    padj_cutoff = 0.05,
    fallback = TRUE,
    order_by = NULL,
    scale_rows = TRUE,
    log = NULL
) {

  .return_null <- function(reason) {
    if (!is.null(log)) log(reason, type = "warn")
    return(NULL)
  }

  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required.", call. = FALSE)
  }

  if (is.null(vsd) || is.null(res_df) || is.null(metadata)) {
    return(.return_null("Heatmap skipped: vsd, res_df or metadata is NULL."))
  }

  mat <- SummarizedExperiment::assay(vsd)

  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    return(.return_null("Heatmap skipped: vst matrix must have rownames and colnames."))
  }

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!"sample" %in% colnames(metadata)) {
    metadata$sample <- rownames(metadata)
  }

  metadata$sample <- as.character(metadata$sample)
  rownames(metadata) <- metadata$sample

  common_samples <- intersect(colnames(mat), rownames(metadata))

  if (length(common_samples) < 2) {
    return(.return_null("Heatmap skipped: fewer than 2 common samples between vst matrix and metadata."))
  }

  mat <- mat[, common_samples, drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

  if (!is.null(sample_subset)) {

    sample_subset <- intersect(sample_subset, colnames(mat))

    if (length(sample_subset) < 2) {
      return(.return_null("Heatmap skipped: fewer than 2 samples after sample_subset."))
    }

    mat <- mat[, sample_subset, drop = FALSE]
    metadata <- metadata[sample_subset, , drop = FALSE]
  }

  if (!is.null(metadata_filter)) {

    keep <- rep(TRUE, nrow(metadata))

    for (nm in names(metadata_filter)) {

      if (!nm %in% colnames(metadata)) {
        stop("metadata_filter column not found: ", nm, call. = FALSE)
      }

      keep <- keep & metadata[[nm]] %in% metadata_filter[[nm]]
    }

    metadata <- metadata[keep, , drop = FALSE]

    if (nrow(metadata) < 2) {
      return(.return_null("Heatmap skipped: fewer than 2 samples after metadata_filter."))
    }

    mat <- mat[, rownames(metadata), drop = FALSE]
  }

  if (!is.null(contrast)) {

    contrast_col <- contrast[1]
    contrast_levels <- contrast[2:3]

    if (!contrast_col %in% colnames(metadata)) {
      stop("Contrast column not found in metadata: ", contrast_col, call. = FALSE)
    }

    keep_samples <- rownames(metadata)[metadata[[contrast_col]] %in% contrast_levels]
    keep_samples <- intersect(keep_samples, colnames(mat))

    if (length(keep_samples) < 2) {
      return(.return_null(
        paste0(
          "Heatmap skipped: fewer than 2 samples found for contrast: ",
          paste(contrast, collapse = " ")
        )
      ))
    }

    mat <- mat[, keep_samples, drop = FALSE]
    metadata <- metadata[keep_samples, , drop = FALSE]
  }

  required_cols <- c("padj", "log2FoldChange")
  missing_cols <- setdiff(required_cols, colnames(res_df))

  if (length(missing_cols)) {
    stop(
      "res_df missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!"ENSEMBL" %in% colnames(res_df)) {
    if ("gene_id" %in% colnames(res_df)) {
      res_df$ENSEMBL <- res_df$gene_id
    } else {
      stop("res_df must contain an ENSEMBL or gene_id column.", call. = FALSE)
    }
  }

  label_col <- .deggo_label_col(res_df)

  rownames(mat) <- sub("\\..*$", "", rownames(mat))
  res_df$ENSEMBL <- sub("\\..*$", "", as.character(res_df$ENSEMBL))

  top <- res_df[
    !is.na(res_df$padj) &
      !is.na(res_df$log2FoldChange) &
      res_df$padj < padj_cutoff,
    ,
    drop = FALSE
  ]

  top <- top[
    order(top$padj, -abs(top$log2FoldChange)),
    ,
    drop = FALSE
  ]

  top <- utils::head(top, top_n_heatmap)

  if (!nrow(top) && isTRUE(fallback)) {

    if (!is.null(log)) {
      log("No significant genes found. Using fallback genes ranked by padj.", type = "warn")
    }

    top <- res_df[
      !is.na(res_df$padj) &
        !is.na(res_df$log2FoldChange),
      ,
      drop = FALSE
    ]

    top <- top[
      order(top$padj, -abs(top$log2FoldChange)),
      ,
      drop = FALSE
    ]

    top <- utils::head(top, top_n_heatmap)
  }

  if (!nrow(top)) {
    return(.return_null("Heatmap skipped: no genes available after filtering."))
  }

  genes_use <- intersect(top$ENSEMBL, rownames(mat))

  if (!length(genes_use)) {
    return(.return_null("Heatmap skipped: no top genes found in vst matrix rownames."))
  }

  top <- top[
    match(genes_use, top$ENSEMBL),
    ,
    drop = FALSE
  ]

  mat_use <- mat[genes_use, , drop = FALSE]

  keep_rows <- apply(mat_use, 1, function(x) all(is.finite(x)))
  keep_cols <- apply(mat_use, 2, function(x) all(is.finite(x)))

  mat_use <- mat_use[keep_rows, keep_cols, drop = FALSE]

  if (nrow(mat_use) < 2 || ncol(mat_use) < 2) {
    return(.return_null("Heatmap skipped: fewer than 2 genes or 2 samples after finite-value filtering."))
  }

  top <- top[
    match(rownames(mat_use), top$ENSEMBL),
    ,
    drop = FALSE
  ]

  gene_labels <- as.character(top[[label_col]])

  gene_labels[is.na(gene_labels) | gene_labels == ""] <- top$ENSEMBL[
    is.na(gene_labels) | gene_labels == ""
  ]

  rownames(mat_use) <- make.unique(gene_labels)

  metadata_use <- metadata[colnames(mat_use), , drop = FALSE]

  if (!is.null(order_by)) {

    missing <- setdiff(order_by, colnames(metadata_use))

    if (length(missing)) {
      stop(
        "order_by column(s) not found: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }

    ord <- do.call(order, metadata_use[, order_by, drop = FALSE])
    metadata_use <- metadata_use[ord, , drop = FALSE]
    mat_use <- mat_use[, rownames(metadata_use), drop = FALSE]
  }

  if (isTRUE(scale_rows)) {

    gene_sd <- apply(mat_use, 1, stats::sd, na.rm = TRUE)
    keep_genes <- is.finite(gene_sd) & gene_sd > 0

    mat_use <- mat_use[keep_genes, , drop = FALSE]

    if (nrow(mat_use) < 2) {
      return(.return_null("Heatmap skipped: fewer than 2 genes after removing zero-variance genes."))
    }

    mat_use <- t(scale(t(mat_use)))

    keep_finite <- apply(mat_use, 1, function(x) all(is.finite(x)))
    mat_use <- mat_use[keep_finite, , drop = FALSE]
  }

  if (nrow(mat_use) < 2 || ncol(mat_use) < 2) {
    return(.return_null("Heatmap skipped: final matrix has fewer than 2 genes or 2 samples."))
  }

  list(
    matrix = mat_use,
    metadata = metadata_use,
    top = top
  )
}


#' Save DEGgo heatmap
#'
#' Internal helper wrapping `pheatmap::pheatmap()`.
#'
#' @param mat Numeric heatmap matrix.
#' @param annotation_col Column annotation data frame or NULL.
#' @param annotation_colors Annotation colors list or NULL.
#' @param main Heatmap title.
#' @param output_dir Output directory.
#' @param filename Output filename without extension.
#' @param cluster_rows Logical.
#' @param cluster_cols Logical.
#' @param fontsize_row Row-label font size.
#' @param fontsize_col Column-label font size.
#' @param width Optional output width in inches.
#' @param height Optional output height in inches.
#' @param show_rownames Logical or NULL.
#' @param show_colnames Logical or NULL.
#'
#' @return Invisibly returns the output file path.
#'
#' @keywords internal
#' @noRd
.save_heatmap <- function(
    mat,
    annotation_col = NULL,
    annotation_colors = NULL,
    main = "Top Differentially Expressed Genes",
    output_dir = "DEGgo_out",
    filename = "Heatmap",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = NULL,
    fontsize_col = NULL,
    show_dendrogram = FALSE,
    width = NULL,
    height = NULL,
    show_rownames = NULL,
    show_colnames = NULL
) {

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  hm_size <- .deggo_heatmap_size(
    n_genes = nrow(mat),
    n_samples = ncol(mat)
  )

  if (is.null(width)) {
    width <- hm_size$width
  }

  if (is.null(height)) {
    height <- hm_size$height
  }

  fontsize_row_use <- fontsize_row %||% if (nrow(mat) <= 20) {
    9
  } else if (nrow(mat) <= 40) {
    7
  } else if (nrow(mat) <= 80) {
    6
  } else {
    5
  }

  fontsize_col_use <- fontsize_col %||% if (ncol(mat) <= 12) {
    9
  } else if (ncol(mat) <= 25) {
    7
  } else if (ncol(mat) <= 50) {
    6
  } else {
    5
  }

  if (is.null(show_rownames)) {
    show_rownames <- nrow(mat) <= 80
  }

  if (is.null(show_colnames)) {
    show_colnames <- ncol(mat) <= 60
  }

  outfile <- file.path(output_dir, paste0(filename, ".png"))

  pheatmap::pheatmap(
    mat,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    fontsize_row = fontsize_row_use,
    fontsize_col = fontsize_col_use,
    fontsize = 9,
    border_color = NA,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    treeheight_row = if (isTRUE(cluster_rows) && isTRUE(show_dendrogram)) 50 else 0,
    treeheight_col = if (isTRUE(cluster_cols) && isTRUE(show_dendrogram)) 50 else 0,
    angle_col = 90,
    main = main,
    color = grDevices::colorRampPalette(
      c("#6497B1", "#F7F7F7", "#740001")
    )(100),
    breaks = seq(-2, 2, length.out = 101),
    filename = outfile,
    width = width,
    height = height
  )

  if (!file.exists(outfile)) {
    warning("Heatmap file was not created by pheatmap.", call. = FALSE)
    return(invisible(NULL))
  }

  invisible(outfile)
}



#' Plot a DEGgo Differential Expression Heatmap
#'
#' @param vsd Variance-stabilized expression object.
#' @param res_df Differential-expression table, preferably a cleaned table
#'   from `de_results$sig_deg_clean`.
#' @param metadata Sample metadata.
#' @param contrast Optional contrast vector.
#' @param sample_subset Optional samples to retain.
#' @param metadata_filter Optional named metadata filtering list.
#' @param top_n_heatmap Number of genes displayed.
#' @param padj_cutoff Adjusted P-value cutoff.
#' @param logfc_cutoff Absolute log2 fold-change cutoff used when selecting
#'   significant genes for the heatmap.
#' @param main Heatmap title.
#' @param output_dir Output directory.
#' @param filename Output filename without extension.
#' @param fallback Use ranked genes if no gene passes the significance cutoff.
#' @param annotation_colors Optional annotation-color list.
#' @param annotation_cols Metadata annotation columns.
#' @param order_by Metadata variables used to order samples.
#' @param order_levels Optional named list defining sample-group order.
#' @param scale_rows Scale expression by gene.
#' @param cluster_rows Cluster genes.
#' @param cluster_cols Cluster samples.
#' @param fontsize_row Row-label size.
#' @param fontsize_col Column-label size.
#' @param show_dendrogram Display dendrograms.
#' @param width Output width.
#' @param height Output height.
#' @param show_rownames Display row names.
#' @param show_colnames Display column names.
#' @param palette DEGgo palette name.
#'
#' @return Invisibly returns the plotted expression matrix.
#'
#' @export
plot_heatmap <- function(
    vsd,
    res_df,
    metadata,
    contrast = NULL,
    sample_subset = NULL,
    metadata_filter = NULL,
    top_n_heatmap = 20,
    padj_cutoff = 0.05,
    logfc_cutoff = 0.25,
    main = "Top Differentially Expressed Genes",
    output_dir = "DEGgo_out",
    filename = "Heatmap",
    fallback = TRUE,
    annotation_colors = NULL,
    annotation_cols = NULL,
    order_by = NULL,
    order_levels = NULL,
    scale_rows = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = NULL,
    fontsize_col = NULL,
    show_dendrogram = FALSE,
    width = NULL,
    height = NULL,
    show_rownames = NULL,
    show_colnames = NULL,
    palette = "default"
) {

  log <- .deggo_msg(
    verbose = TRUE,
    prefix = "DEGgo"
  )

  log(
    "Generating heatmap...",
    type = "info"
  )

  .return_null <- function(reason) {
    log(reason, type = "warn")
    return(NULL)
  }

  # ---------------------------------------------------------- #
  # Validate packages and arguments
  # ---------------------------------------------------------- #

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop(
      "Package 'pheatmap' is required.",
      call. = FALSE
    )
  }

  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop(
      "Package 'SummarizedExperiment' is required.",
      call. = FALSE
    )
  }

  palette <- match.arg(
    palette,
    choices = c(
      "default",
      "nature",
      "jama",
      "nejm",
      "lancet",
      "viridis",
      "okabe"
    )
  )

  if (
    is.null(vsd) ||
    is.null(res_df) ||
    is.null(metadata)
  ) {
    return(
      .return_null(
        "Heatmap skipped: vsd, res_df or metadata is NULL."
      )
    )
  }

  res_df <- as.data.frame(
    res_df,
    stringsAsFactors = FALSE
  )

  if (!nrow(res_df)) {
    return(
      .return_null(
        "Heatmap skipped: the DEG table is empty."
      )
    )
  }

  # ---------------------------------------------------------- #
  # Prepare expression matrix and metadata
  # ---------------------------------------------------------- #

  mat <- SummarizedExperiment::assay(vsd)

  if (
    is.null(rownames(mat)) ||
    is.null(colnames(mat))
  ) {
    return(
      .return_null(
        "Heatmap skipped: VST matrix must have row and column names."
      )
    )
  }

  metadata <- as.data.frame(
    metadata,
    stringsAsFactors = FALSE
  )

  if (!"sample" %in% colnames(metadata)) {
    metadata$sample <- rownames(metadata)
  }

  metadata$sample <- as.character(metadata$sample)
  rownames(metadata) <- metadata$sample

  common_samples <- intersect(
    colnames(mat),
    rownames(metadata)
  )

  if (length(common_samples) < 2L) {
    return(
      .return_null(
        paste0(
          "Heatmap skipped: fewer than 2 common samples ",
          "between VST matrix and metadata."
        )
      )
    )
  }

  mat <- mat[, common_samples, drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

  # ---------------------------------------------------------- #
  # Sample subset
  # ---------------------------------------------------------- #

  if (!is.null(sample_subset)) {

    sample_subset <- intersect(
      as.character(sample_subset),
      colnames(mat)
    )

    if (length(sample_subset) < 2L) {
      return(
        .return_null(
          "Heatmap skipped: fewer than 2 samples after sample subsetting."
        )
      )
    }

    mat <- mat[, sample_subset, drop = FALSE]
    metadata <- metadata[sample_subset, , drop = FALSE]
  }

  # ---------------------------------------------------------- #
  # Metadata filtering
  # ---------------------------------------------------------- #

  if (!is.null(metadata_filter)) {

    if (
      !is.list(metadata_filter) ||
      is.null(names(metadata_filter))
    ) {
      stop(
        "`metadata_filter` must be a named list.",
        call. = FALSE
      )
    }

    keep <- rep(TRUE, nrow(metadata))

    for (nm in names(metadata_filter)) {

      if (!nm %in% colnames(metadata)) {
        stop(
          "Metadata filter column not found: ",
          nm,
          call. = FALSE
        )
      }

      keep <- keep &
        metadata[[nm]] %in% metadata_filter[[nm]]
    }

    metadata <- metadata[keep, , drop = FALSE]

    if (nrow(metadata) < 2L) {
      return(
        .return_null(
          "Heatmap skipped: fewer than 2 samples after metadata filtering."
        )
      )
    }

    mat <- mat[, rownames(metadata), drop = FALSE]
  }

  # ---------------------------------------------------------- #
  # Contrast-specific sample subset
  # ---------------------------------------------------------- #

  if (!is.null(contrast)) {

    if (length(contrast) < 3L) {
      stop(
        "`contrast` must contain a variable and two levels.",
        call. = FALSE
      )
    }

    contrast_col <- contrast[1]
    contrast_levels <- contrast[2:3]

    if (!contrast_col %in% colnames(metadata)) {
      stop(
        "Contrast column not found in metadata: ",
        contrast_col,
        call. = FALSE
      )
    }

    keep_samples <- rownames(metadata)[
      !is.na(metadata[[contrast_col]]) &
        metadata[[contrast_col]] %in% contrast_levels
    ]

    keep_samples <- intersect(
      keep_samples,
      colnames(mat)
    )

    if (length(keep_samples) < 2L) {
      return(
        .return_null(
          paste0(
            "Heatmap skipped: fewer than 2 samples for contrast ",
            paste(contrast_levels, collapse = " vs "),
            "."
          )
        )
      )
    }

    mat <- mat[, keep_samples, drop = FALSE]
    metadata <- metadata[keep_samples, , drop = FALSE]
  }

  # ---------------------------------------------------------- #
  # Validate DEG columns
  # ---------------------------------------------------------- #

  required_cols <- c(
    "padj",
    "log2FoldChange"
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(res_df)
  )

  if (length(missing_cols)) {
    stop(
      "DEG table is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!"ENSEMBL" %in% colnames(res_df)) {

    if ("gene_id" %in% colnames(res_df)) {
      res_df$ENSEMBL <- res_df$gene_id
    } else {
      stop(
        "DEG table must contain `ENSEMBL` or `gene_id`.",
        call. = FALSE
      )
    }
  }

  label_col <- .deggo_label_col(res_df)

  rownames(mat) <- sub(
    "\\..*$",
    "",
    rownames(mat)
  )

  res_df$ENSEMBL <- sub(
    "\\..*$",
    "",
    as.character(res_df$ENSEMBL)
  )

  # ---------------------------------------------------------- #
  # Select top genes
  # ---------------------------------------------------------- #

  top <- res_df[
    !is.na(res_df$padj) &
      !is.na(res_df$log2FoldChange) &
      res_df$padj < padj_cutoff &
      abs(res_df$log2FoldChange) >= logfc_cutoff,
    ,
    drop = FALSE
  ]

  top <- top[
    order(
      top$padj,
      -abs(top$log2FoldChange)
    ),
    ,
    drop = FALSE
  ]

  top <- utils::head(
    top,
    top_n_heatmap
  )

  if (!nrow(top) && isTRUE(fallback)) {

    log(
      paste0(
        "No significant genes found. ",
        "Using fallback genes ranked by adjusted P-value."
      ),
      type = "warn"
    )

    top <- res_df[
      !is.na(res_df$padj) &
        !is.na(res_df$log2FoldChange),
      ,
      drop = FALSE
    ]

    top <- top[
      order(
        top$padj,
        -abs(top$log2FoldChange)
      ),
      ,
      drop = FALSE
    ]

    top <- utils::head(
      top,
      top_n_heatmap
    )
  }

  if (!nrow(top)) {
    return(
      .return_null(
        "Heatmap skipped: no genes available after selection."
      )
    )
  }

  genes_use <- intersect(
    top$ENSEMBL,
    rownames(mat)
  )

  if (length(genes_use) < 2L) {
    return(
      .return_null(
        "Heatmap skipped: fewer than 2 selected genes found in the VST matrix."
      )
    )
  }

  top <- top[
    match(genes_use, top$ENSEMBL),
    ,
    drop = FALSE
  ]

  mat_use <- mat[
    genes_use,
    ,
    drop = FALSE
  ]

  # ---------------------------------------------------------- #
  # Remove non-finite rows and columns
  # ---------------------------------------------------------- #

  keep_rows <- apply(
    mat_use,
    1,
    function(x) all(is.finite(x))
  )

  keep_cols <- apply(
    mat_use,
    2,
    function(x) all(is.finite(x))
  )

  mat_use <- mat_use[
    keep_rows,
    keep_cols,
    drop = FALSE
  ]

  if (
    nrow(mat_use) < 2L ||
    ncol(mat_use) < 2L
  ) {
    return(
      .return_null(
        paste0(
          "Heatmap skipped: fewer than 2 genes or samples ",
          "after finite-value filtering."
        )
      )
    )
  }

  top <- top[
    match(rownames(mat_use), top$ENSEMBL),
    ,
    drop = FALSE
  ]

  gene_labels <- as.character(
    top[[label_col]]
  )

  missing_labels <-
    is.na(gene_labels) |
    !nzchar(gene_labels)

  gene_labels[missing_labels] <-
    top$ENSEMBL[missing_labels]

  rownames(mat_use) <- make.unique(
    gene_labels
  )

  metadata_use <- metadata[
    colnames(mat_use),
    ,
    drop = FALSE
  ]

  # ---------------------------------------------------------- #
  # Sample ordering
  # ---------------------------------------------------------- #

  if (
    !is.null(order_by) &&
    length(order_by)
  ) {

    order_by <- intersect(
      as.character(order_by),
      colnames(metadata_use)
    )

    if (length(order_by)) {

      for (v in order_by) {

        values <- as.character(
          metadata_use[[v]]
        )

        if (
          !is.null(order_levels) &&
          v %in% names(order_levels)
        ) {
          metadata_use[[v]] <- factor(
            values,
            levels = order_levels[[v]]
          )
        } else {
          metadata_use[[v]] <- factor(
            values,
            levels = unique(values)
          )
        }
      }

      ord <- do.call(
        order,
        metadata_use[, order_by, drop = FALSE]
      )

      metadata_use <- metadata_use[
        ord,
        ,
        drop = FALSE
      ]

      mat_use <- mat_use[
        ,
        rownames(metadata_use),
        drop = FALSE
      ]
    }
  }

  # ---------------------------------------------------------- #
  # Annotation
  # ---------------------------------------------------------- #

  annotation_cols <- as.character(
    annotation_cols %||% character(0)
  )

  annotation_cols <- intersect(
    annotation_cols,
    colnames(metadata_use)
  )

  if (length(annotation_cols)) {

    annotation_col <- metadata_use[
      ,
      annotation_cols,
      drop = FALSE
    ]

    for (v in colnames(annotation_col)) {

      values <- as.character(
        annotation_col[[v]]
      )

      annotation_col[[v]] <- factor(
        values,
        levels = unique(values)
      )
    }

  } else {

    annotation_col <- NULL
  }

  # ---------------------------------------------------------- #
  # Row scaling
  # ---------------------------------------------------------- #

  if (isTRUE(scale_rows)) {

    gene_sd <- apply(
      mat_use,
      1,
      stats::sd,
      na.rm = TRUE
    )

    keep_genes <-
      is.finite(gene_sd) &
      gene_sd > 0

    mat_use <- mat_use[
      keep_genes,
      ,
      drop = FALSE
    ]

    if (nrow(mat_use) < 2L) {
      return(
        .return_null(
          "Heatmap skipped: fewer than 2 nonconstant genes."
        )
      )
    }

    mat_use <- t(
      scale(
        t(mat_use)
      )
    )

    mat_use <- mat_use[
      apply(
        mat_use,
        1,
        function(x) all(is.finite(x))
      ),
      ,
      drop = FALSE
    ]
  }

  if (
    nrow(mat_use) < 2L ||
    ncol(mat_use) < 2L
  ) {
    return(
      .return_null(
        "Heatmap skipped: final matrix is too small."
      )
    )
  }

  # ---------------------------------------------------------- #
  # Colors
  # ---------------------------------------------------------- #

  palettes <- .deggo_colors(
    palette = palette
  )

  if (
    is.null(annotation_colors) ||
    is.character(annotation_colors)
  ) {
    annotation_colors <- .deggo_annotation_colors(
      annotation_col = annotation_col,
      palettes = palettes,
      palette = palette
    )
  }

  if (
    !is.null(annotation_colors) &&
    !is.null(annotation_col)
  ) {

    valid_annotations <- intersect(
      names(annotation_colors),
      colnames(annotation_col)
    )

    annotation_colors <- annotation_colors[
      valid_annotations
    ]

    for (v in valid_annotations) {

      levs <- levels(
        annotation_col[[v]]
      )

      annotation_colors[[v]] <-
        annotation_colors[[v]][
          levs[
            levs %in%
              names(annotation_colors[[v]])
          ]
        ]
    }
  }

  heatmap_colors <- grDevices::colorRampPalette(
    palettes$heatmap
  )(100)

  # ---------------------------------------------------------- #
  # Dimensions and text sizes
  # ---------------------------------------------------------- #

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  hm_size <- .deggo_heatmap_size(
    n_genes = nrow(mat_use),
    n_samples = ncol(mat_use)
  )

  width <- width %||% hm_size$width
  height <- height %||% hm_size$height

  fontsize_row_use <- fontsize_row %||%
    if (nrow(mat_use) <= 25L) {
      10
    } else if (nrow(mat_use) <= 50L) {
      8
    } else {
      6
    }

  fontsize_col_use <- fontsize_col %||%
    if (ncol(mat_use) <= 12L) {
      10
    } else if (ncol(mat_use) <= 25L) {
      8
    } else {
      6
    }

  show_rownames <- show_rownames %||%
    nrow(mat_use) <= 80L

  show_colnames <- show_colnames %||%
    ncol(mat_use) <= 60L

  outfile <- file.path(
    output_dir,
    paste0(filename, ".png")
  )

  # ---------------------------------------------------------- #
  # Plot
  # ---------------------------------------------------------- #

  pheatmap::pheatmap(
    mat_use,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    fontsize_row = fontsize_row_use,
    fontsize_col = fontsize_col_use,
    fontsize = 9,
    border_color = NA,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    treeheight_row = if (
      isTRUE(cluster_rows) &&
      isTRUE(show_dendrogram)
    ) {
      50
    } else {
      0
    },
    treeheight_col = if (
      isTRUE(cluster_cols) &&
      isTRUE(show_dendrogram)
    ) {
      50
    } else {
      0
    },
    angle_col = 90,
    main = main,
    color = heatmap_colors,
    breaks = seq(
      -2,
      2,
      length.out = 101
    ),
    filename = outfile,
    width = width,
    height = height
  )

  if (!file.exists(outfile)) {
    return(
      .return_null(
        "Heatmap file was not created by pheatmap."
      )
    )
  }

  invisible(mat_use)
}

# ========================================================= #
# GENE EXPRESSION HEATMAP
# ========================================================= #
#' Plot Expression Heatmap for Selected Genes
#'
#' Generate a publication-ready heatmap for user-defined genes using raw counts
#' and sample metadata. The function matches samples, transforms expression,
#' optionally scales genes by row, adds metadata annotations, and exports a PNG
#' heatmap.
#'
#' By default, genes are clustered but dendrograms are hidden for a clean
#' report-friendly DEGgo visualization. Samples are not clustered by default,
#' preserving metadata-defined order.
#'
#' @param counts Count matrix or count table. If a data frame is provided, one
#'   column must contain gene identifiers.
#' @param metadata Sample metadata data frame.
#' @param genes Character vector of gene IDs or gene symbols to display.
#' @param gene_col Candidate gene ID columns.
#' @param feature_col Candidate gene symbol/name columns.
#' @param sample_col Candidate sample identifier columns in `metadata`.
#' @param assay_transform Expression transformation, either `"log2"` or
#'   `"log2cpm"`.
#' @param annotation_cols Metadata columns shown above the heatmap.
#' @param annotation_colors Optional annotation colors passed to `pheatmap`.
#' @param order_by Optional metadata columns used to order samples.
#' @param output_dir Output directory.
#' @param filename Output filename without extension.
#' @param main Heatmap title.
#' @param color Heatmap color palette.
#' @param breaks Numeric vector of color breaks.
#' @param scale_rows Logical. If `TRUE`, scale expression by gene.
#' @param cluster_rows Logical. If `TRUE`, cluster genes.
#' @param cluster_cols Logical. If `TRUE`, cluster samples.
#' @param fontsize_row Optional row label font size. If `NULL`, chosen
#'   automatically.
#' @param fontsize_col Optional column label font size. If `NULL`, chosen
#'   automatically.
#' @param width Optional plot width in inches. If `NULL`, chosen automatically.
#' @param height Optional plot height in inches. If `NULL`, chosen automatically.
#' @param show_rownames Optional logical. If `NULL`, shown automatically for
#'   heatmaps with 80 genes or fewer.
#' @param show_colnames Optional logical. If `NULL`, shown automatically for
#'   heatmaps with 60 samples or fewer.
#' @param show_dendrogram Logical. If TRUE, display row/column dendrograms
#'   when clustering is enabled.
#'
#' @return Invisibly returns the expression matrix used for plotting.
#'
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
    color = grDevices::colorRampPalette(c("#6497B1", "#F7F7F7", "#740001"))(100),
    breaks = seq(-2, 2, length.out = 101),
    scale_rows = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = NULL,
    fontsize_col = NULL,
    show_dendrogram = FALSE,
    width = NULL,
    height = NULL,
    show_rownames = NULL,
    show_colnames = NULL
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

  if (is.na(gene_col)) {
    stop("No gene ID column found.", call. = FALSE)
  }

  if (is.na(sample_col)) {
    stop("No sample column found.", call. = FALSE)
  }

  metadata[[sample_col]] <- as.character(metadata[[sample_col]])
  rownames(metadata) <- metadata[[sample_col]]

  sample_cols <- metadata[[sample_col]][metadata[[sample_col]] %in% colnames(counts)]

  if (length(sample_cols) < 2) {
    stop("Fewer than 2 matching sample columns found.", call. = FALSE)
  }

  metadata <- metadata[sample_cols, , drop = FALSE]

  mat <- as.matrix(counts[, sample_cols, drop = FALSE])
  rownames(mat) <- as.character(counts[[gene_col]])

  suppressWarnings(storage.mode(mat) <- "numeric")

  if (anyNA(mat)) {
    stop(
      "Expression matrix contains NA after numeric conversion. ",
      "Check that sample columns contain only numeric counts.",
      call. = FALSE
    )
  }

  feature <- if (!is.na(feature_col)) {
    as.character(counts[[feature_col]])
  } else {
    rownames(mat)
  }

  genes_query <- unique(as.character(genes))

  keep <- rownames(mat) %in% genes_query |
    toupper(feature) %in% toupper(genes_query)

  if (!any(keep)) {
    stop("None of the requested genes were found.", call. = FALSE)
  }

  mat_use <- mat[keep, , drop = FALSE]
  feature_use <- feature[keep]

  feature_use[is.na(feature_use) | feature_use == ""] <- rownames(mat_use)[
    is.na(feature_use) | feature_use == ""
  ]

  rownames(mat_use) <- make.unique(feature_use)

  if (assay_transform == "log2") {
    mat_use <- log2(mat_use + 1)
  } else {
    lib <- colSums(mat)
    lib[lib == 0] <- NA_real_
    mat_use <- log2(t(t(mat_use) / lib * 1e6) + 1)
  }

  mat_use <- mat_use[
    apply(mat_use, 1, function(x) all(is.finite(x))),
    apply(mat_use, 2, function(x) all(is.finite(x))),
    drop = FALSE
  ]

  if (nrow(mat_use) < 1 || ncol(mat_use) < 2) {
    stop("Not enough valid genes or samples for heatmap.", call. = FALSE)
  }

  metadata_use <- metadata[colnames(mat_use), , drop = FALSE]

  if (!is.null(order_by)) {
    missing <- setdiff(order_by, colnames(metadata_use))

    if (length(missing)) {
      stop(
        "order_by column(s) not found: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
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
    gene_sd <- apply(mat_use, 1, stats::sd, na.rm = TRUE)
    keep_genes <- is.finite(gene_sd) & gene_sd > 0

    mat_use <- mat_use[keep_genes, , drop = FALSE]

    if (nrow(mat_use) < 1) {
      stop("No valid genes after removing zero-variance genes.", call. = FALSE)
    }

    mat_use <- t(scale(t(mat_use)))

    mat_use <- mat_use[
      apply(mat_use, 1, function(x) all(is.finite(x))),
      ,
      drop = FALSE
    ]
  }

  if (nrow(mat_use) < 1 || ncol(mat_use) < 2) {
    stop("Final heatmap matrix has fewer than 1 gene or 2 samples.", call. = FALSE)
  }

  metadata_use <- metadata_use[colnames(mat_use), , drop = FALSE]

  if (!is.null(annotation_col)) {
    annotation_col <- annotation_col[colnames(mat_use), , drop = FALSE]
  }

  if (is.null(annotation_colors) || is.character(annotation_colors)) {
    annotation_colors <- if (!is.null(annotation_col)) {
      .deggo_annotation_colors(annotation_col)
    } else {
      NULL
    }
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  hm_size <- .deggo_heatmap_size(
    n_genes = nrow(mat_use),
    n_samples = ncol(mat_use)
  )

  if (is.null(width)) width <- hm_size$width
  if (is.null(height)) height <- hm_size$height

  if (is.null(fontsize_row)) {
    fontsize_row_use <- if (nrow(mat_use) <= 20) {
      9
    } else if (nrow(mat_use) <= 40) {
      7
    } else if (nrow(mat_use) <= 80) {
      6
    } else {
      5
    }
  } else {
    fontsize_row_use <- fontsize_row
  }

  if (is.null(fontsize_col)) {
    fontsize_col_use <- if (ncol(mat_use) <= 12) {
      9
    } else if (ncol(mat_use) <= 25) {
      7
    } else if (ncol(mat_use) <= 50) {
      6
    } else {
      5
    }
  } else {
    fontsize_col_use <- fontsize_col
  }

  if (is.null(show_rownames)) {
    show_rownames <- nrow(mat_use) <= 80
  }

  if (is.null(show_colnames)) {
    show_colnames <- ncol(mat_use) <= 60
  }

  outfile <- file.path(output_dir, paste0(filename, ".png"))

  pheatmap::pheatmap(
    mat_use,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    fontsize = 9,
    fontsize_row = fontsize_row_use,
    fontsize_col = fontsize_col_use,
    border_color = NA,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    treeheight_row = if (isTRUE(cluster_rows) && isTRUE(show_dendrogram)) 50 else 0,
    treeheight_col = if (isTRUE(cluster_cols) && isTRUE(show_dendrogram)) 50 else 0,
    angle_col = 90,
    main = main,
    color = color,
    breaks = breaks,
    filename = outfile,
    width = width,
    height = height
  )

  invisible(mat_use)
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
#' @param style Theme style passed to \code{.deggo_theme()}.
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
#'
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
    dpi = 300,
    style = "classic"
) {

  for (pkg in c("ggplot2", "pheatmap")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required.", call. = FALSE)
    }
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  qc_input <- .deggo_qc_prepare(
    counts = counts,
    metadata = metadata,
    gene_col = gene_col,
    feature_col = feature_col,
    sample_col = sample_col
  )

  mat <- qc_input$mat
  metadata <- qc_input$metadata
  counts_df <- qc_input$counts_df
  gene_col <- qc_input$gene_col
  feature_col <- qc_input$feature_col
  sample_col <- qc_input$sample_col

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

  fill_aes <- if (!is.null(color_by)) {
    ggplot2::aes(fill = .data[[color_by]])
  } else {
    ggplot2::aes()
  }

  color_aes <- if (!is.null(color_by)) {
    ggplot2::aes(color = .data[[color_by]], group = sample)
  } else {
    ggplot2::aes(color = sample)
  }

  p_lib <- ggplot2::ggplot(
    qc,
    ggplot2::aes(stats::reorder(sample, library_size), library_size / 1e6)
  ) +
    fill_aes +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    .deggo_theme(style = style, txtsize = 8) +
    ggplot2::labs(x = NULL, y = "Library size (million reads)", title = "Library size")

  p_det <- ggplot2::ggplot(
    qc,
    ggplot2::aes(stats::reorder(sample, detected_genes), detected_genes)
  ) +
    fill_aes +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    .deggo_theme(style = style, txtsize = 8) +
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
    .deggo_theme(style = style, txtsize = 8) +
    ggplot2::labs(x = NULL, y = "log2(count + 1)", title = "Expression distribution")

  p_density <- ggplot2::ggplot(df_long, ggplot2::aes(expression)) +
    color_aes +
    ggplot2::geom_density(linewidth = 0.4, alpha = 0.7) +
    .deggo_theme(style = style, txtsize = 8) +
    ggplot2::labs(x = "log2(count + 1)", y = "Density", title = "Expression density")

  ggplot2::ggsave(file.path(output_dir, "Expression_Boxplot.png"), p_box, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(file.path(output_dir, "Expression_Density.png"), p_density, width = width, height = height, dpi = dpi)

  annotation_col <- metadata[, intersect(annotation_cols, colnames(metadata)), drop = FALSE]

  if (ncol(annotation_col)) {
    annotation_col[] <- lapply(annotation_col, factor)
  } else {
    annotation_col <- NULL
  }

  ann_colors <- .deggo_annotation_colors(annotation_col)

  cor_mat <- stats::cor(logmat, method = "spearman")

  hm_size <- .deggo_heatmap_size(
    n_genes = nrow(mat),
    n_samples = ncol(mat)
  )

  sample_hm_size <- .deggo_heatmap_size(
    n_genes = ncol(mat),
    n_samples = ncol(mat),
    base_width = 5,
    base_height = 5,
    sample_unit = 0.18,
    gene_unit = 0.18,
    min_width = 6,
    max_width = 14,
    min_height = 6,
    max_height = 14
  )

  heatmap_width <- sample_hm_size$width
  heatmap_height <- sample_hm_size$height

  fontsize_sample <- if (ncol(mat) <= 12) {
    8
  } else if (ncol(mat) <= 30) {
    6
  } else {
    5
  }


  pheatmap::pheatmap(
    cor_mat,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    treeheight_row = 0,
    treeheight_col = 0,
    fontsize = 9,
    fontsize_row = fontsize_sample,
    fontsize_col = fontsize_sample,
    border_color = NA,
    show_rownames = ncol(mat) <= 60,
    show_colnames = ncol(mat) <= 60,
    filename = file.path(output_dir, "Sample_Correlation_Heatmap.png"),
    main = "Sample correlation",
    color = grDevices::colorRampPalette(c("#6497B1", "#F7F7F7", "#740001"))(100),
    width = heatmap_width,
    height = heatmap_height
  )

  grDevices::png(
    file.path(output_dir, "Hierarchical_Clustering.png"),
    width = 12,
    height = 6,
    units = "in",
    res = dpi
  )
  plot(stats::hclust(stats::dist(t(logmat))), main = "Hierarchical clustering", xlab = "", sub = "")
  grDevices::dev.off()

  vars <- apply(logmat, 1, stats::var, na.rm = TRUE)
  vars <- vars[is.finite(vars) & vars > 0]

  pca_plots <- NULL
  pca <- NULL

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
        if (length(missing)) {
          tmp <- c(tmp, stats::setNames(grDevices::hcl.colors(length(missing), "Dark 3"), missing))
        }
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
        .deggo_theme(style = style, txtsize = 8) +
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

  top <- character(0)

  if (length(vars) >= 2) {

    top <- names(sort(vars, decreasing = TRUE))[seq_len(min(top_variable_genes, length(vars)))]

    top_mat <- logmat[top, , drop = FALSE]
    top_mat <- t(scale(t(top_mat)))
    top_mat <- top_mat[stats::complete.cases(top_mat), , drop = FALSE]

    if (nrow(top_mat) >= 2) {
      top_hm_size <- .deggo_heatmap_size(
        n_genes = min(nrow(top_mat), 100),
        n_samples = ncol(top_mat),
        base_width = 5.5,
        base_height = 4.5,
        gene_unit = 0.08,
        sample_unit = 0.20,
        min_width = 6,
        max_width = 14,
        min_height = 5,
        max_height = 12
      )

      pheatmap::pheatmap(
        top_mat,
        annotation_col = annotation_col,
        annotation_colors = ann_colors,
        fontsize = 9,
        fontsize_row = 5,
        fontsize_col = fontsize_sample,
        show_rownames = FALSE,
        show_colnames = ncol(top_mat) <= 60,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        treeheight_row = 0,
        treeheight_col = 0,
        border_color = NA,
        color = grDevices::colorRampPalette(c("#6497B1", "#F7F7F7", "#740001"))(100),
        breaks = seq(-2, 2, length.out = 101),
        filename = file.path(output_dir, "TopVariableGenes_Heatmap.png"),
        width = top_hm_size$width,
        height = top_hm_size$height
      )
    }
  }

  marker_heatmap <- NULL

  if (!is.null(markers) && length(markers) && exists("plot_gene_heatmap")) {
    marker_heatmap <- tryCatch(
      plot_gene_heatmap(
        counts = counts_df,
        metadata = metadata,
        genes = markers,
        gene_col = gene_col,
        feature_col = feature_col,
        sample_col = sample_col
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

  utils::write.table(
    qc[qc$recommend_remove, , drop = FALSE],
    file.path(output_dir, "Samples_Recommended_For_Removal.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  p_qc <- ggplot2::ggplot(
    qc,
    ggplot2::aes(stats::reorder(sample, n_flags), n_flags, fill = qc_status)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = .deggo_colors()$qc) +
    .deggo_theme(style = style, txtsize = 8) +
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
#' @param palette Optional color palette used for GO-term visualization.
#' @param x_angle X-axis text angle.
#' @param txtsize Base font size.
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
      Up = "#740001",
      Down = "#6497b1"
    ),
    palette = "default",
    size_range = c(2, 6),
    style = "bw",
    x_angle = 0,
    txtsize = 12
) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.", call. = FALSE)

  Description_wrapped <- FoldEnrichment <- Regulation <- log10FDR <- NULL

  if (!requireNamespace("stringr", quietly = TRUE))
    stop("Package 'stringr' is required.", call. = FALSE)

  if (!requireNamespace("forcats", quietly = TRUE))
    stop("Package 'forcats' is required.", call. = FALSE)

  if (is.null(go_df) || !is.data.frame(go_df) || nrow(go_df) == 0)
    stop("go_df is empty.", call. = FALSE)

  required <- c(
    "Description",
    "FoldEnrichment",
    "p.adjust",
    "Count",
    "Regulation"
  )

  miss <- setdiff(required, colnames(go_df))

  if (length(miss))
    stop(
      "Missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )

  ## ----------------------------- #
  ## Clean
  ## ----------------------------- #

  df <- go_df

  df$Description <- toupper(trimws(df$Description))
  df$Regulation <- factor(df$Regulation,
                          levels = c("Up", "Down"))

  df$FoldEnrichment <- as.numeric(df$FoldEnrichment)
  df$p.adjust <- as.numeric(df$p.adjust)
  df$Count <- as.integer(df$Count)

  keep <-
    !is.na(df$Regulation) &
    is.finite(df$FoldEnrichment) &
    is.finite(df$p.adjust) &
    df$p.adjust > 0

  df <- df[keep, , drop = FALSE]

  if (!nrow(df))
    stop("No valid GO terms after cleaning.", call. = FALSE)

  ## ----------------------------- #
  ## Keyword filtering
  ## ----------------------------- #

  if (!is.null(go_terms_of_interest)) {

    pattern <- paste(
      toupper(go_terms_of_interest),
      collapse = "|"
    )

    keep <- grepl(
      pattern,
      df$Description,
      ignore.case = TRUE
    )

    df <- df[keep, , drop = FALSE]

    if (!nrow(df))
      stop(
        "No GO terms matched the requested keywords.",
        call. = FALSE
      )
  }

  ## ----------------------------- #
  ## Keep best term per regulation
  ## ----------------------------- #

  ord <- order(
    df$Regulation,
    df$Description,
    df$p.adjust,
    -df$FoldEnrichment
  )

  df <- df[ord, ]

  dup <- duplicated(
    paste(df$Regulation, df$Description)
  )

  df <- df[!dup, , drop = FALSE]

  ## ----------------------------- #
  ## Top N
  ## ----------------------------- #

  if (is.null(go_terms_of_interest)) {

    split_df <- split(df, df$Regulation)

    split_df <- lapply(split_df, function(x) {

      x <- x[order(
        x$p.adjust,
        -x$FoldEnrichment
      ), ]

      utils::head(x, top_n)

    })

    df <- do.call(rbind, split_df)
    rownames(df) <- NULL
  }

  ## ----------------------------- #
  ## Plot variables
  ## ----------------------------- #

  df$log10FDR <- -log10(df$p.adjust)

  df$Description_wrapped <-
    stringr::str_wrap(
      df$Description,
      width = wrap_width
    )

  df$Description_wrapped <-
    forcats::fct_reorder(
      df$Description_wrapped,
      df$log10FDR
    )

  regs <- unique(as.character(df$Regulation))

  color_values <- color_values[
    names(color_values) %in% regs
  ]

  ## ----------------------------- #
  ## Plot
  ## ----------------------------- #

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = log10FDR,
      y = Description_wrapped
    )
  ) +

    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = log10FDR,
        yend = Description_wrapped,
        colour = Regulation
      ),
      linewidth = 0.8
    ) +

    ggplot2::geom_point(
      ggplot2::aes(
        size = FoldEnrichment,
        colour = Regulation
      )
    ) +

    ggplot2::scale_colour_manual(
      values = color_values
    ) +

    ggplot2::scale_size(
      range = size_range
    ) +

    ggplot2::labs(
      x = expression(-log[10](FDR)),
      y = NULL,
      colour = "Regulation",
      size = "Fold enrichment",
      title = ifelse(
        is.null(comparison),
        "GO enrichment",
        comparison
      )
    ) +

    .deggo_theme(
      style = style,
      x.ang = x_angle,
      txtsize = txtsize
    )
}


# ========================================================= #
# PLOT ALL GO TERMS
# ========================================================= #
#' Plot GO enrichment terms across all DEGgo comparisons
#'
#' Generate Gene Ontology enrichment plots for all comparisons stored in a
#' DEGgo results object. This function iterates over `results$go_results`,
#' creates one GO plot per comparison using `plot_go_terms()`, and optionally
#' saves each plot as PNG and PDF.
#'
#' @param results A DEGgo results object returned by `run_deggo()`. Must contain
#'   a `go_results` list.
#' @param top_n Integer. Number of top enriched GO terms to display per plot.
#'   Default is `10`.
#' @param txtsize Numeric. Base text size used in the GO plots. Default is `8`.
#' @param style Character. Plot theme style passed to `plot_go_terms()`.
#'   Default is `"bw"`.
#' @param output_dir Optional character. Directory where GO plots are saved.
#'   If `NULL`, plots are returned but not written to disk.
#' @param width Numeric. Plot width in inches for saved files. Default is `8`.
#' @param height Numeric. Plot height in inches for saved files. Default is `6`.
#' @param dpi Integer. Resolution in dots per inch for saved PNG files.
#'   Default is `300`.
#'
#' @return A named list of `ggplot` objects, one per comparison. Comparisons
#'   without enriched GO terms return `NULL`.
#'
#' @examples
#' \dontrun{
#' go_plots <- plot_all_go_terms(
#'   results = results,
#'   top_n = 10,
#'   output_dir = "GO_plots"
#' )
#' }
#'
#' @export
#'
plot_all_go_terms <- function(
    results,
    top_n = 10,
    txtsize = 8,
    style = "bw",
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
      txtsize = txtsize,
      style = style
    )

    if (!is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

      ggplot2::ggsave(
        filename = file.path(output_dir, paste0(comp, "_GO_terms.png")),
        plot = p,
        width = width,
        height = height,
        dpi = dpi,
        bg = "white"
      )

      ggplot2::ggsave(
        filename = file.path(output_dir, paste0(comp, "_GO_terms.pdf")),
        plot = p,
        width = width,
        height = height,
        bg = "white"
      )
    }

    p
  })

  names(plots) <- names(results$go_results)

  plots
}

#' @keywords internal
#' @noRd
#'
.deggo_pca_vars <- function(
    metadata,
    sample_col = "sample",
    analysis_mode = c("single", "pairwise"),
    design_formula = NULL,
    contrast = NULL,
    pairwise_group_cols = NULL,
    pca_vars = NULL
) {
  analysis_mode <- match.arg(analysis_mode)

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!is.null(pca_vars)) {
    vars <- pca_vars
  } else if (analysis_mode == "pairwise") {
    vars <- pairwise_group_cols
  } else {
    vars <- unique(c(
      all.vars(design_formula),
      contrast[1]
    ))
  }

  vars <- unique(vars)
  vars <- vars[!is.na(vars) & vars != ""]
  vars <- setdiff(vars, sample_col)
  vars <- intersect(vars, colnames(metadata))

  vars <- vars[vapply(vars, function(v) {
    n <- length(unique(stats::na.omit(metadata[[v]])))
    n >= 2 && n <= 20
  }, logical(1))]

  vars
}


#' @keywords internal
#' @noRd
#'
.deggo_make_pca_plots <- function(
    dds,
    metadata,
    output_dir,
    analysis_mode,
    sample_col = "sample",
    pca_vars = NULL,
    txtsize = 12
) {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!sample_col %in% colnames(metadata)) {
    metadata[[sample_col]] <- rownames(metadata)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  plots <- list()

  plots$sample <- plot_pca(
    dds = dds,
    metadata = metadata,
    output_dir = output_dir,
    filename = paste0(analysis_mode, "_PCA_by_", sample_col),
    intgroup = sample_col,
    color_by = sample_col,
    title = paste("PCA by", sample_col),
    txtsize = txtsize
  )

  if (!is.null(pca_vars) && length(pca_vars)) {
    for (v in pca_vars) {
      plots[[v]] <- plot_pca(
        dds = dds,
        metadata = metadata,
        output_dir = output_dir,
        filename = paste0(analysis_mode, "_PCA_by_", v),
        intgroup = v,
        color_by = v,
        title = paste("PCA by", v),
        txtsize = txtsize
      )
    }
  }

  if (!is.null(pca_vars) && length(pca_vars) >= 2) {
    combs <- utils::combn(pca_vars, 2, simplify = FALSE)

    for (cc in combs) {
      v1 <- cc[1]
      v2 <- cc[2]
      nm <- paste(v1, v2, sep = "_")

      plots[[nm]] <- plot_pca(
        dds = dds,
        metadata = metadata,
        output_dir = output_dir,
        filename = paste0(analysis_mode, "_PCA_by_", nm),
        intgroup = c(v1, v2),
        color_by = v1,
        shape_by = v2,
        title = paste("PCA by", v1, "+", v2),
        txtsize = txtsize
      )
    }
  }

  plots
}


#' @keywords internal
#' @noRd
#'
.deggo_make_plots <- function(
    de_results,
    counts,
    metadata,
    dirs,
    analysis_mode,
    method,
    padj_cutoff,
    logfc_cutoff,
    top_n_heatmap,
    top_n_labels,
    txtsize,
    log,
    contrast = NULL,
    order_by = NULL,
    pairwise_contrast_col = NULL,
    pairwise_contrasts = NULL,
    sample_col = "sample",
    pca_vars = NULL,
    plot_order_vars = NULL,
    design_formula = NULL,
    pairwise_group_cols = NULL,
    heatmap_annotation_cols = "auto",
    palette = "default"
) {

  metadata <- as.data.frame(
    metadata,
    stringsAsFactors = FALSE
  )

  if (!sample_col %in% colnames(metadata)) {
    metadata[[sample_col]] <- rownames(metadata)
  }

  metadata[[sample_col]] <- as.character(
    metadata[[sample_col]]
  )

  rownames(metadata) <- metadata[[sample_col]]

  # ---------------------------------------------------------- #
  # PCA and heatmaps
  # ---------------------------------------------------------- #

  if (
    identical(method, "DESeq2") &&
    !is.null(de_results$dds)
  ) {

    log(
      "[10/11] Generating PCA",
      type = "step"
    )

    de_results$pca <- .deggo_make_pca_plots(
      dds = de_results$dds,
      metadata = metadata,
      output_dir = dirs$pca,
      analysis_mode = analysis_mode,
      sample_col = sample_col,
      pca_vars = pca_vars,
      txtsize = txtsize
    )

    if (
      is.null(de_results$pca) ||
      is.null(de_results$pca$sample) ||
      is.null(de_results$pca$sample$vsd)
    ) {
      log(
        "Heatmaps skipped: VST object was not produced by PCA.",
        type = "warn"
      )
    } else {

      vsd <- de_results$pca$sample$vsd

      log(
        "[10/11] Generating heatmaps",
        type = "step"
      )

      sig_tables <- de_results$sig_deg_clean

      if (
        is.null(sig_tables) ||
        !is.list(sig_tables) ||
        !length(sig_tables)
      ) {
        log(
          "Heatmaps skipped: no cleaned significant DEG tables available.",
          type = "warn"
        )
      } else {

        for (nm in names(sig_tables)) {

          sig_df <- sig_tables[[nm]]

          if (
            is.null(sig_df) ||
            !is.data.frame(sig_df) ||
            nrow(sig_df) <= 1L
          ) {
            log(
              paste0(
                "Heatmap skipped for ",
                nm,
                ": fewer than 2 genes."
              ),
              type = "warn"
            )
            next
          }

          heatmap_contrast <- NULL
          heatmap_sample_subset <- NULL
          heatmap_filter <- NULL
          heatmap_order_by <- order_by
          heatmap_order_levels <- NULL
          current_contrast <- NULL

          # -------------------------------------------------- #
          # Single-comparison heatmap settings
          # -------------------------------------------------- #

          if (identical(analysis_mode, "single")) {

            current_contrast <- contrast
            heatmap_contrast <- contrast

            if (
              !is.null(plot_order_vars) &&
              length(plot_order_vars)
            ) {
              heatmap_order_by <- intersect(
                plot_order_vars,
                colnames(metadata)
              )
            }

            if (
              !is.null(contrast) &&
              length(contrast) >= 3L
            ) {

              contrast_var <- contrast[1]
              contrast_levels <- contrast[2:3]

              if (contrast_var %in% colnames(metadata)) {
                heatmap_order_levels <- stats::setNames(
                  list(rev(contrast_levels)),
                  contrast_var
                )
              }
            }
          }

          # -------------------------------------------------- #
          # Pairwise heatmap settings
          # -------------------------------------------------- #

          if (
            identical(analysis_mode, "pairwise") &&
            !is.null(pairwise_contrasts) &&
            nm %in% names(pairwise_contrasts)
          ) {

            current_contrast <- pairwise_contrasts[[nm]]

            if (
              is.null(current_contrast) ||
              length(current_contrast) < 3L
            ) {
              log(
                paste0(
                  "Heatmap skipped for ",
                  nm,
                  ": invalid pairwise contrast."
                ),
                type = "warn"
              )
              next
            }

            contrast_var <- current_contrast[1]
            contrast_levels <- current_contrast[2:3]

            if (!contrast_var %in% colnames(metadata)) {
              log(
                paste0(
                  "Heatmap skipped for ",
                  nm,
                  ": contrast column '",
                  contrast_var,
                  "' was not found in metadata."
                ),
                type = "warn"
              )
              next
            }

            heatmap_sample_subset <- rownames(metadata)[
              !is.na(metadata[[contrast_var]]) &
                metadata[[contrast_var]] %in% contrast_levels
            ]

            if (length(heatmap_sample_subset) < 2L) {
              log(
                paste0(
                  "Heatmap skipped for ",
                  nm,
                  ": fewer than 2 samples matched ",
                  paste(contrast_levels, collapse = " vs "),
                  "."
                ),
                type = "warn"
              )
              next
            }

            heatmap_order_by <- contrast_var

            heatmap_order_levels <- stats::setNames(
              list(rev(contrast_levels)),
              contrast_var
            )
          }

          # -------------------------------------------------- #
          # Resolve annotation columns
          # -------------------------------------------------- #

          annotation_cols_use <-
            .deggo_heatmap_annotation_cols(
              metadata = metadata,
              annotation_cols = heatmap_annotation_cols,
              analysis_mode = analysis_mode,
              design_formula = design_formula,
              contrast = contrast,
              pairwise_group_cols = pairwise_group_cols,
              current_contrast = current_contrast,
              sample_col = sample_col
            )

          if (
            !is.null(heatmap_sample_subset) &&
            length(heatmap_sample_subset)
          ) {

            metadata_annotation_use <- metadata[
              intersect(
                heatmap_sample_subset,
                rownames(metadata)
              ),
              ,
              drop = FALSE
            ]

          } else {

            metadata_annotation_use <- metadata
          }

          # Remove annotation variables that become completely empty
          # after contrast-specific subsetting.
          if (length(annotation_cols_use)) {
            annotation_cols_use <- annotation_cols_use[
              vapply(
                annotation_cols_use,
                function(v) {
                  values <- metadata_annotation_use[[v]]
                  values <- values[
                    !is.na(values) &
                      nzchar(as.character(values))
                  ]

                  length(unique(values)) >= 1L
                },
                logical(1)
              )
            ]
          }

          tryCatch(
            plot_heatmap(
              vsd = vsd,
              res_df = sig_df,
              metadata = metadata,
              contrast = heatmap_contrast,
              sample_subset = heatmap_sample_subset,
              metadata_filter = heatmap_filter,
              annotation_cols = annotation_cols_use,
              order_by = heatmap_order_by,
              order_levels = heatmap_order_levels,
              output_dir = dirs$heatmaps,
              filename = paste0(nm, "_Heatmap"),
              main = nm,
              top_n_heatmap = top_n_heatmap,
              padj_cutoff = padj_cutoff,
              logfc_cutoff = logfc_cutoff,
              cluster_rows = TRUE,
              cluster_cols = FALSE,
              show_dendrogram = FALSE,
              palette = palette
            ),
            error = function(e) {
              log(
                paste0(
                  "Heatmap failed for ",
                  nm,
                  ": ",
                  conditionMessage(e)
                ),
                type = "warn"
              )

              NULL
            }
          )
        }
      }
    }
  }

  # ---------------------------------------------------------- #
  # Volcano plots
  # ---------------------------------------------------------- #

  log(
    "[10/11] Generating volcano plots",
    type = "step"
  )

  result_tables <- de_results$results

  if (
    is.null(result_tables) ||
    !is.list(result_tables) ||
    !length(result_tables)
  ) {
    log(
      "Volcano plots skipped: no differential-expression results available.",
      type = "warn"
    )

    return(de_results)
  }

  for (nm in names(result_tables)) {

    res_df <- result_tables[[nm]]

    if (
      is.null(res_df) ||
      !is.data.frame(res_df) ||
      nrow(res_df) == 0L
    ) {
      log(
        paste0(
          "Volcano plot skipped for ",
          nm,
          ": empty result table."
        ),
        type = "warn"
      )

      next
    }

    tryCatch(
      plot_volcano(
        res_df = res_df,
        title = nm,
        output_dir = dirs$volcano,
        filename = paste0(
          nm,
          "_Volcano"
        ),
        padj_cutoff = padj_cutoff,
        logfc_cutoff = logfc_cutoff,
        top_n_labels = top_n_labels,
        txtsize = txtsize
      ),
      error = function(e) {
        log(
          paste0(
            "Volcano plot failed for ",
            nm,
            ": ",
            conditionMessage(e)
          ),
          type = "warn"
        )

        NULL
      }
    )
  }

  de_results
}




#' @keywords internal
#' @noRd
#'
.deggo_make_go <- function(
    de_results,
    dirs,
    ontology,
    orgdb,
    txtsize,
    log
) {

  log("[11/11] Running GO enrichment", type = "step")

  go_results <- list()

  go_plot_dir <- file.path(dirs$go, "GO_plots")
  dir.create(dirs$go, recursive = TRUE, showWarnings = FALSE)
  dir.create(go_plot_dir, recursive = TRUE, showWarnings = FALSE)

  for (nm in names(de_results$sig_deg_clean)) {

    sig_df <- de_results$sig_deg_clean[[nm]]

    if (is.null(sig_df) || nrow(sig_df) == 0) {
      go_results[[nm]] <- NULL
      next
    }

    go_obj <- tryCatch(
      run_go_enrichment(
        sig_deg = sig_df,
        comparison = nm,
        ontology = ontology,
        orgdb = orgdb,
        output_dir = dirs$go,
        entrez_col = "entrez_id"
      ),
      error = function(e) {
        log(paste("GO skipped for", nm, ":", conditionMessage(e)), type = "warn")
        NULL
      }
    )

    go_results[[nm]] <- go_obj

    if (is.null(go_obj)) {
      next
    }

    go_df <- if (is.data.frame(go_obj)) {
      go_obj
    } else if (!is.null(go_obj$go_results) && is.data.frame(go_obj$go_results)) {
      go_obj$go_results
    } else {
      NULL
    }

    if (!is.null(go_df) && nrow(go_df) > 0) {
      p_go <- tryCatch(
        plot_go_terms(
          go_df = go_df,
          comparison = nm,
          top_n = 10,
          txtsize = txtsize,
          style = "bw"
        ),
        error = function(e) {
          log(paste("GO plot skipped for", nm, ":", conditionMessage(e)), type = "warn")
          NULL
        }
      )

      if (!is.null(p_go)) {
        ggplot2::ggsave(
          filename = file.path(go_plot_dir, paste0(nm, "_GO_terms.png")),
          plot = p_go,
          width = 8,
          height = 6,
          dpi = 300,
          bg = "white"
        )

        ggplot2::ggsave(
          filename = file.path(go_plot_dir, paste0(nm, "_GO_terms.pdf")),
          plot = p_go,
          width = 8,
          height = 6,
          bg = "white"
        )
      }
    }
  }

  de_results$go_results <- go_results
  de_results
}
