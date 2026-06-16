# =========================================================
# FORMAT RUNTIME
# =========================================================

#' Format runtime
#' @keywords internal
#' @noRd
.format_duration <- function(seconds) {

  if (is.null(seconds) || is.na(seconds)) {
    return(NA_character_)
  }

  seconds <- as.numeric(seconds)

  hrs  <- floor(seconds / 3600)
  mins <- floor((seconds %% 3600) / 60)
  secs <- seconds %% 60

  if (hrs > 0) {
    sprintf("%dh %dm %.1fs", hrs, mins, secs)
  } else if (mins > 0) {
    sprintf("%dm %.1fs", mins, secs)
  } else {
    sprintf("%.1fs", secs)
  }
}


# =========================================================
# LOGGER
# =========================================================

#' Create package logger
#'
#' Lightweight logger used throughout GnRHcell.
#'
#' @param verbose Logical.
#' @param prefix Package prefix.
#'
#' @return Logging function.
#'
#' @keywords internal
#' @noRd
.msg <- function(
    verbose = TRUE,
    prefix = "INFO"
) {

  t0 <- Sys.time()

  function(
    ...,
    type = c(
      "info",
      "step",
      "done",
      "warn",
      "error",
      "header"
    ),
    duration = NULL
  ) {

    if (!isTRUE(verbose)) {
      return(invisible(NULL))
    }

    type <- match.arg(type)

    txt <- paste(..., collapse = " ")

    if (type == "done") {

      if (is.null(duration)) {
        duration <- as.numeric(
          difftime(
            Sys.time(),
            t0,
            units = "secs"
          )
        )
      }

      txt <- paste0(
        txt,
        " Duration: ",
        .format_duration(duration)
      )
    }

    tag <- switch(
      type,
      info   = "[INFO]",
      step   = "[STEP]",
      done   = "[DONE]",
      warn   = "[WARN]",
      error  = "[ERROR]",
      header = paste0("[", prefix, "]")
    )

    cat(
      tag,
      " ",
      txt,
      "\n",
      sep = ""
    )

    invisible(NULL)
  }
}


# =========================================================
# SAVE SESSION INFO
# =========================================================

#' Save session info
#'
#' @param output_dir outdir directory
#'
#' @keywords internal
#' @noRd
.safe_write_session_info <- function(output_dir) {
  x <- tryCatch(
    utils::capture.output(utils::sessionInfo()),
    error = function(e) {
      c(
        "sessionInfo() failed.",
        paste("Reason:", conditionMessage(e)),
        "",
        "Fallback:",
        utils::capture.output(Sys.info()),
        utils::capture.output(R.version)
      )
    }
  )

  writeLines(
    x,
    file.path(output_dir, "sessionInfo.txt")
  )
}


# =========================================================
# .flag_low_mad
# =========================================================

#' Flag low outliers using MAD
#'
#' @param x Numeric vector.
#' @param nmad Number of MADs below the median.
#'
#' @return Logical vector.
#' @keywords internal
#' @noRd
.flag_low_mad <- function(x, nmad = 3) {
  med <- stats::median(x, na.rm = TRUE)
  madv <- stats::mad(x, na.rm = TRUE)
  x < med - nmad * madv
}

# =========================================================
# .flag_high_mad
# =========================================================

#' Flag high outliers using MAD
#'
#' @param x Numeric vector.
#' @param nmad Number of MADs above the median.
#'
#' @return Logical vector.
#' @keywords internal
#' @noRd
.flag_high_mad <- function(x, nmad = 3) {
  med <- stats::median(x, na.rm = TRUE)
  madv <- stats::mad(x, na.rm = TRUE)
  x > med + nmad * madv
}

