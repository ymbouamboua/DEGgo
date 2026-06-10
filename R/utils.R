#' Format elapsed runtime into human-readable text
#'
#' Converts a duration in seconds into a compact human-readable string
#' suitable for console logging.
#'
#' Formatting rules:
#' \itemize{
#'   \item seconds < 60 → \code{"7.6s"}
#'   \item minutes < 60 → \code{"1m 0.6s"}
#'   \item hours ≥ 1 → \code{"1h 3m 0.4s"}
#' }
#'
#' @param seconds Numeric duration in seconds.
#'
#' @return A character string containing formatted elapsed time.
#'
#' @keywords internal
#' @noRd
.format_duration <- function(seconds) {
  
  if (is.null(seconds) || is.na(seconds)) {
    return(NA_character_)
  }
  
  seconds <- as.numeric(seconds)
  
  hrs <- floor(seconds / 3600)
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



#' Styled GnRHcell logger
#'
#' CRAN-safe console logger using ASCII labels and optional ANSI colors.
#'
#' @param verbose Logical; print messages.
#' @param color Logical; use ANSI colors.
#'
#' @return A logging function.
#' @keywords internal
#' @noRd
.msg <- function(
    verbose = TRUE,
    color = interactive()
) {
  
  t0 <- Sys.time()
  
  col_reset  <- if (color) "\033[0m" else ""
  col_bold   <- if (color) "\033[1m" else ""
  
  # devtools/cli-like palette
  col_blue   <- if (color) "\033[34m" else ""
  col_green  <- if (color) "\033[32m" else ""
  col_yellow <- if (color) "\033[33m" else ""
  col_red    <- if (color) "\033[31m" else ""
  col_gray   <- if (color) "\033[90m" else ""
  col_cyan   <- if (color) "\033[36m" else ""
  col_purple <- if (color) "\033[35m" else ""
  
  function(
    ...,
    type = c("info", "step", "done", "warn", "error", "header"),
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
          difftime(Sys.time(), t0, units = "secs")
        )
      }
      
      txt <- sprintf(
        "%s Duration: %s",
        txt,
        .format_duration(duration)
      )
    }
    
    prefix <- switch(
      type,
      info   = "[INFO]",
      step   = "[STEP]",
      done   = "[DONE]",
      warn   = "[WARN]",
      error  = "[ERROR]",
      header = "[DEGgo]"
    )
    
    prefix_col <- switch(
      type,
      info   = col_gray,
      step   = col_blue,
      done   = col_green,
      warn   = col_yellow,
      error  = col_red,
      header = paste0(col_purple, col_bold)
    )
    
    txt_col <- switch(
      type,
      info   = col_gray,
      step   = col_cyan,
      done   = col_green,
      warn   = col_yellow,
      error  = col_red,
      header = paste0(col_purple, col_bold)
    )
    
    cat(
      prefix_col,
      prefix,
      col_reset,
      " ",
      txt_col,
      txt,
      col_reset,
      "\n",
      sep = ""
    )
    
    invisible(NULL)
  }
}