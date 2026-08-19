# rasch :: shared display formatting
# ===========================================================================
# One formatting vocabulary for every surface that prints results (console
# methods, the HTML report, and the app): probabilities below the display
# resolution read "< 0.001" rather than "0.000"; numeric columns carry
# fixed decimals; integer-valued columns print as integers; logical flags
# print as "*" or blank.
# ===========================================================================

.fmt_p <- function(p, digits = 3) {
  lim <- 10^-digits
  out <- ifelse(is.finite(p) & p < lim, paste0("< ", format(lim, scientific = FALSE)),
                ifelse(is.finite(p), sprintf(paste0("%.", digits, "f"), p), ""))
  out
}

# obs_p and est_p are observed and expected category PROPORTIONS, not
# probabilities from a test: a category nobody chose has an observed
# proportion of exactly zero, which "< 0.001" would misreport.
.is_pcol <- function(nm)
  grepl("^p$|^p_|_p$|^prob$|^p\\.", nm) & !nm %in% c("obs_p", "est_p")

# Format a data frame for display: fixed decimals, "< 0.001" probabilities,
# clean integers, starred logicals. Returns a character data frame.
.fmt_df <- function(d, digits = 3) {
  out <- d
  for (j in seq_along(d)) {
    v <- d[[j]]; nm <- names(d)[j]
    if (is.logical(v)) out[[j]] <- ifelse(is.na(v), "", ifelse(v, "*", ""))
    else if (is.numeric(v)) {
      if (.is_pcol(nm)) out[[j]] <- .fmt_p(v, digits)
      else if (all(is.na(v) | v == round(v))) out[[j]] <-
          ifelse(is.na(v), "", format(v, scientific = FALSE, trim = TRUE))
      else out[[j]] <- ifelse(is.na(v), "",
                              sprintf(paste0("%.", digits, "f"), v))
    }
  }
  out
}

# A results table prints through the vocabulary above rather than through
# base R's, whose default switches to scientific notation for the small
# probabilities a large sample produces: a p of 4e-83 reads "< 0.001" here
# and "4.00e-83" otherwise. Only plain data frames are tagged, so tables
# carrying a class of their own keep their own printing.
.tag_tables <- function(x) {
  if (is.data.frame(x)) {
    if (identical(class(x), "data.frame"))
      class(x) <- c("rasch_table", "data.frame")
    return(x)
  }
  if (is.list(x) && length(x)) x[] <- lapply(x, .tag_tables)
  x
}

#' @export
print.rasch_table <- function(x, ..., digits = 3, n = getOption("rasch.print_rows", 50L)) {
  d <- as.data.frame(x)
  if (!nrow(d)) {
    cat("<empty table>\n")
    return(invisible(x))
  }
  more <- if (is.finite(n) && nrow(d) > n) nrow(d) - n else 0L
  shown <- if (more) utils::head(d, n) else d
  # extra arguments are absorbed rather than forwarded: callers pass
  # row.names = FALSE, which the inner call already sets
  print(.fmt_df(shown, digits), row.names = FALSE)
  if (more)
    cat(sprintf("... %d more row%s; as.data.frame() for the unrounded values\n",
                more, if (more == 1L) "" else "s"))
  invisible(x)
}

# write.csv encodes a small double in scientific notation: a p of 2.3e-17
# lands in the file as "2.3e-17". Raising scipen for the call keeps the file
# in the plain decimals the screen shows, at full precision -- rounding here
# would lose information a saved table should keep.
.write_csv_plain <- function(d, path) {
  op <- options(scipen = 999)
  on.exit(options(op), add = TRUE)
  utils::write.csv(d, path, row.names = FALSE)
}
