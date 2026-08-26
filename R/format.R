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

# Crossed factor cells must be keyed by the factor codes, not by pasting the
# labels and refactoring the result. Pasted labels can collide (for example,
# A = "a:b", B = "c" and A = "a", B = "b:c"). The ordinary labels are
# retained where unique; only genuinely ambiguous labels receive a cell tag.
.factor_keys <- function(x) {
  x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  if (!ncol(x)) stop("at least one factor is needed")
  parts <- lapply(x, function(v) {
    z <- as.character(v)
    ifelse(is.na(z), "N;", paste0("S", nchar(z, type = "bytes"), ":", z, ";"))
  })
  do.call(paste0, parts)
}

.factor_cells <- function(x, sep = ":") {
  x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  if (!ncol(x)) stop("at least one factor is needed")
  fs <- lapply(x, function(v) {
    if (is.factor(v)) droplevels(v) else factor(v)
  })
  code <- rep(1, nrow(x)); mult <- 1
  for (f in fs) {
    code <- code + (as.integer(f) - 1) * mult
    mult <- mult * max(nlevels(f), 1)
  }
  present <- sort(unique(code[!is.na(code)]))
  first <- match(present, code)
  labels <- vapply(first, function(i)
    paste(vapply(fs, function(f) as.character(f[i]), ""), collapse = sep), "")
  clash <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  if (any(clash))
    labels[clash] <- paste0(labels[clash], " [cell ", present[clash], "]")
  # A literal user label can itself equal a generated disambiguation label.
  # In that exceptional case the stable mixed-radix code is unambiguous.
  if (anyDuplicated(labels)) labels <- paste0("cell ", present, ": ", labels)
  factor(code, levels = present, labels = labels)
}

# A column-selecting argument must resolve to exactly one existing column
# before it is dereferenced: an empty or multiple name otherwise fails with
# a base subscript error that names neither the argument nor the problem.
.check_reshape_column <- function(data, x, name) {
  if (length(x) != 1L || is.na(x) || !(is.character(x) || is.numeric(x)))
    stop("`", name, "` must name exactly one column", call. = FALSE)
  x <- as.character(x)
  if (!x %in% names(data)) stop("column not found: ", x, call. = FALSE)
  invisible(x)
}

.check_column_names <- function(x) {
  if (!is.data.frame(x)) return(invisible(NULL))
  nm <- names(x)
  if (anyNA(nm) || any(!nzchar(nm)))
    stop("data column names must be non-missing and non-empty")
  if (anyDuplicated(nm))
    stop("data column names must be unique: ",
         paste(unique(nm[duplicated(nm)]), collapse = ", "))
  invisible(NULL)
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
