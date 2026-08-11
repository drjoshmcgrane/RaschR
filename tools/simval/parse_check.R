# Parse-check every script in the simulation battery. Run from the package
# root:  Rscript tools/simval/parse_check.R
# Exits non-zero if any file fails to parse, printing the parse error.
files <- list.files("tools/simval", pattern = "[.]R$", recursive = TRUE,
                    full.names = TRUE)
bad <- 0L
for (f in files) {
  e <- tryCatch({ parse(f); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(e)) {
    bad <- bad + 1L
    cat(sprintf("PARSE FAIL %s\n  %s\n", f, sub("\n.*", "", e)))
  }
}
cat(sprintf("%d files checked, %d parse failures\n", length(files), bad))
if (bad > 0L) quit(status = 1L)
