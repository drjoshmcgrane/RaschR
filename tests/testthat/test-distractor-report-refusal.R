.report_keyed_fixture <- function(repeated = TRUE) {
  set.seed(9712)
  n <- 80L
  raw <- matrix(sample(c("A", "B"), n * 6, TRUE), n, 6,
                 dimnames = list(NULL, paste0("I", 1:6)))
  rows <- if (repeated) rep(seq_len(n), each = 2L) else seq_len(n)
  rasch(raw[rows, ], id = rows,
        key = setNames(rep("A", 6), colnames(raw)))
}

test_that("distractor reports distinguish design refusals from faults", {
  repeated <- .report_keyed_fixture()
  z <- .distractors_or_note(repeated)
  expect_true(repeated$est$converged)
  expect_null(z$table)
  expect_match(z$note, "one response row per person")
  independent <- .report_keyed_fixture(FALSE)
  ordinary <- .distractors_or_note(independent)
  expect_null(ordinary$note)
  expect_identical(ordinary$table, distractor_analysis(independent))
  expect_null(.distractors_or_note(rasch(independent$X)))
  testthat::with_mocked_bindings({
    expect_error(.distractors_or_note(repeated), "unexpected calculation failure")
    out <- tempfile("unexpected-distractor-")
    expect_error(save_outputs(repeated, out, item_plots = FALSE),
                 "unexpected calculation failure")
    expect_false(dir.exists(out))
  }, distractor_analysis = function(...) stop("unexpected calculation failure"),
    .package = "rasch")
})

test_that("keyed repeated-person exports complete and retain the refusal", {
  fit <- .report_keyed_fixture()
  out <- tempfile("repeated-keyed-output-")
  html <- tempfile(fileext = ".html")
  on.exit(unlink(c(out, html), recursive = TRUE), add = TRUE)
  # Parallel residual inference is separately unavailable for repeated IDs;
  # its existing export warning is not the distractor result under test.
  expect_no_error(suppressWarnings(save_outputs(
    fit, out, formats = "png", item_plots = TRUE, dpi = 40)))
  note <- utils::read.csv(file.path(out, "tables", "distractor_analysis.csv"))
  expect_identical(names(note), "note")
  expect_match(note$note, "one response row per person")
  summary <- paste(readLines(file.path(out, "summary.txt")), collapse = "\n")
  expect_match(summary, "Distractor analysis unavailable:")
  expect_match(summary, "one response row per person")
  plots <- list.files(file.path(out, "plots"), recursive = TRUE)
  expect_true(any(grepl("_icc\\.png$", plots)))
  expect_false(any(grepl("_options\\.png$", plots)))
  expect_no_error(suppressWarnings(report_html(fit, html, dpi = 40)))
  document <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(document, "<h2>Distractor analysis</h2>", fixed = TRUE)
  expect_match(document, "Not available: distractor analysis needs one response row per person",
               fixed = TRUE)
  expect_match(document, "<h2>Person estimates</h2>", fixed = TRUE)
})

test_that("rendered reports retain keyed repeated-person refusal notes", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  fit <- .report_keyed_fixture()
  html <- tempfile(fileext = ".html")
  on.exit(unlink(html), add = TRUE)
  expect_no_error(report_document(fit, html))
  document <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(document, "Distractor analysis", fixed = TRUE)
  expect_match(document, "one response row per person", fixed = TRUE)
  expect_match(document, "Person estimates", fixed = TRUE)
})
