test_that("every Shiny result card has one canonical explainer", {
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  help <- testthat::test_path("..", "..", "inst", "shiny", "help.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  if (!file.exists(help)) help <- system.file("shiny", "help.R", package = "rasch")

  h <- new.env(parent = baseenv())
  sys.source(help, envir = h)
  # duplicate names silently shadow later entries (APP_HELP["x"] returns
  # only the first match), so the registry must be unique by construction
  expect_false(anyDuplicated(names(h$APP_HELP)) > 0,
               info = paste("duplicate APP_HELP names:",
                            paste(unique(names(h$APP_HELP)[
                              duplicated(names(h$APP_HELP))]), collapse = ", ")))

  src <- paste(readLines(app, warn = FALSE), collapse = "\n")
  card_ids <- unique(unlist(lapply(c("tableCard", "plotCard", "statCard"),
    function(fun) {
      hit <- regmatches(src, gregexpr(
        sprintf("%s\\(\\\"[^\\\"]+", fun), src, perl = TRUE))[[1]]
      if (!length(hit) || identical(hit, character(0))) return(character(0))
      sub(sprintf("^%s\\(\\\"", fun), "", hit)
    })))

  expect_true(all(card_ids %in% names(h$APP_HELP)),
              info = paste("missing:",
                           paste(setdiff(card_ids, names(h$APP_HELP)),
                                 collapse = ", ")))

  hand_built_tables <- c(
    "chisq_int_tbl", "chisq_cat_tbl", "ctt_tbl", "rescore_tbl",
    "dif_posthoc_tbl", "dif_size_tbl", "resolve_tbl", "contr_tbl",
    "eq_tbl", "dm_tbl", "dep_tbl", "spread_tbl", "cmp_tbl",
    "sim_recovery_tbl", "sim_preview", "preview"
  )
  expect_true(all(hand_built_tables %in% names(h$APP_HELP)))
  expect_true(all(vapply(hand_built_tables, function(id)
    grepl(sprintf('app_help\\("%s"\\)', id), src), logical(1))),
    info = "a hand-built table is missing its information control")
  expect_true(all(nzchar(h$APP_HELP)))
})

test_that("Shiny explainers remain succinct", {
  help <- testthat::test_path("..", "..", "inst", "shiny", "help.R")
  if (!file.exists(help)) help <- system.file("shiny", "help.R", package = "rasch")
  h <- new.env(parent = baseenv())
  sys.source(help, envir = h)
  words <- lengths(strsplit(trimws(h$APP_HELP), "[[:space:]]+"))
  expect_true(max(words) <= 60,
              info = paste(names(which.max(words)), "has", max(words), "words"))
})

test_that("every analytical result has an R-code disclosure", {
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  tree <- parse(app)

  calls_named <- function(x, names) {
    out <- list()
    walk <- function(z) {
      if (is.call(z)) {
        fun <- tryCatch(as.character(z[[1]])[1], error = function(e) "")
        if (fun %in% names) out[[length(out) + 1L]] <<- z
        if (length(z) > 1L)
          for (i in 2:length(z)) if (length(z[[i]])) walk(z[[i]])
      } else if (is.expression(z) || is.pairlist(z)) {
        for (i in seq_along(z)) if (length(z[[i]])) walk(z[[i]])
      }
    }
    walk(x)
    out
  }
  first_string <- function(z)
    if (length(z) >= 2L && is.character(z[[2]]) && length(z[[2]]) == 1L)
      z[[2]] else NA_character_

  ui_calls <- calls_named(tree,
    c("plotCard", "tableCard", "statCard", "rcode_details"))
  server_calls <- calls_named(tree,
    c("register_plot", "register_table", "register_stat_box", "register_code"))
  ui_ids <- unique(na.omit(vapply(ui_calls, first_string, character(1))))
  server_ids <- unique(na.omit(vapply(server_calls, first_string, character(1))))
  expect_setequal(ui_ids, server_ids)

  managed <- calls_named(tree,
    c("register_plot", "register_table", "register_stat_box"))
  has_code <- vapply(managed, function(z) "code" %in% names(as.list(z)),
                     logical(1))
  expect_true(all(has_code),
    info = paste("missing code argument:",
      paste(vapply(managed[!has_code], first_string, character(1)),
            collapse = ", ")))
})

test_that("bundled app examples are exactly reconstructible", {
  examples <- testthat::test_path("..", "..", "inst", "shiny", "examples.R")
  if (!file.exists(examples))
    examples <- system.file("shiny", "examples.R", package = "rasch")
  e <- new.env(parent = asNamespace("rasch"))
  sys.source(examples, envir = e)
  expect_identical(.app_example_data("pcm"), e$.demo_data())
  expect_identical(.app_example_data("dich"), e$.demo_dich())
  expect_identical(.app_example_data("rsm"), e$.demo_rsm())
  expect_identical(.app_example_data("mfrm"), e$.demo_mfrm())
  expect_identical(.app_example_data("efrm"), e$.demo_efrm())
  expect_identical(.app_example_data("btl"), e$.demo_btl())
})
