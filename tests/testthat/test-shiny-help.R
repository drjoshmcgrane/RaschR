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

test_that("the app retains the agreed model labels and frame safeguards", {
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  src <- paste(readLines(app, warn = FALSE), collapse = "\n")

  expect_match(src, '"Multiple Ratings \\(MFRM\\)" = "mfrm"')
  expect_match(src, '"Extended Frames \\(EFRM\\)" = "efrm"')
  expect_false(grepl('"Many-facet \\(MFRM\\)" = "mfrm"', src))
  expect_false(grepl('"Rated / many-facet \\(MFRM\\)" = "mfrm"', src))

  # Adding BTL frames must not silently discard specifications unsupported by
  # btl_efrm(): each case is stopped before the frame fit is called.
  for (text in c("uses external anchors", "uses aggregated count weights",
                 "has ordered response categories", "assigns half a win to ties",
                 "includes within-judge order effects"))
    expect_match(src, text, fixed = TRUE)

  expect_match(src, 'if (inherits(f, "rasch_efrm")) f$item_arbitrary',
               fixed = TRUE)
  expect_match(src, "before <- app_item_estimates(analysis())", fixed = TRUE)
  expect_match(src, "after <- app_item_estimates(fit())", fixed = TRUE)
  expect_match(src, "rho = phi / alpha", fixed = TRUE)
  expect_false(grepl("rho = phi x alpha", src, fixed = TRUE))
  expect_match(src, 'style_lo_red(num_dt(d), d, "p_adj", 0.05)',
               fixed = TRUE)
  expect_match(src, '!all(c("object", "location") %in% names(a))',
               fixed = TRUE)
  expect_match(src, 'if (!"se" %in% names(a)) a$se <- NA_real_',
               fixed = TRUE)
  expect_match(src, 'identical(r$shift_method, "unweighted")',
               fixed = TRUE)
  expect_match(src, 'p_bold = c("p_adj", "p_anova_adj", "chisq_p_boot_adj",',
               fixed = TRUE)
  expect_false(grepl('p_bold = c("p", "p_adj")', src, fixed = TRUE))
  expect_match(src, ': judge-group DIF needs judge-constant factors',
               fixed = TRUE)
  expect_match(src, 'judge_factor <- function(x, label)', fixed = TRUE)
  expect_match(src, '"p_adj_kappa"', fixed = TRUE)
  expect_match(src, "maxit = eo$maxit, tol = eo$tol", fixed = TRUE)
  expect_match(src, 'input$spread_alpha %||% 0.05', fixed = TRUE)
  expect_match(src, 'run_adjust <- "holm"', fixed = TRUE)
  expect_false(grepl("dif_padj|spread_padj|inv_screen", src))
  expect_match(src, 'spread_test(fit, alpha = %s, p_adjust = %s)',
               fixed = TRUE)
  expect_match(src, 'callr::r_bg', fixed = TRUE)
  expect_match(src, 'Cancel EFRM estimation', fixed = TRUE)
  expect_match(src, 'input$ef_workers', fixed = TRUE)
  expect_match(src, 'selected = max(.efrm_worker_values)', fixed = TRUE)
  expect_match(src, 'input$ef_seed', fixed = TRUE)
  expect_match(src, 'process$kill_tree()', fixed = TRUE)
  expect_match(src, 'paste0("workers = ", workers)', fixed = TRUE)
  expect_match(src, 'paste0("seed = ", seed)', fixed = TRUE)
  expect_match(src, '0.10 + 0.84 * fraction', fixed = TRUE)
  expect_match(src, '"full person bootstrap" = 0.50 + 0.44 * fraction',
               fixed = TRUE)
  expect_false(grepl('"Hybrid (fast)"', src, fixed = TRUE))

  # Inference defaults and structural remedies are deliberately conservative.
  expect_match(src, '"Judge bootstrap (recommended)" = "judge_bootstrap"',
               fixed = TRUE)
  expect_match(src, '"Parametric bootstrap" = "bootstrap"', fixed = TRUE)
  expect_match(src, 'input$btlef_se %||% "judge_bootstrap"', fixed = TRUE)
  expect_match(src, 'input$btlef_workers', fixed = TRUE)
  expect_match(src, 'input$btlef_seed', fixed = TRUE)
  expect_match(src, 'identical(se_method, "judge_bootstrap")', fixed = TRUE)
  expect_match(src, 'Cancel frame estimation', fixed = TRUE)
  expect_match(src, 'btlef_job <- reactiveVal(NULL)', fixed = TRUE)
  expect_match(src, 'process$kill_tree()', fixed = TRUE)
  expect_match(src, '"boot_reps = %s, workers = %s, seed = %s, "',
               fixed = TRUE)
  expect_false(grepl('input$pc_id', src, fixed = TRUE))
  expect_match(src, 'unique(f$virtual_map$item)', fixed = TRUE)
  expect_match(src, "A location split cannot model", fixed = TRUE)
})

test_that("the app uses structurally stable responsive control layouts", {
  skip_if_not_installed("xml2")
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  old_cache <- Sys.getenv("R_USER_CACHE_DIR", unset = NA_character_)
  Sys.setenv(R_USER_CACHE_DIR = file.path(tempdir(), "r-cache"))
  on.exit(if (is.na(old_cache)) Sys.unsetenv("R_USER_CACHE_DIR") else
    Sys.setenv(R_USER_CACHE_DIR = old_cache), add = TRUE)
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  doc <- xml2::read_html(htmltools::renderTags(e$ui)$html)
  bad <- xml2::xml_find_all(doc,
    "//*[contains(concat(' ',normalize-space(@class),' '),' bslib-sidebar-layout ')]")
  expect_length(bad, 0L)
  top <- xml2::xml_find_all(doc,
    "//body/div[contains(@class,'container-fluid')]/div[contains(@class,'tab-content')]")
  expect_length(top, 1L)
  # Accordion headings are themselves buttons. An information button inside
  # one is invalid HTML and causes Chromium to close the tab container early.
  expect_length(xml2::xml_find_all(doc, "//button//button"), 0L)
  values <- xml2::xml_attr(xml2::xml_find_all(top,
    "./div[contains(concat(' ',normalize-space(@class),' '),' tab-pane ')]"),
    "data-value")
  expect_true(all(c("p_data", "p_targeting", "p_dif", "p_frames", "p_export")
                  %in% values))
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
