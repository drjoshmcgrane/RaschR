# Screenshots of the graphical interface, for the vignettes and the README.
#
# Run from the package root:  Rscript tools/screenshots.R [name ...]
# With no arguments every shot is retaken; naming shots takes only those.
#
# The app is served from the in-tree inst/shiny by a background R process and
# driven through Chrome's debugging protocol, so a shot is a real render of
# the current source rather than a picture pasted in and forgotten. Each shot
# follows the path a reader would: pick the example dataset, estimate, open
# the panel, and only then capture. That is slower than poking values into
# the session, and it is the reason the images stay honest when the interface
# moves under them.

for (p in c("chromote", "callr", "shiny", "jsonlite"))
  if (!requireNamespace(p, quietly = TRUE))
    stop("tools/screenshots.R needs the suggested package '", p, "'")

`%||%` <- function(a, b) if (is.null(a)) b else a

# The app contains background processes of its own, so loading the source tree
# only in this coordinator is not enough: a child can otherwise find an older
# installed rasch on .libPaths(). Install the current tree into an isolated
# library and launch that exact installation. A failed installation stops the
# run before any existing screenshots can be overwritten.
ROOT <- normalizePath(".", mustWork = TRUE)
SHOT_LIB <- tempfile("rasch-screenshot-library-")
dir.create(SHOT_LIB, recursive = TRUE)
on.exit(unlink(SHOT_LIB, recursive = TRUE, force = TRUE), add = TRUE)
message("installing the current source tree for the screenshot run")
install_result <- callr::rcmd_safe(
  "INSTALL",
  c("--no-multiarch", "--with-keep.source", "-l", SHOT_LIB, ROOT),
  wd = ROOT, stdout = "|", stderr = "|", fail_on_status = FALSE)
if (!identical(install_result$status, 0L))
  stop("could not install the current source tree for screenshots:\n",
       paste(c(install_result$stdout, install_result$stderr), collapse = "\n"))

PORT   <- 7841L
WIDTH  <- 1440L
HEIGHT <- 900L
SCALE  <- 1                     # 1440 px remains sharp at vignette display width

# --- the shots -------------------------------------------------------------
# panel is the nav value; open names accordion panels to expand before the
# capture, now that they all start closed; before is any extra JS to run.
SHOTS <- list(
  list(name = "app-data", demo = "pcm", panel = "p_data", fit = FALSE,
       dir = "vignettes/figures",
       alt = "The Data panel: the sidebar assigns the person identifier, person factors and item columns, and the main area previews the responses."),
  list(name = "app-summary", demo = "pcm", panel = "p_summary", height = 1400L,
       open = "test_tcc", before = "document.getElementById('run_lr').click();",
       expect = c("#lr_txt" = "Adjusted chi-square"),
       expect_selector = "#tcc img",
       dir = "vignettes/figures"),
  list(name = "app-items", demo = "pcm", panel = "p_items", height = 1150L,
       table = "items_tbl", row = "I02", expect = c("#sel_item_title" = "I02"),
       dir = "vignettes/figures"),
  list(name = "app-items-chisq", demo = "pcm", panel = "p_items",
       tab = "Chi-square", table = "items_tbl", row = "I02",
       expect = c("#chisq_caption" = "I02"), dir = "vignettes/figures"),
  list(name = "app-persons", demo = "pcm", panel = "p_persons",
       open = "persons_pfit", dir = "vignettes/figures"),
  list(name = "app-targeting", demo = "pcm", panel = "p_targeting",
       height = 1600L,
       before = "var e=document.getElementById('tg_information'); if(e && !e.checked) e.click();",
       expect_checked = "tg_information",
       expect_selector = c("#pim_p img", "#wright img"),
       dir = "vignettes/figures"),
  list(name = "app-dif", demo = "pcm", panel = "p_dif",
       open = "dif_anova", table = "dif_tbl", row = "I08",
       expect = c("#dif_tbl table tbody tr.selected" = "I08"),
       expect_selector = "#dif_icc img", dir = "vignettes/figures"),
  list(name = "app-local", demo = "pcm", panel = "p_ld",
       open = "ld_cormat", height = 950L, dir = "vignettes/figures"),
  list(name = "app-rcode", demo = "pcm", panel = "p_data",
       open = "rcode_acc", focus = "#rcode_acc", focus_pad = 10L, height = 640L,
       expect = c("#rcode_fit" = "fit <-"), dir = "vignettes/figures"),
  list(name = "app-export", demo = "pcm", panel = "p_export",
       dir = "vignettes/figures"),
  # The README introduces the multiple-choice workflow separately. Keep its
  # discoverable miskey and ICC out of the polytomous vignette walkthrough.
  list(name = "readme-items", file = "app-items", demo = "dich",
       panel = "p_items", height = 1150L, table = "items_tbl", row = "I07",
       expect = c("#sel_item_title" = "I07"), dir = "man/figures"),
  # comparative judgement
  list(name = "app-cj-data", demo = "btl", panel = "p_data", fit = FALSE,
       dir = "vignettes/figures"),
  list(name = "app-cj-summary", demo = "btl", panel = "p_summary",
       dir = "vignettes/figures"),
  list(name = "app-cj-items", demo = "btl", panel = "p_items",
       open = "btl_caterpillar", focus = "[data-value='btl_caterpillar']",
       focus_pad = 15L, height = 950L, dir = "vignettes/figures"),
  list(name = "app-cj-judges", demo = "btl", panel = "p_persons",
       open = "btl_judge_fit",
       before = "var t=$('#btl_judges_tbl table').DataTable(); t.search('J36').draw();",
       wait_text = c("#btl_judges_tbl table tbody" = "J36"),
       table = "btl_judges_tbl", row = "J36",
       expect = c("#btl_judges_tbl table tbody tr.selected" = "J36"),
       expect_selector = "#btl_judge_map img", dir = "vignettes/figures"),
  list(name = "app-cj-dif", demo = "btl", panel = "p_dif",
       open = "bdif_anova",
       before = "document.getElementById('bdif_run').click();",
       wait_selector = "#bdif_anova_tbl table tbody tr",
       table = "bdif_anova_tbl", row = "O1",
       expect = c("#bdif_anova_tbl table tbody tr.selected" = "panel"),
       expect_selector = "#bdif_occ img", height = 1150L,
       dir = "vignettes/figures"))

# The model each example implies. A shot asserts this before capturing:
# switching example while an analysis is on screen swaps the data under a
# fitted model, and a capture taken during that swap shows the PREVIOUS
# analysis, mid-fade, with nothing to say it is the wrong one.
DEMO_MODEL <- c(dich = "rasch", pcm = "rasch", rsm = "rasch",
                mfrm = "mfrm", efrm = "efrm", btl = "btl")

# --- driving ---------------------------------------------------------------
js <- function(b, expr, wait = TRUE)
  b$Runtime$evaluate(expr, awaitPromise = wait, returnByValue = TRUE)$result$value

# Shiny sets a body class while any output is recomputing. Waiting on that,
# rather than on a fixed sleep, is what keeps a shot from catching a spinner.
# Only VISIBLE outputs count: Shiny marks an output .recalculating when it is
# first bound and clears it when it renders, so outputs in panels that have
# never been opened keep the class indefinitely -- 150 of them here, which is
# never zero and turns every wait into its own timeout.
# The timeout is generous because it now measures real work: a comparative
# judgement fit runs for minutes, and a capture taken while it is still going
# is a blank panel under a correct-looking header.
settle <- function(b, timeout = 600) {
  deadline <- Sys.time() + timeout
  quiet <- 0L
  while (Sys.time() < deadline) {
    busy <- isTRUE(js(b, "document.body.classList.contains('shiny-busy') ||
      Array.from(document.querySelectorAll('.recalculating'))
           .some(function(e) { return e.offsetParent !== null })"))
    quiet <- if (busy) 0L else quiet + 1L
    if (quiet >= 8L) return(invisible(TRUE))
    Sys.sleep(0.25)
  }
  stop("the application was still busy after ", timeout, "s", call. = FALSE)
}

# Estimate is an input_task_button: the fit runs in an ExtendedTask, which does
# NOT raise Shiny's busy flag, so waiting on that alone returns to an
# apparently idle page and captures an empty panel under a header that already
# looks right. Wait for the analysis itself to appear -- and say so when it
# does not, rather than sitting out the whole timeout in silence.
rendered <- function(b, what, timeout = 600, chars = 600) {
  deadline <- Sys.time() + timeout
  repeat {
    # the whole page's text, not a pane's: every selector for "the panel" in
    # this layout matched something empty. A blank panel leaves only the
    # navigation bar behind, a few hundred characters below anything real.
    n <- js(b, "document.body.innerText.trim().length")
    if (isTRUE(n >= chars)) return(invisible(TRUE))
    if (Sys.time() > deadline) {
      stop(what, ": the page still holds only ", n, " characters after ",
           timeout, "s", call. = FALSE)
    }
    Sys.sleep(0.5)
  }
}

wait_for <- function(b, selector, timeout = 60) {
  deadline <- Sys.time() + timeout
  while (Sys.time() < deadline) {
    if (isTRUE(js(b, sprintf("document.querySelector(%s) !== null",
                             shQuote(selector, type = "cmd"))))) return(TRUE)
    Sys.sleep(0.25)
  }
  stop("timed out waiting for ", selector)
}

wait_for_text <- function(b, selector, expected, timeout = 60) {
  deadline <- Sys.time() + timeout
  while (Sys.time() < deadline) {
    got <- js(b, sprintf(
      "(function(){var e=document.querySelector(%s); return e ? e.innerText : ''})()",
      shQuote(selector, type = "cmd")))
    if (is.character(got) && grepl(expected, got, fixed = TRUE)) return(TRUE)
    Sys.sleep(0.25)
  }
  stop("timed out waiting for ", shQuote(expected), " in ", selector)
}

current_demo <- NULL
current_fit_demo <- NULL

assert_clean <- function(b, shot) {
  errors <- js(b, "Array.from(document.querySelectorAll('.shiny-output-error'))
    .filter(function(e){return e.offsetParent !== null})
    .map(function(e){return e.innerText.trim()}).filter(Boolean)")
  if (length(errors))
    stop(shot$name, ": visible Shiny error: ", paste(errors, collapse = " | "))
  busy <- js(b, "document.body.classList.contains('shiny-busy') ||
    Array.from(document.querySelectorAll('.recalculating'))
      .some(function(e){return e.offsetParent !== null})")
  if (isTRUE(busy)) stop(shot$name, ": capture attempted while outputs were busy")
  invisible(TRUE)
}

assert_text <- function(b, selector, expected, shot) {
  got <- js(b, sprintf(
    "(function(){var e=document.querySelector(%s); return e ? e.innerText.trim() : null})()",
    shQuote(selector, type = "cmd")))
  if (is.null(got) || !grepl(expected, got, fixed = TRUE))
    stop(shot$name, ": expected ", shQuote(expected), " in ", selector,
         "; found ", shQuote(got %||% "<missing>"))
}

take <- function(b, shot) {
  message("  ", shot$name)
  want <- DEMO_MODEL[[shot$demo]]
  # Changing example on top of a fitted analysis leaves the old fit on screen
  # while the new data loads. Reload instead: the session comes back at the
  # welcome screen with nothing to mistake for the new analysis.
  if (!identical(current_demo, shot$demo)) {
    if (!is.null(current_demo)) {
      js(b, "location.reload()", wait = FALSE)
      current_fit_demo <<- NULL
      Sys.sleep(2)
      wait_for(b, "#demo_choice", timeout = 90)
      # the input exists before the websocket does; wait for Shiny itself
      deadline <- Sys.time() + 60
      while (Sys.time() < deadline &&
             !isTRUE(js(b, "typeof Shiny !== 'undefined' && Shiny.shinyapp &&
                            Shiny.shinyapp.$socket &&
                            Shiny.shinyapp.$socket.readyState === 1")))
        Sys.sleep(0.5)
      settle(b)
    }
    changed <- js(b, sprintf(
      "(function(){var e=document.getElementById('demo_choice'); if(!e) return false;
        if(e.selectize) e.selectize.setValue('%s');
        else {e.value='%s'; e.dispatchEvent(new Event('change',{bubbles:true}));}
        return true})()", shot$demo, shot$demo))
    if (!isTRUE(changed)) stop(shot$name, ": the example selector was not found")
    settle(b)
    # the example selects its model; wait for that to arrive before estimating
    deadline <- Sys.time() + 30
    repeat {
      got <- js(b, "(function(){var e = document.querySelector('input[name=model_type]:checked'); return e ? e.value : ''})()")
      if (identical(got, want) || Sys.time() > deadline) break
      Sys.sleep(0.25)
    }
    current_demo <<- shot$demo
  }
  if (!identical(shot$fit, FALSE) &&
      !identical(current_fit_demo, shot$demo)) {
    js(b, "document.getElementById('run').click()")
    Sys.sleep(1)
    # Estimate is an input_task_button: the fit runs in an ExtendedTask, which
    # does NOT raise Shiny's busy flag, so waiting on that alone returns to an
    # apparently idle page and captures an empty panel under a header that
    # already looks right. Wait for the analysis itself to appear.
    wait_for(b, "#nav_status .rasch-nav-summary", timeout = 600)
    settle(b)
    current_fit_demo <<- shot$demo
  }
  got <- js(b, "(function(){var e = document.querySelector('input[name=model_type]:checked'); return e ? e.value : ''})()")
  if (!identical(got, want))
    stop(shot$name, ": the application is showing the ", got, " model, not ",
         want, " -- the example did not take")
  shown_demo <- js(b, "(function(){var e=document.getElementById('demo_choice');
    if(!e) return null; return e.selectize ? e.selectize.getValue() : e.value})()")
  if (!identical(shown_demo, shot$demo))
    stop(shot$name, ": the visible example selector did not retain ", shot$demo)
  js(b, sprintf(
    "var l = document.querySelector('a[data-value=\"%s\"]'); if (l) l.click();",
    shot$panel))
  # a panel reached through a nav menu leaves that menu hanging open over the
  # page, and neither a body click nor a click elsewhere dismisses it
  js(b, "document.querySelectorAll('.dropdown-menu.show').forEach(function(m){ m.classList.remove('show'); });
         document.querySelectorAll('.dropdown-toggle.show, [aria-expanded=\"true\"].dropdown-toggle').forEach(function(t){ t.classList.remove('show'); t.setAttribute('aria-expanded', 'false'); });")
  settle(b)
  # `open` names either an accordion (by id) or one of its panels (by the
  # data-value bslib puts on the item) -- resolve both, and expand the first
  # collapsed button found inside
  # the panel just navigated to must have rendered before anything is opened
  # inside it or captured from it
  # Plot-led panels contain little text even when fully drawn. Output activity
  # is handled by settle(); panel-specific image assertions below verify the
  # plots that must be present in the captured state.
  rendered(b, shot$panel, chars = 200)
  if (!is.null(shot$open)) {
    js(b, sprintf(
      "var p = document.querySelector('#%s') || document.querySelector('[data-value=\"%s\"]');
       if (p) { var btn = p.querySelector('.accordion-button.collapsed'); if (btn) btn.click(); }",
      shot$open, shot$open))
    settle(b)
  }
  if (!is.null(shot$before)) js(b, shot$before)
  if (!is.null(shot$wait_text))
    for (selector in names(shot$wait_text))
      wait_for_text(b, selector, unname(shot$wait_text[[selector]]),
                    timeout = 600)
  if (!is.null(shot$wait_selector)) {
    wait_for(b, shot$wait_selector, timeout = 600)
    settle(b)
  }
  # a table that drives a plot beside it needs the interesting row selected,
  # or the figure shows whichever row the panel opened on
  if (!is.null(shot$row)) {
    target <- js(b, sprintf(
      "(function(){var root=document.getElementById(%s); if(!root) return false;
       var r=Array.from(root.querySelectorAll('table tbody tr')).find(function(t){
         return Array.from(t.querySelectorAll('td')).some(function(c){return c.textContent.trim() === %s;});
       }); if(!r) return null; r.scrollIntoView({block:'center'});
       var q=r.getBoundingClientRect();
       return {selected:r.classList.contains('selected'),
               x:q.left + Math.min(20, q.width/2), y:q.top + q.height/2}})()",
      shQuote(shot$table, type = "cmd"), shQuote(shot$row, type = "cmd")))
    if (is.null(target) || !is.finite(target$x) || !is.finite(target$y))
      stop(shot$name, ": could not select ", shot$row, " in ", shot$table)
    if (!isTRUE(target$selected)) {
      b$Input$dispatchMouseEvent(type = "mousePressed", x = target$x,
                                 y = target$y, button = "left", clickCount = 1)
      b$Input$dispatchMouseEvent(type = "mouseReleased", x = target$x,
                                 y = target$y, button = "left", clickCount = 1)
    }
    settle(b)
  }
  # Select the item before changing the detail tab: switching tabs can move
  # the table while a coordinate-based browser click is in flight.
  if (!is.null(shot$tab)) {
    js(b, sprintf(
      "var t = Array.from(document.querySelectorAll('.nav-link')).find(function(e){return e.textContent.trim() === %s}); if (t) t.click();",
      shQuote(shot$tab, type = "cmd")))
    settle(b)
  }
  settle(b)
  if (!is.null(shot$expect))
    for (selector in names(shot$expect))
      assert_text(b, selector, unname(shot$expect[[selector]]), shot)
  if (!is.null(shot$expect_selector)) {
    for (selector in shot$expect_selector) {
      present <- js(b, sprintf("document.querySelector(%s) !== null",
                               shQuote(selector, type = "cmd")))
      if (!isTRUE(present))
        stop(shot$name, ": expected selector ", selector)
    }
  }
  if (!is.null(shot$expect_checked)) {
    checked <- js(b, sprintf(
      "(function(){var e=document.getElementById(%s); return !!(e && e.checked)})()",
      shQuote(shot$expect_checked, type = "cmd")))
    if (!isTRUE(checked))
      stop(shot$name, ": expected checked input ", shot$expect_checked)
  }
  assert_clean(b, shot)
  # a fit raises toast notifications ("15 item(s) scored 0/1 against the key")
  # that float over the panel and would sit in the middle of the figure
  js(b, "var n = document.getElementById('shiny-notification-panel'); if (n) n.remove();")
  js(b, "window.scrollTo(0, 0);")
  Sys.sleep(0.75)
  # The whole page is far too tall to read at a vignette's width -- the Items
  # panel alone runs to nearly 6,000 pixels. Capture what a reader sees: the
  # viewport, or the taller frame a shot asks for.
  h <- shot$height %||% HEIGHT
  focus_pad <- shot$focus_pad %||% 80L
  y <- if (is.null(shot$focus)) 0 else js(b, sprintf(
    "(function(){var e=document.querySelector(%s); return e ? Math.max(0, e.getBoundingClientRect().top + window.scrollY - %d) : 0})()",
    shQuote(shot$focus, type = "cmd"), focus_pad))
  raw <- b$Page$captureScreenshot(
    format = "png", captureBeyondViewport = TRUE,
    clip = list(x = 0, y = y, width = WIDTH, height = h, scale = SCALE))$data
  bin <- jsonlite::base64_dec(raw)
  for (d in shot$dir) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    writeBin(bin, file.path(d, paste0(shot$file %||% shot$name, ".png")))
  }
  invisible(TRUE)
}

# --- run -------------------------------------------------------------------
args <- commandArgs(TRUE)
shots <- if (length(args)) Filter(function(s) s$name %in% args, SHOTS) else SHOTS
if (!length(shots)) stop("no shot matched: ", paste(args, collapse = ", "))

app <- callr::r_bg(function(lib, port) {
  .libPaths(c(lib, .libPaths()))
  library(rasch, lib.loc = lib)
  dir <- system.file("shiny", package = "rasch", lib.loc = lib)
  if (!nzchar(dir)) stop("the installed package does not contain the Shiny app")
  options(shiny.port = port, shiny.host = "127.0.0.1",
          shiny.launch.browser = FALSE, shiny.maxRequestSize = 100 * 1024^2)
  shiny::runApp(dir)
}, args = list(lib = SHOT_LIB, port = PORT),
  supervise = TRUE, stdout = "|", stderr = "|")
on.exit({ try(app$kill(), silent = TRUE) }, add = TRUE)

message("starting the application on port ", PORT)
for (i in 1:120) {
  if (!app$is_alive())
    stop("the application exited:\n", paste(app$read_all_error_lines(), collapse = "\n"))
  ok <- tryCatch({
    con <- suppressWarnings(socketConnection("127.0.0.1", PORT, open = "r+",
                                             blocking = TRUE, timeout = 1))
    close(con); TRUE
  }, error = function(e) FALSE)
  if (ok) break
  Sys.sleep(0.5)
}

options(chromote.timeout = 120)
b <- chromote::ChromoteSession$new(width = WIDTH, height = HEIGHT)
on.exit(try(b$close(), silent = TRUE), add = TRUE)
b$Page$navigate(sprintf("http://127.0.0.1:%d", PORT))
b$Page$loadEventFired()
wait_for(b, ".shiny-bound-input, #demo_choice")
settle(b)

message("taking ", length(shots), " shot(s)")
for (s in shots) take(b, s)
message("done")
