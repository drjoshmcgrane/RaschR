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

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
for (p in c("chromote", "callr", "shiny", "jsonlite"))
  if (!requireNamespace(p, quietly = TRUE))
    stop("tools/screenshots.R needs the suggested package '", p, "'")

PORT   <- 7841L
WIDTH  <- 1440L
HEIGHT <- 900L
SCALE  <- 2                     # retina: the images are read at 90-100% width

# --- the shots -------------------------------------------------------------
# panel is the nav value; open names accordion panels to expand before the
# capture, now that they all start closed; before is any extra JS to run.
SHOTS <- list(
  list(name = "app-data", demo = "dich", panel = "p_data", fit = FALSE,
       dir = "vignettes/figures",
       alt = "The Data panel: the sidebar assigns the person identifier, person factors and item columns, and the main area previews the responses."),
  list(name = "app-summary", demo = "dich", panel = "p_summary", height = 1150L,
       dir = "vignettes/figures"),
  list(name = "app-items", demo = "dich", panel = "p_items", height = 1150L,
       row = "I07", dir = c("vignettes/figures", "man/figures")),
  list(name = "app-items-chisq", demo = "dich", panel = "p_items",
       tab = "Chi-square", row = "I07", dir = "vignettes/figures"),
  list(name = "app-persons", demo = "dich", panel = "p_persons",
       open = "persons_pfit", dir = "vignettes/figures"),
  list(name = "app-targeting", demo = "dich", panel = "p_targeting",
       dir = "vignettes/figures"),
  list(name = "app-dif", demo = "dich", panel = "p_dif",
       open = "dif_anova", row = "I05", dir = "vignettes/figures"),
  list(name = "app-local", demo = "dich", panel = "p_ld",
       open = "ld_cormat", dir = "vignettes/figures"),
  list(name = "app-rcode", demo = "dich", panel = "p_summary",
       open = "rcode_acc", dir = "vignettes/figures"),
  list(name = "app-export", demo = "dich", panel = "p_export",
       dir = "vignettes/figures"),
  # comparative judgement
  list(name = "app-cj-data", demo = "btl", panel = "p_data", fit = FALSE,
       dir = "vignettes/figures"),
  list(name = "app-cj-summary", demo = "btl", panel = "p_summary",
       dir = "vignettes/figures"),
  list(name = "app-cj-items", demo = "btl", panel = "p_items",
       open = "btl_caterpillar", dir = "vignettes/figures"),
  list(name = "app-cj-judges", demo = "btl", panel = "p_persons",
       open = "btl_judge_fit", dir = "vignettes/figures"),
  list(name = "app-cj-dif", demo = "btl", panel = "p_dif",
       open = "bdif_anova", dir = "vignettes/figures"))

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
  warning("still busy after ", timeout, "s", call. = FALSE)
  invisible(FALSE)
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
      warning(what, ": the page still holds only ", n, " characters after ",
              timeout, "s", call. = FALSE)
      return(invisible(FALSE))
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

current_demo <- NULL

take <- function(b, shot) {
  message("  ", shot$name)
  want <- DEMO_MODEL[[shot$demo]]
  # Changing example on top of a fitted analysis leaves the old fit on screen
  # while the new data loads. Reload instead: the session comes back at the
  # welcome screen with nothing to mistake for the new analysis.
  if (!identical(current_demo, shot$demo)) {
    if (!is.null(current_demo)) {
      js(b, "location.reload()", wait = FALSE)
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
    js(b, sprintf("Shiny.setInputValue('demo_choice', '%s', {priority: 'event'})",
                  shot$demo))
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
  if (!identical(shot$fit, FALSE)) {
    js(b, "document.getElementById('run').click()")
    Sys.sleep(1)
    # Estimate is an input_task_button: the fit runs in an ExtendedTask, which
    # does NOT raise Shiny's busy flag, so waiting on that alone returns to an
    # apparently idle page and captures an empty panel under a header that
    # already looks right. Wait for the analysis itself to appear.
    rendered(b, "the fit")
    settle(b)
  }
  got <- js(b, "(function(){var e = document.querySelector('input[name=model_type]:checked'); return e ? e.value : ''})()")
  if (!identical(got, want))
    stop(shot$name, ": the application is showing the ", got, " model, not ",
         want, " -- the example did not take")
  js(b, sprintf(
    "var l = document.querySelector('a[data-value=\"%s\"]'); if (l) l.click();",
    shot$panel))
  # a panel reached through a nav menu leaves that menu hanging open over the
  # page, and neither a body click nor a click elsewhere dismisses it
  js(b, "document.querySelectorAll('.dropdown-menu.show').forEach(function(m){ m.classList.remove('show'); });
         document.querySelectorAll('.dropdown-toggle.show, [aria-expanded=\"true\"].dropdown-toggle').forEach(function(t){ t.classList.remove('show'); t.setAttribute('aria-expanded', 'false'); });")
  settle(b)
  if (!is.null(shot$tab))
    js(b, sprintf(
      "var t = Array.from(document.querySelectorAll('.nav-link')).find(function(e){return e.textContent.trim() === %s}); if (t) t.click();",
      shQuote(shot$tab, type = "cmd")))
  # `open` names either an accordion (by id) or one of its panels (by the
  # data-value bslib puts on the item) -- resolve both, and expand the first
  # collapsed button found inside
  # the panel just navigated to must have rendered before anything is opened
  # inside it or captured from it
  rendered(b, shot$panel)
  if (!is.null(shot$open))
    js(b, sprintf(
      "var p = document.querySelector('#%s') || document.querySelector('[data-value=\"%s\"]');
       if (p) { var btn = p.querySelector('.accordion-button.collapsed'); if (btn) btn.click(); }",
      shot$open, shot$open))
  # a table that drives a plot beside it needs the interesting row selected,
  # or the figure shows whichever row the panel opened on
  if (!is.null(shot$row)) {
    js(b, sprintf(
      "var r = Array.from(document.querySelectorAll('table tbody tr')).find(function(t){ var c = t.querySelector('td'); return c && c.textContent.trim() === %s; });
       if (r) r.click();", shQuote(shot$row, type = "cmd")))
    settle(b)
  }
  settle(b)
  # a fit raises toast notifications ("15 item(s) scored 0/1 against the key")
  # that float over the panel and would sit in the middle of the figure
  js(b, "var n = document.getElementById('shiny-notification-panel'); if (n) n.remove();")
  js(b, "window.scrollTo(0, 0);")
  Sys.sleep(0.75)
  # The whole page is far too tall to read at a vignette's width -- the Items
  # panel alone runs to nearly 6,000 pixels. Capture what a reader sees: the
  # viewport, or the taller frame a shot asks for.
  h <- shot$height %||% HEIGHT
  raw <- b$Page$captureScreenshot(
    format = "png", captureBeyondViewport = TRUE,
    clip = list(x = 0, y = 0, width = WIDTH, height = h, scale = SCALE))$data
  bin <- jsonlite::base64_dec(raw)
  for (d in shot$dir) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    writeBin(bin, file.path(d, paste0(shot$name, ".png")))
  }
  invisible(TRUE)
}

# --- run -------------------------------------------------------------------
args <- commandArgs(TRUE)
shots <- if (length(args)) Filter(function(s) s$name %in% args, SHOTS) else SHOTS
if (!length(shots)) stop("no shot matched: ", paste(args, collapse = ", "))

app <- callr::r_bg(function(dir, port) {
  options(shiny.port = port, shiny.host = "127.0.0.1",
          shiny.launch.browser = FALSE, shiny.maxRequestSize = 100 * 1024^2)
  shiny::runApp(dir)
}, args = list(dir = normalizePath("inst/shiny"), port = PORT),
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
