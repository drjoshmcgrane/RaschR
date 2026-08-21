# Launch the rasch point-and-click graphical interface

Opens the Shiny application for fitting models and examining their
tables, plots and diagnostics. The R code for each result is available
in the app. Analyses can be saved and reopened, or exported as HTML,
Word or PDF reports.

## Usage

``` r
run_app(...)
```

## Arguments

- ...:

  Passed to
  [`shiny::runApp`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Value

Called for its side effect of launching the app.

## Details

The app's display packages ('shiny', 'bslib', 'DT', 'bsicons') are
suggested rather than required by the package. If any are missing,
`run_app` lists them all and, in an interactive session, offers to
install them before launching.

## Examples

``` r
if (interactive()) run_app()
```
