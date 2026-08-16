# Launch the rasch point-and-click graphical interface

Opens the guided Shiny application for users who prefer a graphical
workflow. The app supports data import; assignment of item, person,
group, rater, and comparison roles; model selection and fitting;
interactive diagnostics, tables, plots, and export. The corresponding R
call is shown for each analysis, but users do not need to write R
analysis code to use the principal workflows.

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
