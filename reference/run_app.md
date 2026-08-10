# Launch the rasch point-and-click graphical interface

Opens the guided Shiny application for users who prefer a graphical
workflow. The app supports data import; assignment of item, person,
group, rater, and comparison roles; model selection and fitting;
interactive diagnostics, tables, and plots; and one-click export. The
corresponding R call is shown for each analysis, but users do not need
to write R analysis code to use the principal workflows.

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

## Examples

``` r
if (interactive()) run_app()
```
