# A test that draws without choosing a file should use a null device. Leaving
# R's default file device in place creates tests/testthat/Rplots.pdf during an
# otherwise read-only test run.
options(device = function(...) grDevices::pdf(file = NULL))
