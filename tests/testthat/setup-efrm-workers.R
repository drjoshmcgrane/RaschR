# Keep routine tests serial. Parallel behaviour has its own installed-package
# integration test.
options(rasch.max_workers = 1L)
