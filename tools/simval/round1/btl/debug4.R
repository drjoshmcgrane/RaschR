suppressWarnings(pkgload::load_all(".", quiet=TRUE))
set.seed(42)
N <- 3000
delta <- c(I1 = -0.9, I2 = -0.1, I3 = 0.3, I4 = 0.7)
theta <- rnorm(N, 0, 1.4)
Y <- sapply(delta, function(dd) rbinom(N, 1, plogis(theta - dd)))
colnames(Y) <- names(delta)
pc <- pcml(Y)
print(pc$thr)
