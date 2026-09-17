# Two panels of two outcomes each, measured on disjoint halves of the sites (no site has both panels).
# Truth: outcome 2 loads on the factor of outcome 1 (same panel), outcome 4 on the factor of outcome 3 (same panel),
# and outcome 3 on the factor of outcome 1 (across panels). Outcome 4 does not borrow from panel 1.
library(meshedghs)
set.seed(1)
n <- 900; q <- 4
coords <- cbind(runif(n), runif(n))                       # coordinates must be an unnamed matrix
D <- as.matrix(dist(coords))
V <- sapply(1:q, function(h) as.vector(t(chol(exp(-8 * D))) %*% rnorm(n)))
Lambda <- diag(q); Lambda[2, 1] <- 0.8; Lambda[4, 3] <- -0.7; Lambda[3, 1] <- 0.9
Y <- V %*% t(Lambda) + sqrt(0.05) * matrix(rnorm(n * q), n, q)
colnames(Y) <- c("a1", "a2", "b1", "b2")
panel_a <- sample(n, n / 2)
Y[panel_a, 3:4] <- NA; Y[-panel_a, 1:2] <- NA             # misalignment: NA marks "not measured here"

design_groups(c(1, 1, 2, 2))                              # the groups that spmeshed_ghs() will build

fit <- spmeshed_ghs(Y, x = matrix(1, n, 1), coords = coords,
                    axis_partition = c(5, 5), n_samples = 1000, n_burnin = 1500, n_thin = 1, n_threads = 2, verbose = 0,
                    prior = list(phi = c(0.5, 60)),
                    settings = list(adapting = TRUE, forced_grid = FALSE, cache = NULL, ps = TRUE, saving = FALSE, low_mem = TRUE, hmc = 4))
fit                                                       # blocks inferred from the missingness pattern; borrowing shares
borrowing_map(fit)                                        # decoupled selection: which blocks each outcome needs
round(zero_distance_cor(fit)$mean, 2)                     # posterior mean of cov2cor(Lambda Lambda')
round(cov2cor(tcrossprod(Lambda)), 2)                     # truth
