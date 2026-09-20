# meshedghs

**Adaptive borrowing across misaligned outcomes in meshed Gaussian process coregionalization models.**

`meshedghs` is an R package for Bayesian multivariate spatial and space-depth regression when different outcomes are measured at different locations. Typical cases are serial tissue sections stained with different antibody panels, a 3D volume in which marker panels alternate across depth, or monitoring networks that record different pollutants at different stations.

A joint model can predict an outcome where it was not measured by borrowing from the outcomes that were. Whether borrowing helps depends on the design. When two panels are measured far apart relative to the range of the process they share, the data carry almost no information about their dependence. The cross-loadings of an unrestricted joint model are then weakly identified, and their posterior dispersion adds to the predictive variance of every outcome. This package lets the data decide, block by block, how much each outcome borrows.

Status: version 0.1.3, pre-release. It accompanies a manuscript in preparation (S. Bhadury, *Design-adaptive borrowing across misaligned outcomes in Bayesian coregionalization models*). The simulation study of that manuscript is still running, so this README makes no claim about predictive gains. See [Limitations](#limitations).

## What the package contributes

1. **A design-based grouped horseshoe prior on the cross-loadings** of a linear model of coregionalization (LMC) whose factors are meshed Gaussian processes. The loadings of outcome *j* on the factors anchored in observation block *r* form one group with one local scale, so an outcome borrows from a whole block or reverts to univariate prediction with respect to it.
2. **The prior is placed on the identified loadings, under the parameter expansion used for sampling.** The meshed sampler draws expanded loadings whose scale is confounded with the factor variance. Placing a shrinkage prior on those would shrink a quantity with no meaning. The package derives the prior precision of the expanded loadings, the extra factor in the Metropolis ratio of the covariance parameters, and the sufficient statistic for the scales, for both the 2D Matérn and the 3D space-depth covariance.
3. **Exact draws of the local and global scales from discretized full conditionals** on a logarithmic grid, without auxiliary variables.
4. **Posterior inference on borrowing that is not based on shrinkage weights.** `borrowing_share()` gives the posterior of an effect size: the fraction of an outcome's latent spatial variance that is shared with a block. `borrowing_map()` gives a sparse map by decoupled shrinkage and selection (Hahn and Carvalho, 2015). Thresholding horseshoe shrinkage weights saturates when the number of locations is large, and is provided as a diagnostic only.
5. **An optional proper prior on the diagonal loadings**, for standardized outcomes on small domains, handled consistently with the parameter expansion.

The prior, its theory and the summaries of borrowing involve only the loading matrix and apply to any LMC; this package is the scalable implementation, in which the factors are meshed Gaussian processes. The meshed Gaussian process machinery (partitioning, the block Gibbs sampler for the latent factors, caching, prediction, non-Gaussian outcomes) is derived from the `meshed` package, version 0.2.3, by Michele Peruzzi (GPL >= 3), whose functions are retained. Michele Peruzzi was not involved in the development of this package. With no groups and no diagonal prior supplied, the code path is that of `meshed` 0.2.3.

Details of the model and of every change to the sampler are in [IMPLEMENTATION.md](IMPLEMENTATION.md).

## Installation

Requirements: R >= 4.1, a C++17 compiler with OpenMP (Rtools on Windows), and the packages `Rcpp`, `RcppArmadillo`, `dplyr`, `glue`, `rlang`, `magrittr`, `FNN`.

```r
# install.packages("remotes")
remotes::install_github("sagnikbhadury/meshedghs")   # while the repository is private, this needs a GitHub token (GITHUB_PAT)
```

or from a clone:

```
R CMD INSTALL meshedghs
```

On Windows, install to a library outside synchronized folders (Dropbox, OneDrive); they lock the compiled DLL.

## Quick start

Two panels of two outcomes each, measured on disjoint halves of 900 sites, so that no site has both panels. In the truth, `b1` loads on the factor of `a1` (across panels) and `b2` does not borrow from panel A. The script is `inst/examples/two_panels.R` and runs in about three minutes on two threads.

```r
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

fit <- spmeshed_ghs(Y, x = matrix(1, n, 1), coords = coords,
                    axis_partition = c(5, 5), n_samples = 1000, n_burnin = 1500, n_thin = 1, n_threads = 2, verbose = 0,
                    prior = list(phi = c(0.5, 60)),
                    settings = list(adapting = TRUE, forced_grid = FALSE, cache = NULL, ps = TRUE,
                                    saving = FALSE, low_mem = TRUE, hmc = 4))
fit
```

The blocks are inferred from the missingness pattern of the columns of `Y`. Printing the fit gives the posterior of the borrowing shares:

```
Meshed LMC with design horseshoe loadings: 4 outcomes, 4 groups, blocks 1 1 2 2
 group outcome borrows_from size     2.5%   50% 97.5% pr_above_eps
     1      a2           a1    1 1.10e-01 0.260 0.418        1.000
     2      b1       a1, a2    2 2.22e-01 0.442 0.692        1.000
     3      b2       a1, a2    2 4.24e-05 0.011 0.133        0.171
     4      b2           b1    1 2.04e-01 0.344 0.543        1.000
```

`b1` takes 44% of its latent variance from the factors of panel A (truth: 0.81 / 1.81 = 0.45), and `b2` takes 1% (truth: 0). The zero-distance correlations between panels are recovered although no site was measured on both:

```r
round(zero_distance_cor(fit)$mean, 2)     # posterior mean          # truth: cov2cor(Lambda %*% t(Lambda))
#       a1    a2    b1    b2              #       a1   a2    b1    b2
# a1  1.00  0.50  0.65 -0.10              # a1  1.00 0.62  0.67  0.00
# a2  0.50  1.00  0.32 -0.04              # a2  0.62 1.00  0.42  0.00
# b1  0.65  0.32  1.00 -0.51              # b1  0.67 0.42  1.00 -0.43
# b2 -0.10 -0.04 -0.51  1.00              # b2  0.00 0.00 -0.43  1.00
```

The sparse map applies a conservative rule: a group is dropped only if the posterior probability that dropping it gives up at most `eps = 0.05` of the latent variance is at least `level = 0.9`.

```r
borrowing_map(fit)
#   outcome group borrows_from retained path_score
# 1      a2     1           a1     TRUE 0.76773424
# 2      b1     2       a1, a2     TRUE 0.76319031
# 3      b2     3       a1, a2     TRUE 0.09219496
# 4      b2     4           b1     TRUE 0.74903936
```

In this short run the null group (`b2` from panel A) has `pr(share <= 0.05) = 0.83`, which is below 0.9, so it is retained; `borrowing_map(fit, level = 0.8)` drops it. The shares, and not the map, are the primary summary.

`fit$yhat_mcmc` holds posterior predictive draws at all input locations, including those where an outcome is `NA`, so predicting a panel where it was not measured requires no further step. Prediction at new locations uses `predict()` as in `meshed`.

## The model in brief

For outcomes $j = 1, \dots, q$ at locations $\ell$,

$$y_j(\ell) = x_j(\ell)^\top \beta_j + \lambda_{[j,:]}\, v(\ell) + \varepsilon_j(\ell), \qquad \varepsilon_j(\ell) \sim N(0, \tau_j^2),$$

where $v = (v_1, \dots, v_q)$ are independent meshed Gaussian processes with unit variance and $\Lambda$ is lower triangular. Factor $h$ is anchored at outcome $h$. Outcomes are ordered so that the observation blocks $B_1, \dots, B_R$ (outcomes measured at the same locations) are contiguous. For outcome $j$ and block $r$, the group $g = (j, r)$ collects the loadings $\lambda_{ji}$ with $i \in B_r$, $i < j$, and

$$\lambda_{ji} \mid \gamma_g, \zeta \sim N(0, \zeta^2 \gamma_g^2), \qquad \gamma_g \sim C^+(0, 1), \qquad \zeta \sim C^+(0, 1).$$

If all the loadings of outcome $j$ on other blocks are zero, its predictive distribution depends on the data of its own block only, so the model nests blockwise and univariate prediction. The borrowing share is

$$\omega_{j,r} = \frac{\sum_{i \in B_r,\, i<j} \lambda_{ji}^2}{\sum_{h \le j} \lambda_{jh}^2}.$$

Shares are directed: outcome $j$ can load only on the factors of outcomes that precede it. Put the panel you regard as the reference first.

## Functions

| Function | Purpose |
|---|---|
| `spmeshed_ghs(y, x, coords, blocks = NULL, grouping = "design", prior, ...)` | Fits the LMC with `k = q` factors and the grouped horseshoe. `grouping = "element"` gives one scale per loading (ordinary horseshoe), `"none"` the unrestricted LMC. All other arguments are those of `spmeshed()`. |
| `design_groups(blocks, grouping)` | The `q x q` integer matrix of group labels (0 = not shrunk). |
| `borrowing_share(fit, eps = 0.05)` | Posterior quantiles and draws of the borrowing shares, and `pr(share > eps)`. Works for any fit with `k = q`, under any prior. |
| `borrowing_map(fit, eps = 0.05, level = 0.9)` | Sparse map of borrowing by decoupled shrinkage and selection. |
| `zero_distance_cor(fit)` | Posterior of the cross-correlations at zero distance, `cov2cor(Lambda Lambda')`. |
| `shrinkage_weights(fit)` | Horseshoe shrinkage weights per group. Diagnostic only. |
| `spmeshed()`, `predict()`, `rmeshedgp()`, `summary_list_mean()`, `summary_list_q()` | Retained from `meshed` 0.2.3. `spmeshed()` accepts the new prior entries below. |

New entries of `prior`:

| Entry | Meaning | Default |
|---|---|---|
| `lambda_group` | `q x k` integer matrix of group labels, as returned by `design_groups()`. Set by `spmeshed_ghs()`. | all zero (no shrinkage) |
| `lambda_diag_sd` | Standard deviation `s0` of a `N(0, s0^2)` prior on the diagonal loadings, on the identified scale. Use about 3 for standardized outcomes on small domains, where the marginal variance of a long-range factor is weakly identified. | 0 (flat, as in `meshed`) |

New entries of the output: `ghs_lam_mcmc` (groups x draws, local scales $\gamma_g$), `ghs_tau_mcmc` (global scale $\zeta$), and `ghs` (group matrix, blocks, outcome names).

## Practical notes

- **Coordinates** must be an unnamed numeric matrix (`unname(as.matrix(coords))`). Named columns break the internal joins.
- **Three coordinates** select the space-depth (or space-time) covariance of Gneiting type, $\rho(h, u) = (a|u| + 1)^{-1} \exp\{-\phi \|h\| (a|u| + 1)^{-b/2}\}$. The prior on the depth decay `a` defaults to the bounds given for `phi`. That default is rarely appropriate: with depth in micrometers and `a <= 30`, adjacent sections are forced to be strongly correlated. Set `prior$a` explicitly on the scale of your depth coordinate.
- **Gridded locations are about ten times faster** than irregular ones, because the block geometry is cached. For dense imaging data, consider aggregating to a regular voxel grid.
- **Standardize the outcomes** if you use `lambda_diag_sd`, and in general for the half-Cauchy scales to be on a sensible scale.
- **Non-Gaussian outcomes** (`family = "poisson"`, `"binomial"`, `"beta"`): the scales are held at one during the first half of the burn-in. The gradient-based update of the loadings starts at zero and moves slowly, and the scales would otherwise collapse onto the starting value. This path has had limited testing.
- **Threads.** The sampler scales sublinearly beyond a few threads. For many fits, run several processes with two threads each.
- **Identification.** $\Lambda$ is identified up to the signs of its columns. All summaries provided here are invariant to them.

## Limitations

- The grouped horseshoe cannot create information. When the design carries none about the dependence between two panels, the best available outcome is reversion to univariate prediction. On serial sections of a colorectal cancer specimen with negligible dependence across sections, a univariate meshed model predicted better than every joint model, including this one.
- On a 3D breast tumor volume with panels alternating across depth, joint models improved on univariate prediction only when some layers were measured on both panels. With no overlap, the grouped horseshoe did not improve on the univariate model in prediction, and its 90% intervals undercovered (0.84). It did recover cross-panel correlations that were never observed at a common location.
- The theory behind the prior (posterior dispersion of order $\zeta\, I^{-1/2}$ for a weakly identified group, against $I^{-1}$ under a flat prior, where $I$ is the design information) concerns a local Gaussian experiment, with the global scale fixed. No contraction rate is claimed.
- `k = q` factors and a lower-triangular $\Lambda$ are assumed; the shares depend on the ordering of the outcomes.
- A controlled simulation study against the unrestricted LMC, univariate models, NNGP and the inside-out cross-covariance model is in progress. Results will be added here when it is complete.

## Citation

Until the manuscript is available, please cite the package and the methods it builds on:

> Bhadury, S. (2026). meshedghs: Adaptive borrowing across misaligned outcomes in meshed Gaussian process coregionalization models. R package version 0.1.3. https://github.com/sagnikbhadury/meshedghs

## References

- Peruzzi, M., Banerjee, S. and Finley, A. O. (2022). Highly scalable Bayesian geostatistical modeling via meshed Gaussian processes on partitioned domains. *Journal of the American Statistical Association*, 117, 969-982. doi:10.1080/01621459.2020.1833889
- Peruzzi, M., Banerjee, S., Dunson, D. B. and Finley, A. O. (2026). Gridding and parameter expansion for scalable latent Gaussian models of spatial multivariate data. *Bayesian Analysis*, 21, 1219-1245. doi:10.1214/25-BA1515
- Peruzzi, M. and Dunson, D. B. (2024). Spatial meshing for general Bayesian multivariate models. *Journal of Machine Learning Research*, 25(87), 1-49.
- Carvalho, C. M., Polson, N. G. and Scott, J. G. (2010). The horseshoe estimator for sparse signals. *Biometrika*, 97, 465-480.
- Xu, Z., Schmidt, D. F., Makalic, E., Qian, G. and Hopper, J. L. (2016). Bayesian grouped horseshoe regression with application to additive models. *AI 2016: Advances in Artificial Intelligence*, 229-240.
- Hahn, P. R. and Carvalho, C. M. (2015). Decoupling shrinkage and selection in Bayesian linear models: a posterior summary perspective. *Journal of the American Statistical Association*, 110, 435-448. doi:10.1080/01621459.2014.993077
- Bolfarine, H., Carvalho, C. M., Lopes, H. F. and Murray, J. S. (2024). Decoupling shrinkage and selection in Gaussian linear factor analysis. *Bayesian Analysis*, 19, 181-203. doi:10.1214/22-BA1349

## License

GPL (>= 3). The package is a derivative work of `meshed` 0.2.3, copyright Michele Peruzzi, distributed under GPL (>= 3).
