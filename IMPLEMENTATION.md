# Implementation notes

This document states the model, derives what changes in the sampler when a grouped horseshoe prior is placed on the loadings of a meshed coregionalization model, and maps each piece to the source code. Symbols follow the manuscript. The global scale is $\zeta$ and the local scales are $\gamma_g$, because $\lambda$ and $\tau$ are taken by the loadings and the measurement error.

Contents: [1. Model](#1-model) · [2. Groups](#2-design-based-groups) · [3. Parameter expansion](#3-parameter-expansion-and-why-it-matters) · [4. Sampler](#4-the-sampler) · [5. Scale updates](#5-scale-updates-on-a-grid) · [6. Diagonal prior](#6-prior-on-the-diagonal-loadings) · [7. Non-Gaussian outcomes](#7-non-gaussian-outcomes) · [8. Posterior summaries](#8-posterior-summaries-of-borrowing) · [9. Source map](#9-source-map) · [10. Validation](#10-validation) · [11. Known issues](#11-known-issues)

## 1. Model

For outcomes $j = 1, \dots, q$ observed at sets of locations $T_j \subset D \subset \mathbb{R}^d$, $d \in \{2, 3\}$, which need not coincide,

$$y_j(\ell) = x_j(\ell)^\top \beta_j + \sum_{h \le j} \lambda_{jh} v_h(\ell) + \varepsilon_j(\ell), \qquad \varepsilon_j(\ell) \sim N(0, \tau_j^2).$$

The factors $v_1, \dots, v_k$, with $k = q$, are independent meshed Gaussian processes (Peruzzi, Banerjee and Finley, 2022) on an axis-parallel partition of the domain with a cubic mesh, each with unit variance and correlation parameters $\theta_h$. $\Lambda$ is lower triangular, the Cholesky-type identification of the zero-distance covariance $C(0) = \Lambda\Lambda^\top$, up to the signs of the columns. Factor $h$ is anchored at outcome $h$: it is the first factor that outcome $h$ loads on and no earlier outcome does.

Correlation functions. In $d = 2$, Matérn with smoothness $\nu \in \{1/2, 3/2, 5/2\}$ and decay $\phi_h$. In $d = 3$, with the third coordinate $u$ read as depth or time,

$$\rho_h\{(s, u), (s', u')\} = \frac{1}{a_h |u - u'| + 1} \exp\Big\{-\frac{\phi_h \|s - s'\|}{(a_h |u - u'| + 1)^{b_h/2}}\Big\}.$$

Priors: $\beta_j$ Gaussian, $\tau_j^2$ inverse gamma, $\theta_h$ uniform on a box, as in `meshed`.

## 2. Design-based groups

Let $B_1, \dots, B_R$ partition the outcomes into observation blocks: outcomes in the same block are measured at the same locations. Outcomes are ordered so that blocks are contiguous. For each outcome $j \ge 2$ and each block $r$ that contains an index below $j$, the group $g = (j, r)$ is

$$\{\lambda_{ji} : i \in B_r,\ i < j\}, \qquad d_g = \#\{i \in B_r : i < j\}.$$

Groups within the block of outcome $j$ itself are included; they are strongly identified and the prior has little effect on them. The prior is

$$\lambda_{ji} \mid \gamma_g, \zeta \sim N(0, \zeta^2\gamma_g^2), \qquad \gamma_g \sim C^+(0,1), \qquad \zeta \sim C^+(0,1).$$

`design_groups(blocks, grouping)` returns the $q \times q$ integer matrix $G$ with $G_{ji} = g$ for a shrunk loading and $0$ otherwise. For `blocks = c(1, 1, 2, 2)`:

```
     [,1] [,2] [,3] [,4]
[1,]    0    0    0    0
[2,]    1    0    0    0      group 1: a2 on the factor of a1 (within panel A)
[3,]    2    2    0    0      group 2: b1 on the factors of panel A (across panels)
[4,]    3    3    4    0      group 3: b2 on panel A; group 4: b2 on the factor of b1
```

`grouping = "element"` gives every loading its own group, which is the ordinary horseshoe; `"none"` gives the unrestricted LMC.

If every loading of outcome $j$ on other blocks is zero, then conditionally on the parameters the likelihood factorizes over blocks and the predictive distribution of outcome $j$ depends only on the data of its block. The model therefore nests blockwise and univariate prediction, and a group whose local scale collapses takes outcome $j$ back to it.

## 3. Parameter expansion, and why it matters

`meshed` samples in an expanded parametrization (Peruzzi, Banerjee, Dunson and Finley, 2025). Each factor carries a working variance $\omega_h = \omega(\theta_h)$ that is not identified, the sampler works with $v_h^* = \omega_h^{1/2} v_h$ and

$$\lambda^*_{jh} = \lambda_{jh}\, \omega_h^{-1/2},$$

and the identified quantities are recovered after sampling. The working variance is $\omega_h = \sigma_h^2$ for the 3D covariance and $\omega_h = \sigma_h^2 \phi_h^{-2\nu}$ for the 2D Matérn; in the code it is the squared diagonal of `ps_back(theta, d, twonu, using_ps)`.

In `meshed`, the prior on the loadings is flat (precision $10^{-6}$), so the distinction between $\lambda$ and $\lambda^*$ does not matter. With a shrinkage prior it does. A horseshoe on $\lambda^*$ would shrink a quantity whose scale is arbitrary and changes along the chain. We place the prior on the identified $\lambda$ and transform it. With $\theta$ a priori independent of $\Lambda$:

**(i) Prior of the expanded loadings.** $\lambda^*_{ji} \mid \theta_i, \gamma_g, \zeta \sim N(0, \zeta^2\gamma_g^2/\omega_i)$. The Jacobian $\omega_i^{1/2}$ of the change of variables is absorbed into this density, and the marginal prior of $\theta$ is unchanged.

**(ii) Update of the loadings.** For a Gaussian outcome, the full conditional of $(\beta_j, \lambda^*_{[j,:]})$ is Gaussian with precision

$$\tau_j^{-2} [X_j\ V_j^*]^\top [X_j\ V_j^*] + \mathrm{diag}(V_\beta^{-1}, D_j), \qquad (D_j)_{ii} = \omega_i / (\zeta^2 \gamma_{g(j,i)}^2)\ \text{for}\ i < j.$$

**(iii) Update of the covariance parameters.** Because the prior of $\lambda^*$ depends on $\theta$ through $\omega$, the Metropolis ratio for $\theta$ acquires the factor

$$\prod_{i} \omega_i^{m_i/2} \exp\Big\{-\frac{\omega_i}{2} \sum_{j > i} \frac{(\lambda^*_{ji})^2}{\zeta^2 \gamma_{g(j,i)}^2}\Big\},$$

where $m_i$ is the number of shrunk loadings in column $i$. Without this factor the prior implied on the identified loadings is no longer exactly the grouped horseshoe. In a sampler for the prior alone, omitting it displaced the quantiles of $|\lambda|$ by about 5%.

**(iv) Sufficient statistic for the scales.** The scales depend on the loadings through the identified quantity $S_g = \sum_{(j,i) \in g} (\lambda^*_{ji})^2 \omega_i = \sum_{(j,i) \in g} \lambda_{ji}^2$.

An error that is easy to make, and that an early version of this code made: using $\sigma^2$ as the working variance in $d = 2$. It is $\sigma^2\phi^{-2\nu}$. All four items above now obtain $\omega$ from the same function, `ghs_sigmasq()`.

## 4. The sampler

One iteration, with new or changed steps in bold:

1. For each outcome $j$: draw $(\beta_j, \lambda^*_{[j,:]})$ from the Gaussian full conditional, **with prior precision (ii)**; for non-Gaussian outcomes, the gradient-based update of `meshed` with the same prior precision.
2. **Draw every $\gamma_g$, then $\zeta$, from their discretized full conditionals given $S_g$ (Section 5).**
3. Draw $\tau_j^2$.
4. Draw the latent factors $v^*$ block by block over the mesh, in parallel within each color of the partition.
5. Update $\theta$ by robust adaptive Metropolis, **with the factor (iii) in the acceptance ratio**.

Steps 3 to 5 are those of `meshed` apart from the acceptance ratio. The cost of step 2 is $O(G \times 2001)$ per iteration, which is negligible next to step 4.

## 5. Scale updates on a grid

With $u = \log\gamma_g$,

$$\pi(u \mid \cdot) \propto (\zeta^2 e^{2u})^{-d_g/2} \exp\Big\{-\frac{S_g}{2\zeta^2 e^{2u}}\Big\} \frac{e^{u}}{1 + e^{2u}},$$

where the last term is the half-Cauchy density of $\gamma_g$ times the Jacobian $e^u$. With $t = \log\zeta$,

$$\pi(t \mid \cdot) \propto e^{-t\sum_g d_g} \exp\Big\{-e^{-2t} \sum_g \frac{S_g}{2\gamma_g^2}\Big\} \frac{e^{t}}{1 + e^{2t}}.$$

Both are evaluated on 2001 equally spaced points of $[-12, 12]$ and sampled by inversion of the discrete distribution function (`ghs_grid_draw()`). Log weights are shifted by their maximum before exponentiation. The grid spacing is 0.012 on the log scale, and the range $e^{\pm 12}$ covers scales from $6 \times 10^{-6}$ to $1.6 \times 10^{5}$. The draw is exact for the discretized distribution; there are no auxiliary variables and no rejection step, and mixing of the scales does not depend on a slice or an inverse-gamma augmentation.

## 6. Prior on the diagonal loadings

The diagonal loadings carry the marginal variances and are not shrunk toward zero. With a flat prior, and a factor whose range is long relative to the domain, the marginal variance is weakly identified (Zhang, 2004) and can drift to large values; this happened on a 25-section tissue dataset, with implied variances up to 27 for standardized outcomes. `prior$lambda_diag_sd = s0` places $\lambda_{jj} \sim N(0, s_0^2)$ on the identified scale. The same argument as in Section 3 applies with $s_0^2$ in place of $\zeta^2\gamma_g^2$: the prior precision of $\lambda^*_{jj}$ is $\omega_j/s_0^2$, and the $\theta$ ratio gains $\omega_j^{1/2}\exp\{-\omega_j (\lambda^*_{jj})^2/(2 s_0^2)\}$. The default `0` keeps the flat prior of `meshed`. The option can be used with or without the horseshoe groups.

What has been established about this option is limited: on a case in which the flat prior behaved well, the proper prior gave the same answer, so it does no harm there. It has not been shown to cure the drift on the dataset where the drift was observed.

## 7. Non-Gaussian outcomes

For Poisson, binomial and beta outcomes the loadings are updated by the gradient-based samplers of Peruzzi and Dunson (2024), started at zero. They move slowly at first. If the scales are updated from the first iteration, they see $S_g \approx 0$, collapse to the lower end of the grid, and the prior precision then holds the loadings at zero: every cross-loading ended near $10^{-5}$ in our first Poisson run. The scales are therefore held at one for the first `n_burnin / 2` iterations when any outcome is non-Gaussian (`ghs_warmup`). For Gaussian outcomes the update of the loadings is an exact draw and no warm-up is used. The non-Gaussian path has been tested on one dataset only.

## 8. Posterior summaries of borrowing

**Borrowing share.** Because the factors are independent with unit variance a priori, the latent process of outcome $j$ has variance $\sum_{h \le j}\lambda_{jh}^2$, of which $\sum_{i \in B_r, i<j}\lambda_{ji}^2$ is carried by the factors anchored in block $r$:

$$\omega_{j,r} = \frac{\sum_{i \in B_r,\, i<j} \lambda_{ji}^2}{\sum_{h \le j} \lambda_{jh}^2} \in [0, 1].$$

`borrowing_share()` evaluates it draw by draw and reports quantiles and $\mathrm{pr}(\omega_{j,r} > \epsilon \mid \text{data})$, with $\epsilon = 0.05$ by default. It is an effect size, free of $n_j$ and $\tau_j^2$. The denominator contains the diagonal loading, so a share inherits any weakness in the identification of the marginal variance.

**Why not shrinkage weights.** The usual horseshoe summary is $\kappa_g = \{1 + n_j\zeta^2\gamma_g^2/\tau_j^2\}^{-1}$ thresholded at one half. It is a detection rule. With tens of thousands of locations and small measurement error, $n_j/\tau_j^2$ is so large that any loading that is not exactly zero is declared a signal: in one analysis the rule gave probability above 0.95 to borrowing for groups whose loadings were below 0.06 in absolute value on unit-variance outcomes. `shrinkage_weights()` is kept as a diagnostic.

**Decoupled selection.** Following Hahn and Carvalho (2015), and Bolfarine et al. (2024) for factor loadings, selection is a decision made after the posterior is obtained. For a sparse candidate $g_j$ for row $j$, the loss

$$E\big(\|\lambda_{[j,:]} - g_j\|^2 \mid \text{data}\big) + \mu \sum_{g \in \text{row } j} d_g^{1/2} \|g_{j,g}\|$$

is minimized by group soft-thresholding of the posterior mean, $g_{j,g} = (1 - \mu d_g^{1/2}/\{2\|\bar\lambda_{j,g}\|\})_+ \bar\lambda_{j,g}$, so the path drops groups in increasing order of $\|\bar\lambda_{j,g}\|/d_g^{1/2}$ (`path_score`). For a retained set $S$, the share of latent variance given up is $\delta_j(S) = \sum_{g \notin S}\omega_{j,g}$. `borrowing_map()` reports the sparsest $S$ on the path with $\mathrm{pr}\{\delta_j(S) \le \epsilon \mid \text{data}\} \ge$ `level`.

Both summaries are functions of the draws of $\Lambda$ and can be computed for an unrestricted fit with `blocks` supplied, which is how priors are compared on the same footing.

## 9. Source map

New files and the changes to files inherited from `meshed` 0.2.3.

| File | What it contains |
|---|---|
| `R/ghs.R` | `design_groups()`, `spmeshed_ghs()` (infers blocks from the missingness pattern, sets `k = q`, passes the group matrix), `zero_distance_cor()`, `shrinkage_weights()`, `print` method. New. |
| `R/borrowing.R` | `borrowing_share()`, `borrowing_map()`. New. |
| `R/spmeshed.r` | Reads `prior$lambda_group` (may mark only strictly lower entries inside the loading mask; otherwise an error) and `prior$lambda_diag_sd`, and passes them to the C++ entry point. Returns `ghs_lam_mcmc`, `ghs_tau_mcmc`. |
| `src/meshed.h` | Members `use_ghs`, `Lambda_group`, `n_ghs_groups`, `ghs_lam` ($\gamma$), `ghs_tau` ($\zeta$), `ghs_grid`, `ghs_warmup`, `ghs_iter`, `lambda_diag_sd`; methods `init_ghs()`, `ghs_sigmasq()`, `sample_ghs_scales()`, `ghs_theta_prior_logratio()`. |
| `src/meshed_betalambdatau_sampler.cpp` | Prior precision (ii) for grouped entries and for the diagonal, in both the Gaussian and the non-Gaussian update of the loadings; call to `sample_ghs_scales()` after the loadings; definitions of the four methods and of `ghs_grid_draw()`. |
| `src/meshed_theta_sampler.cpp` | Adds `ghs_theta_prior_logratio(theta_new, theta_old)` to the log prior ratio of the Metropolis step: factor (iii), including the diagonal term. Returns 0 when parameter expansion is off. |
| `src/meshed_mcmc.cpp` | Arguments `lambda_group` and `lambda_diag_sd`; initialization; warm-up length for non-Gaussian outcomes; storage and return of the scale draws. |
| `src/Makevars`, `src/Makevars.win` | `CXX_STD = CXX17`, required by current RcppArmadillo. |

In the C++ code the names follow the horseshoe literature (`ghs_lam` for local, `ghs_tau` for global scales); in the documentation they are $\gamma_g$ and $\zeta$.

## 10. Validation

Checks carried out during development. The scripts live with the manuscript materials and will be released with it.

| Check | Result |
|---|---|
| Prior-only sampler in the expanded coordinates, one loading, $\zeta = 1$, $\omega \sim$ Inv-Gamma(2, 1), 400,000 iterations | With factor (iii): median, 0.75 and 0.9 quantiles of $\lvert\lambda\rvert$ within 0.8% of the horseshoe, and $\omega$ keeps its prior (0.5%). Without it: quantiles displaced by 5%. |
| Grid full conditional of $\log\gamma_g$ against importance sampling from the prior (4 million draws) | 0.1, 0.5, 0.9 quantiles agree to 0.006 on the log scale. |
| Group soft-thresholding against a numerical optimizer of the decoupled-selection loss | Same minimum to $10^{-9}$. |
| Recovery of a known sparse loading matrix in 2D, two panels on disjoint sites | Non-null shares recovered, null share near zero (the example in the README). |
| 2D versus 3D working variance | Development version 0.2.3.9000 used $\sigma^2$ in both cases, which mis-scaled the prior for 2D data only. Corrected by derivation; all uses now go through `ghs_sigmasq()`, which calls the same `ps_back()` that `meshed` uses to recover the identified loadings. Runs made with the faulty version were discarded. |
| No groups, no diagonal prior | `use_ghs` is false, the prior precision is that of `meshed`, and `ghs_theta_prior_logratio()` returns 0, so the code path is that of `meshed` 0.2.3. |

## 11. Known issues

- Convergence of the decay of long-range factors can be slow. Trace plots of `theta_mcmc` and of the implied marginal variances `rowSums(lambda^2)` should be inspected; 5,000 iterations were not always enough in our simulations.
- The depth-decay prior defaults to the bounds of `phi`. Set `prior$a`.
- The initial value of the depth decay is the midpoint of its prior support, which can be far from the posterior when the support is wide; allow for it in the burn-in.
- Irregular 3D locations are slow (no geometry caching). About 6 seconds per iteration at 40,000 irregular locations and six outcomes, against 0.6 for 45,000 gridded voxels.
- Vignettes of `meshed` were removed; they have not been revised for this package.

## References

See [README.md](README.md#references), and

- Zhang, H. (2004). Inconsistent estimation and asymptotically equal interpolations in model-based geostatistics. *Journal of the American Statistical Association*, 99, 250-261.
- He, V. X. and Wand, M. P. (2025). The grouped horseshoe distribution and its statistical properties. *Stat*. doi:10.1002/sta4.70050
