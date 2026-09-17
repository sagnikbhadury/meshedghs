# meshedghs 0.1.3 (17 Sep 2026)

 * Optional proper prior on the diagonal loadings, `prior$lambda_diag_sd`, on the identified scale, with the matching term in the Metropolis ratio of the covariance parameters.

# meshedghs 0.1.2

 * Non-Gaussian outcomes: the horseshoe scales are held at one for the first half of the burn-in. Without this the scales collapsed onto the zero starting value of the gradient-based loading update.

# meshedghs 0.1.1

 * `borrowing_share()` and `borrowing_map()`: posterior inference on borrowing through an effect size and decoupled shrinkage and selection. `shrinkage_weights()` is now a diagnostic only.

# meshedghs 0.1.0

 * First version as a package of its own, derived from `meshed` 0.2.3 (M. Peruzzi, GPL >= 3).
 * Design-based grouped horseshoe prior on the strictly lower loadings of the LMC, specified on the identified loadings under parameter expansion; grid sampler for the local and global scales; `spmeshed_ghs()`, `design_groups()`, `zero_distance_cor()`.
 * Fixed during development: the working variance of the expansion is `sigmasq * phi^(-2 nu)` for the 2D Matern covariance, not `sigmasq`.
 * Builds with C++17. The vignettes of `meshed` were removed.

The entries below are the change log of `meshed`, kept for reference.

# 19 Sep 22 - `meshed` version 0.2.3

 * Fixed a bug in beta regression models for some gradient-based sampling methods
 * misc cleanup

# 28 Jan 22 - `meshed` version 0.2.1

 * Updated to reflect progress detailed in Peruzzi & Dunson (2022) <arXiv:2201.10080>
 * `spmeshed` can now use less memory by only storing essentials. 
 * Other changes and cleanup

#  6 Oct 21 - `meshed` version 0.1.4

 * version bump after fixing vignette bugs for CRAN

# 23 Sep 21 - `meshed` version 0.1.3

 * better memory management when exiting MCMC
 * fixed `M_PI` for upcoming Rcpp updates
 * other minor bugfixes

# 21 Jun 21 - `meshed` version 0.1.2

 * fixed a couple of typos in docs
 * more OMP options in prior sampling
 * fix compile error on Solaris

# 10 Jun 21 - `meshed` version 0.1.1

Addressing comments from CRAN:

 * fixed `par()` setting on vignette
 * set `verbose=0` in all vignettes and examples

# 9 Jun 21 - `meshed` version 0.1.0

First version submitted to CRAN.