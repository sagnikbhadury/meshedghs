\name{borrowing_share}
\alias{borrowing_share}
\alias{borrowing_map}
\title{Posterior inference on borrowing between blocks of outcomes}
\description{
\code{borrowing_share} computes the posterior distribution of the borrowing share, the fraction of the latent spatial variance of an outcome that is carried by the factors anchored in a block of other outcomes.
\code{borrowing_map} reports, for each outcome, the sparsest set of blocks without which its latent variation cannot be reproduced to within a tolerance, by decoupled shrinkage and selection.
}
\usage{
borrowing_share(fit, blocks = NULL, grouping = "design", eps = 0.05,
                probs = c(0.025, 0.5, 0.975))

borrowing_map(fit, blocks = NULL, grouping = "design", eps = 0.05, level = 0.9)
}
\arguments{
\item{fit}{output of \code{\link{spmeshed_ghs}}, or of \code{\link{spmeshed}} with \code{k = q} factors, in which case \code{blocks} is required. The summaries are functions of the posterior draws of the loadings and are defined under any prior.}
\item{blocks}{integer vector with one entry per outcome; outcomes that share an observation pattern share a label. Must be nondecreasing. Taken from \code{fit} when it was produced by \code{spmeshed_ghs}.}
\item{grouping}{\code{"design"} for one group per (outcome, block) pair, \code{"element"} for one group per loading.}
\item{eps}{share below which borrowing is regarded as practically irrelevant.}
\item{probs}{posterior quantiles to report.}
\item{level}{posterior probability required by the selection rule.}
}
\details{
With \eqn{\Lambda} lower triangular and independent unit-variance factors, the latent process of outcome \eqn{j} has variance \eqn{\sum_h \lambda_{jh}^2}. The borrowing share of outcome \eqn{j} from block \eqn{r} is
\deqn{\omega_{j,r} = \sum_{i \in B_r, i < j} \lambda_{ji}^2 / \sum_{h \le j} \lambda_{jh}^2,}
evaluated draw by draw. It is an effect size: it does not depend on the number of observations or on the measurement error variance. The shares are directed, because outcome \eqn{j} can load only on the factors of outcomes that precede it.

\code{borrowing_map} follows Hahn and Carvalho (2015). Group soft-thresholding of the posterior mean loadings orders the groups of each outcome; groups are dropped in that order for as long as the variance share given up, \eqn{\delta_j(S) = \sum_{g \notin S} \omega_{j,g}}, satisfies \eqn{pr(\delta_j(S) \le eps | data) \ge level}.

Thresholding horseshoe shrinkage weights (\code{\link{shrinkage_weights}}) is not recommended for this purpose: with tens of thousands of locations and small measurement error, every loading that is not exactly zero is declared a signal.
}
\value{
\code{borrowing_share} returns a list with \code{table} (one row per group: outcome, the outcomes it borrows from, group size, posterior quantiles and \eqn{pr(\omega > eps)}), \code{draws} (posterior draws, one column per group) and \code{eps}.

\code{borrowing_map} returns a data frame with one row per group: outcome, group, the outcomes it borrows from, whether the group is retained, and its score on the group soft-thresholding path.
}
\references{
Hahn, P. R. and Carvalho, C. M. (2015) Decoupling shrinkage and selection in Bayesian linear models: a posterior summary perspective. \emph{Journal of the American Statistical Association}, 110, 435--448. \doi{10.1080/01621459.2014.993077}

Bolfarine, H., Carvalho, C. M., Lopes, H. F. and Murray, J. S. (2024) Decoupling shrinkage and selection in Gaussian linear factor analysis. \emph{Bayesian Analysis}, 19, 181--203. \doi{10.1214/22-BA1349}
}
\seealso{\code{\link{spmeshed_ghs}}, \code{\link{design_groups}}}
\examples{
\dontrun{
fit <- spmeshed_ghs(Y, x, coords, n_samples = 1000, n_burnin = 1500)
borrowing_share(fit)$table
borrowing_map(fit)
}
}
