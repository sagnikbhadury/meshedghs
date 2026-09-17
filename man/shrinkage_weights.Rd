\name{shrinkage_weights}
\alias{shrinkage_weights}
\alias{zero_distance_cor}
\title{Posterior summaries of borrowing}
\description{\code{shrinkage_weights} summarizes, for each group, the weight \eqn{\kappa_g = (1 + n_j \zeta^2\gamma_g^2/\tau_j^2)^{-1}}, which is close to one when outcome \eqn{j} does not borrow from the block. \code{zero_distance_cor} returns the posterior of the cross-correlations at zero distance implied by the loadings.}
\usage{
shrinkage_weights(fit, probs = c(0.025, 0.5, 0.975))
zero_distance_cor(fit, probs = c(0.025, 0.5, 0.975))
}
\arguments{
\item{fit}{output of \code{\link{spmeshed_ghs}} (or of \code{\link{spmeshed}} for \code{zero_distance_cor}).}
\item{probs}{posterior quantiles to report.}
}
\value{A data frame with one row per group, or a list with the posterior mean, quantiles and draws of the correlation matrix.}
\author{Sagnik Bhadury}
