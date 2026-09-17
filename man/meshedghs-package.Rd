\name{meshedghs-package}

\alias{meshedghs}
\alias{meshedghs-package}

\docType{package}

\title{
Adaptive borrowing across misaligned outcomes in meshed Gaussian process coregionalization models
}

\description{
	\code{meshedghs} fits multivariate spatial regressions for misaligned outcomes in which a design-based grouped horseshoe prior on the off-diagonal loadings of a linear model of coregionalization lets each outcome borrow from a block of other outcomes only when the data support it; see \code{\link{spmeshed_ghs}}, \code{\link{design_groups}}, \code{\link{shrinkage_weights}} and \code{\link{zero_distance_cor}}. The meshed Gaussian process sampler is derived from the \code{meshed} package (version 0.2.3) by Michele Peruzzi, whose functions are retained: \code{meshed} is a flexible package for Bayesian regression analysis on spatial or spatiotemporal datasets. The main function for fitting regression models is \code{\link{spmeshed}}, which outputs posterior samples obtained from Markov chain Monte Carlo which can be summarised using standard tools. The package also provides a function \code{\link{rmeshedgp}} for quickly simulating correlated spatial or spatiotemporal data at a very large number of locations.
}

\details{
  The functions \code{\link{rmeshedgp}} and \code{\link{spmeshed}} are provided for prior and posterior sampling (respectively) of Bayesian spatial or spatiotemporal multivariate regression models based on Meshed Gaussian Processes as introduced by Peruzzi, Banerjee, and Finley (2020). Posterior sampling via \code{\link{spmeshed}} proceeds by default via GriPS as detailed in Peruzzi, Banerjee, Dunson, and Finley (2021). When at least one outcome is not modeled with Gaussian errors, sampling proceeds taking advantage of Metropolis-adjusted Langevin dynamics as detailed in Peruzzi and Dunson (2022).
}

\author{
Sagnik Bhadury. The meshed Gaussian process sampler is derived from the meshed package (version 0.2.3, GPL >= 3) by Michele Peruzzi, who was not involved in the development of this package.
}

\references{
	Peruzzi, M., Banerjee, S., and Finley, A.O. (2022)
  Highly Scalable Bayesian Geostatistical Modeling via Meshed Gaussian Processes on Partitioned Domains. \emph{Journal of the American Statistical Association}, 117(538):969-982. \doi{10.1080/01621459.2020.1833889}

  Peruzzi, M., Banerjee, S., Dunson, D.B., and Finley, A.O.  (2021)
  Grid-Parametrize-Split (GriPS) for Improved Scalable Inference in Spatial Big Data Analysis. \url{https://arxiv.org/abs/2101.03579}
  
  Peruzzi, M., Dunson, D.B. (2022)
  Spatial meshing for general Bayesian multivariate models. \url{https://arxiv.org/abs/2201.10080}
}

\keyword{ package }

\seealso{
	\code{\link{spmeshed}}, \code{\link{rmeshedgp}}
}
