\name{spmeshed_ghs}
\alias{spmeshed_ghs}
\alias{print.spmeshed_ghs}
\title{Meshed coregionalization model with design-based grouped horseshoe loadings}
\description{Fits a multivariate meshed Gaussian process regression with \eqn{k = q} latent factors and a grouped horseshoe prior on the strictly lower loadings, with groups given by blocks of outcomes that share an observation pattern.}
\usage{
spmeshed_ghs(y, x, coords, blocks = NULL, grouping = c("design", "element", "none"),
             prior = list(), ...)
}
\arguments{
\item{y}{matrix of outcomes with \eqn{n} rows and \eqn{q} columns; \code{NA} marks outcomes that are not observed at a location and are predicted within the sampler.}
\item{x}{matrix of covariates with \eqn{n} rows.}
\item{coords}{unnamed matrix of coordinates with 2 (space) or 3 (space and depth or time) columns.}
\item{blocks}{integer vector of length \eqn{q}, nondecreasing, giving the block (for example the assay panel) of each outcome. If \code{NULL}, blocks are inferred from the missingness patterns of the columns of \code{y}.}
\item{grouping}{\code{"design"}: one local scale for the loadings of each outcome on the factors of each block; \code{"element"}: one local scale per loading; \code{"none"}: the diffuse prior of \code{\link{spmeshed}}.}
\item{prior}{list passed to \code{\link{spmeshed}}; \code{prior$phi} is required, and \code{prior$a} should be set for three-dimensional coordinates.}
\item{...}{further arguments passed to \code{\link{spmeshed}}, such as \code{axis_partition}, \code{n_samples}, \code{n_burnin}, \code{n_threads}, \code{starting} and \code{settings}.}
}
\value{An object of class \code{c("spmeshed_ghs", "spmeshed")}: the output of \code{\link{spmeshed}} with the chains \code{ghs_lam_mcmc} (local scales) and \code{ghs_tau_mcmc} (global scale), and a list \code{ghs} holding the group matrix, the blocks and the number of observations of each outcome.}
\details{The horseshoe is placed on the identified loadings under the parameter expansion of the sampler, and the local and global scales are drawn exactly from their discretized full conditionals on a fixed grid. The gradient-based updates for non-Gaussian outcomes use the diffuse prior.}
\seealso{\code{\link{design_groups}}, \code{\link{shrinkage_weights}}, \code{\link{zero_distance_cor}}}
\author{Sagnik Bhadury}
