\name{design_groups}
\alias{design_groups}
\title{Groups of loadings implied by a measurement design}
\description{Builds the \eqn{q \times q} integer matrix that assigns each strictly lower loading of a lower-triangular coregionalization matrix to a local horseshoe scale.}
\usage{design_groups(blocks, grouping = c("design", "element", "none"))}
\arguments{
\item{blocks}{nondecreasing integer vector of length \eqn{q} with the block of each outcome.}
\item{grouping}{\code{"design"}, \code{"element"} or \code{"none"}; see \code{\link{spmeshed_ghs}}.}
}
\value{An integer matrix with 0 for entries that are not shrunk and \eqn{g \ge 1} for the entries of group \eqn{g}.}
\examples{
design_groups(c(1, 1, 1, 2, 2, 2))
design_groups(c(1, 1, 2, 3, 3, 3), "element")
}
\author{Sagnik Bhadury}
