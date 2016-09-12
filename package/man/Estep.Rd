\name{Estep}
\alias{Estep}
\title{
Expectation step for the EM algorithm in gene classification.
}
\description{
In the gene classification step, we adopt an EM algorithm to classify each cell
into being silent (Off), monoallelically expressed ('A' or 'B'), or biallelically
expressed ('AB').
}
\usage{
Estep(a, b, phi, epsilon, n, nA, nB)
}
\arguments{
    \item{a}{
    Hyperparameter for the beta-distribution with gene categorization.}
    \item{b}{
    Hyperparameter for the beta-distribution with gene categorization.}
    \item{phi}{
    Vector of proportions of each state.}
    \item{epsilon}{
    Sequencing errors, set at 0.001 by default.}
    \item{n}{
    Number of cells.}
    \item{nA}{
    A-allele read counts.}
    \item{nB}{
    B-allele read counts.}
}
\value{
    \item{Zhat }{Matrix specifying the posterior probability of assignment for
    each cell.}
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\seealso{
    \code{\link{Mstep}},
    \code{\link{gene_classify}}
}
\examples{
# Estep(a, b, phi, epsilon, n, nA, nB)
}
\keyword{ package }