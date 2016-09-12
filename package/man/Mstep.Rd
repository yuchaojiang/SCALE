\name{Mstep}
\alias{Mstep}
\title{
Maximization step for the EM algorithm in gene classification.
}
\description{
In the gene classification step, we adopt an EM algorithm to classify each cell
into being silent (Off), monoallelically expressed ('A' or 'B'), or biallelically
expressed ('AB').
}
\usage{
Mstep(Zhat,cura,n,nA,nB)
}
\arguments{
    \item{Zhat}{
    Matrix specifying the posterior probability of assignment for
    each cell.}
    \item{cura}{
    Current parameter for a and b to be updated.}
    \item{n}{
    Number of cells.}
    \item{nA}{
    A-allele read counts.}
    \item{nB}{
    B-allele read counts.}
}
\value{
    \item{par.list }{List of updated parameters.}
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\seealso{
    \code{\link{Estep}},
    \code{\link{gene_classify}}
}
\examples{
# Mstep(Zhat,cura,n,nA,nB)
}
\keyword{ package }