\name{gene_classify}
\alias{gene_classify}
\title{
Gene classification based on allele-specific expression at single-cell level.
}
\description{
In the gene classification step, we adopt an EM algorithm to classify each cell
into being silent (Off), monoallelically expressed ('A' or 'B'), or biallelically
expressed ('AB').
}
\usage{
gene_classify(alleleA, alleleB, epsilon=NULL)
}
\arguments{
    \item{alleleA}{
    Matrix of A-allele expression, of dimension number of genes times number
    of cells.}
    \item{alleleB}{
    Matrix of B-allele expression, of dimension number of genes times number
    of cells.}
    \item{epsilon}{
    Sequencing errors, set at 0.001 by default.}
}
\value{
    \item{gene.class.obj }{List of gene categorization results.}
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\seealso{
    \code{\link{Estep}},
    \code{\link{Mstep}}
}
\examples{
data(mouse.blastocyst)
alleleA=mouse.blastocyst$alleleA
alleleB=mouse.blastocyst$alleleB
gene.class.obj=gene_classify(alleleA=alleleA[1:5,],alleleB=alleleB[1:5,])
}
\keyword{ package }