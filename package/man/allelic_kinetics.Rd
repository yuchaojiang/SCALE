\name{allelic_kinetics}
\alias{allelic_kinetics}
\title{
Inferring allele-specific bursting kinetics by single-cell RNA sequencing with
adjustment of technical variability and complicating factors such as cell size.
}
\description{
The two alleles of a gene have two Poisson-Beta distributions with respective 
parameters. These two Poisson-Beta distributions share the same cell-size
factor. A hierarchical model is used to account for technical variability that
is introduced by sequencing and library prep. Histogram repiling method is used
to adjust for technical variability (bandwidth is optimized based on correlations
of the inferred kinetic parameters between the two alleles). Moment estimator is
used to estimate bursting kinetics.
}
\usage{
allelic_kinetics(alleleA,alleleB,abkt,gene.category,cellsize,pdf=NULL)
}
\arguments{
    \item{alleleA}{
    Matrix of A-allele expression, of dimension number of genes times number
    of cells.}
    \item{alleleB}{
    Matrix of B-allele expression, of dimension number of genes times number
    of cells.}
    \item{abkt}{
    Inferred parameters associated with technical variability.
    \code{\link{tech_bias}}}
    \item{gene.category}{
    Vector of gene categories from the Bayes framework of gene classification.
    \code{\link{gene_classify}}}
    \item{cellsize}{
    Vector of cell sizes.}
    \item{pdf}{
    Whether a pdf should be generated -- scatterplot of inferred bursting
    kinetics between the two alleles. TRUE by default.}
}
\value{
    \item{allelic.kinetics.obj }{List of inferred kinetic parameters.}
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\seealso{
    \code{\link{tech_bias}},
    \code{\link{diff_allelic_bursting}},
    \code{\link{non_ind_bursting}}
}
\examples{
data(mouse.blastocyst)
alleleA=mouse.blastocyst$alleleA
alleleB=mouse.blastocyst$alleleB
spikein_input=mouse.blastocyst$spikein_input
data(gene.class.obj)
gene.category=gene.class.obj$gene.category
abkt=tech_bias(spikein_input=spikein_input, alleleA = alleleA, 
               alleleB = alleleB, readlength = 50, pdf=FALSE)
cellsize=rep(1,ncol(alleleA))  # cell size input
allelic.kinetics.obj=allelic_kinetics(alleleA = alleleA[1:100,], 
                                      alleleB = alleleB[1:100,], 
                                      abkt = abkt, 
                                      gene.category = gene.category[1:100], 
                                      cellsize = cellsize, pdf=FALSE)
}
\keyword{ package }