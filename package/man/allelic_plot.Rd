\name{allelic_plot}
\alias{allelic_plot}
\title{
Plot adjusted allele-specific read counts for each gene with inferred kinetic
parameters.
}
\description{
Main plotting function for SCALE. For each gene, a pdf plot is generated with
adjusted allelic read counts as well as inferred kinetic parameters and 
hypothesis testing results.
}
\usage{
allelic_plot(alleleA,alleleB,gene.class.obj, allelic.kinetics.obj,
                      diff.allelic.obj,non.ind.obj,i)
}
\arguments{
    \item{alleleA}{
    Matrix of A-allele expression, of dimension number of genes times number
    of cells.}
    \item{alleleB}{
    Matrix of B-allele expression, of dimension number of genes times number
    of cells.}
    \item{gene.class.obj}{
    List returned from gene categorization.
    \code{\link{gene_classify}}}
    \item{allelic.kinetics.obj}{
    List returned from estimating allele-specific bursting kinetics.
    \code{\link{allelic_kinetics}}}
    \item{diff.allelic.obj}{
    List returned from testing two alleles sharing different bursting kinetics.
    \code{\link{diff_allelic_bursting}}}
    \item{non.ind.obj}{
    List returned from testing two alleles firing in a dependent fashion.
    \code{\link{non_ind_bursting}}}
    \item{i}{
    Index of gene to be plotted.}
}
\value{
Plot of gene-specific output by SCALE (pdf format).
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\examples{
# i=which(genename=='Btf3l4')
# allelic_plot(alleleA = alleleA, alleleB = alleleB,
#              gene.class.obj = gene.class.obj,
#              allelic.kinetics.obj = allelic.kinetics.obj,
#              diff.allelic.obj = diff.allelic.obj,
#              non.ind.obj = non.ind.obj, i= i)
}
\keyword{ package }