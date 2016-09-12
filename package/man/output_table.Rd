\name{output_table}
\alias{output_table}
\title{
Generate output for SCALE.
}
\description{
A tab-delimited text file will be generated as final output of SCALE, where each
row corresponds to a gene.
}
\usage{
output_table(alleleA, alleleB, gene.class.obj, allelic.kinetics.obj,
             diff.allelic.obj, non.ind.obj)
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
}
\value{
    \item{SCALE.output }{Matrix output of SCALE, where rows correspond to genes,
    and columns correspond to inferred parameters and summary statistics.}
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\examples{
data(mouse.blastocyst)
alleleA=mouse.blastocyst$alleleA
alleleB=mouse.blastocyst$alleleB
data(gene.class.obj)
gene.category=gene.class.obj$gene.category
results.list=gene.class.obj$results.list
data(allelic.kinetics.obj)
data(diff.allelic.obj)
non.ind.obj=non_ind_bursting(alleleA = alleleA, alleleB = alleleB,
                             gene.category = gene.category,
                             results.list = results.list)
SCALE.output=output_table(alleleA=alleleA, alleleB=alleleB,
                          gene.class.obj = gene.class.obj,
                          allelic.kinetics.obj = allelic.kinetics.obj,
                          diff.allelic.obj = diff.allelic.obj,
                          non.ind.obj = non.ind.obj)
}
\keyword{ package }