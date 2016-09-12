\name{non_ind_bursting}
\alias{non_ind_bursting}
\title{
chi-square testing on whether the two alleles of a gene burst independently.
}
\description{
For each gene, a chi-square test of independence is carried out on whether the
two alleles of a gene fire independently. The number of cells being silent ('Off'),
monoallelically expressed ('A' or 'B'), or biallelically expressed ('AB') is
based on the posterior assignment from the Bayes categorization framework.
}
\usage{
non_ind_bursting(alleleA, alleleB, gene.category, results.list)
}
\arguments{
     \item{alleleA}{
    Matrix of A-allele expression, of dimension number of genes times number
    of cells.}
    \item{alleleB}{
    Matrix of B-allele expression, of dimension number of genes times number
    of cells.}
    \item{gene.category}{
    Vector of gene categories from the Bayes framework of gene classification.
    \code{\link{gene_classify}}}
    \item{results.list}{
    List of cell posterior assignments for all genes.
    \code{\link{gene_classify}}}
}
\value{
    \item{non.ind.obj }{List of p-values and type of dependent bursting.}
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
non.ind.obj=non_ind_bursting(alleleA = alleleA, alleleB = alleleB,
                             gene.category = gene.category,
                             results.list = results.list)
}
\keyword{ package }