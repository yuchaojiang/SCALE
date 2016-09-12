\name{diff_allelic_bursting}
\alias{diff_allelic_bursting}
\title{
Testing whether the two alleles of a gene share different bursting kinetic
parameters.
}
\description{
Nonparametric bootstrap hypothesis test on whether the two alleles of a gene
share the same burst frequency and burst size. Two versions of bootstrap are
implemented: (i) bootstrap the raw allele-specific expression and adjust
technical variability for each bootstrap sample; (ii) bootstrap the corrected
allele-specific expression where the technical variability is adjusted from the
raw expression data. The first version is implemented in mode 'raw' while the
second version is implemented in 'corrected'. Both methods return very similar
results while the 'corrected' mode runs faster.}
\usage{
diff_allelic_bursting(alleleA, alleleB, cellsize, gene.category, abkt,
                      allelic.kinetics.obj, nboot=NULL, mode=NULL)
}
\arguments{
    \item{alleleA}{
    Matrix of A-allele expression, of dimension number of genes times number
    of cells.}
    \item{alleleB}{
    Matrix of B-allele expression, of dimension number of genes times number
    of cells.}
    \item{cellsize}{
    Vector of cell sizes.}
    \item{gene.category}{
    Vector of gene categories from the Bayes framework of gene classification.
    \code{\link{gene_classify}}}
    \item{abkt}{
    Inferred parameters associated with technical variability.
    \code{\link{tech_bias}}}
    \item{allelic.kinetics.obj}{
    List returned from estimating allele-specific bursting kinetics.
    \code{\link{allelic_kinetics}}}
    \item{nboot}{
    Number of bootstrap samples. 10,000 by default.}
    \item{mode}{
    Running mode. Can be either 'raw' or 'corrected'. 'corrected' by default.
    }
}
\value{
    \item{diff.allelic.obj }{List of p-values from testing of shared kinetics 
    between the two alleles.}
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\seealso{
    \code{\link{gene_classify}},
    \code{\link{tech_bias}},
    \code{\link{allelic_kinetics}}
}
\examples{
# diff.allelic.obj=diff_allelic_bursting(alleleA = alleleA,
#                                        alleleB = alleleB,
#                                        cellsize = cellsize,
#                                        gene.category = gene.category,
#                                        abkt = abkt,
#                                        allelic.kinetics.obj = allelic.kinetics.obj,
#                                        mode='corrected')
}
\keyword{ package }