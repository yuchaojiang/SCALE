\name{tech_bias}
\alias{tech_bias}
\title{
Adjusting for technical variability that is introduced by single-cell 
RNA-sequencing.
}
\description{
Dropouts, amplification and sequencing bias are adjusted for using spike-ins as
internal controls.
}
\usage{
tech_bias(spikein_input, alleleA, alleleB, pdf=NULL)
}
\arguments{
    \item{spikein_input}{
    Input for spike-ins. Including true number of molecules, lengths, and 
    observed read counts for spike-ins across cells.}
   \item{alleleA}{
    Matrix of A-allele expression, of dimension number of genes times number
    of cells.}
    \item{alleleB}{
    Matrix of B-allele expression, of dimension number of genes times number
    of cells.}
    \item{pdf}{
    Whether a pdf plot should be generated. TRUE by default.}
}
\value{
    \item{abkt }{Estimated paramters associated with technical variability.}
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\examples{
data(mouse.blastocyst)
alleleA=mouse.blastocyst$alleleA
alleleB=mouse.blastocyst$alleleB
spikein_input=mouse.blastocyst$spikein_input
abkt=tech_bias(spikein_input=spikein_input, alleleA = alleleA, 
               alleleB = alleleB, pdf=FALSE)
}
\keyword{ package }