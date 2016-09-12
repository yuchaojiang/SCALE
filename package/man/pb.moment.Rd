\name{pb.moment}
\alias{pb.moment}
\title{
Moment estimator for Poisson-Beta distribution with cell size factor.
}
\description{
To estimate the allele-specific kinetic parameters, we propose a moment estimator 
method that also integrates and adjusts for cell size factor.
}
\usage{
pb.moment(Y.temp,cellsize)
}
\arguments{
    \item{Y.temp}{
    Vector of allele-specific read counts adjusted for technical variability.}
    \item{cellsize}{
    Vector of cell sizes.}
}
\value{
    \item{kinetic.estimates }{Estimated activation, deactivation and transcription rates.}
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\seealso{
    \code{\link{allelic_kinetics}},
    \code{\link{diff_allelic_bursting}}
}
\examples{
# pb.moment(Y.temp,cellsize)
}
\keyword{ package }