\name{bwd.nrd}
\alias{bwd.nrd}
\title{
 Silverman's rule optimal for the estimation of a kernel bandwidth}
\description{Computes the optimal bandwidth with the Silverman's rule of thumb,
to be used for a kernel estimator with given points and weights.}


\usage{bwd.nrd(x, w, d = 2)}
\arguments{
  \item{x}{numeric vector: sample points to be used for a normal kernel estimator.}
  \item{w}{numeric vector of the same length of \code{x}: weights to give to the elements of \code{x}.}
  \item{d}{number of dimensions of the kernel estimator.}
}
\details{Computes the optimal bandwidth with the Silverman rule, for a kernel
estimator with points \code{x} and weights \code{w}. 
If a multivariate kernel is used, (i.e. \code{d} > 1),   
\code{bwd.nrd} must be called for each variable. It computes dispersion only with the weighted standard deviation, with no robust alternative. Called by \code{kde2dnew.fortran}.}
\value{The value of the bandwidth for a sample \code{x} and weights \code{w}.
}
\note{It is used in connection with the the declustering method of \code{etasFLP}. Points with an higher probability of being part of the background seismicity will weight more in the estimation the background seismicity.}
\references{
Silverman, B.W. (1986). \emph{Density Estimation for Statistics and Data Analysis}. Chapman and Hall: London.}
\author{
Marcello Chiodi}
\note{This is a slight  modification of \code{\link{bw.nrd}}.}
%\seealso{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
## The function is currently defined as
function (x,w,d=2) 
{
    if (length(x) < 2L) 
        stop("need at least 2 data points")
     m<-weighted.mean(x,w) 
     return(sqrt(weighted.mean((x-m)^2,w)) * 
	      (length(x)*(d+2)/4)^(-1/(d+4)))
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{bandwidth}
\keyword{kernel}