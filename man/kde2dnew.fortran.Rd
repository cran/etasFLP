\name{kde2dnew.fortran}
\alias{kde2dnew.fortran}
\alias{kde2d.integral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A 2-d normal kernel estimator}
\description{A simple and quick 2-d weighted normal kernel estimator, with fixed bandwidth 
and relative integral.}

\usage{
kde2dnew.fortran(xkern,ykern,gx,gy,h,factor.xy=1,eps=0,w=replicate(length(xkern),1),
hvarx=replicate(length(xkern),1),hvary=replicate(length(xkern),1)
)

kde2d.integral(xkern,ykern,gx=xkern,gy=ykern,eps=0,factor.xy=1,
           h  = c(   bwd.nrd(xkern,w),bwd.nrd(ykern,w)),w=replicate(length(xkern),1),
           hvarx=replicate(length(xkern),1),hvary=replicate(length(xkern),1)
           )

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xkern}{x-values of kernel points of length \code{n} (\code{n}=\code{length(xkern)}).}
  \item{ykern}{y-values of kernel points of length \code{n}.}
  \item{gx}{x-values of the points where densities must be estimated.}
  \item{gy}{y-values of the points where densities must be estimated.}
  \item{h}{bandwidths: a length 2 numerical vector.}
  \item{eps}{enlargment factor for the region of interest.}
  \item{factor.xy}{expansion factor for bandwidths (density will be smoother if  \code{factor.xy}>1).}
  \item{w}{vector of weights to give to observed points (length \code{n}).}[]
    \item{hvarx}{Longitude bandwidths adjustement used in the kernel estimator of background seismicity. The length must be equal to the number of events of the catalog after event selection (can be less than \code{nrow(cat.orig)}). Default value = \code{replicate(length(xkern),1)} }
  \item{hvary}{Longitude bandwidths adjustement used in the kernel estimator of background seismicity. The length must be equal to the number of events of the catalog after event selection (can be less than \code{nrow(cat.orig)}).Default value = \code{replicate(length(xkern),1)}}
 }
\details{A standard bivariate normal kernel estimator.}
\value{grid values and estimated densities.
}
\references{
Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. 
Fourth edition. Springer.
Wand, M.P and Jones, M.C. (1995). \emph{Kernel Smoothing}. London: Chapman & Hall/CRC. 
}
\author{Marcello Chiodi.}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{kernel}
