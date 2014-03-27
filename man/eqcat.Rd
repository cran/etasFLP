\name{eqcat}
\alias{eqcat}
\alias{summary.eqcat}
\alias{plot.eqcat}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check earthquake catalog}
\description{
Preliminary check of the names of an earthquake catalog. \code{summary} and \code{plot} methods for earthquake catalogs.}
\usage{
eqcat(x)
\method{plot}{eqcat}(x,...)
\method{summary}{eqcat}(object,extended=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an earthquake catalog.}
  \item{object}{an \code{eqcat} object.}
  \item{extended}{if \code{TRUE} some extra summary functions are computed.}
    \item{...}{other arguments.}
}
\details{
Minimal check of an earthquake catalog; checks only if it is suitable for the use as argument of the functions of \code{etasFLP} (mainly \code{etasclass}); checks only the presence of
variables with the names \code{time},  \code{lat},   \code{long},  \code{z},     \code{magn1}.
\code{summary} and \code{plot} methods are defined for earthquake catalogs.
}
\value{
If the catalog passes the check, then the catalog is returned with the new class name \code{eqcat}; otherwise an error message is printed.
\item{cat}{a catalog is returned. If the check is ok, this is an \code{eqcat} class object.}
\item{ok}{A flag: \code{TRUE} if the check is ok; \code{FALSE} elsewhere.}
}
\seealso{ \code{\link{etasclass}}}
\author{
Marcello Chiodi}
\note{
In this first version if you have a catalog without the depth (\code{z}), please insert however a constant column. The depth can be used only in some plot and not in the estimation routines of the package etasFLP; \code{etasclass} uses only
\code{time},  \code{lat},   \code{long},  \code{magn1}.}

%% ~Make other sections like Warning with \section{Warning }{....} ~
%\references{}

%\seealso{}
\examples{
data(italycatalog)
f=eqcat(italycatalog)
print(f$ok)

}
\keyword{ catalog }
\keyword{ earthquake }

