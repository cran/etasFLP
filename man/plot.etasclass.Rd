\name{plot.etasclass}
\alias{plot.etasclass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot method for etasclass objects}
\description{This is the main method  to  visualize graphically the output of an object of class \code{etasclass}.

By default the space-time region is the same used for the estimation of the ETAS model. 
Background, triggered and total space intensities
are also plotted for a grid of values.
}
\usage{
\method{plot}{etasclass}(x,pdf=FALSE,file ="etasplot", ngrid=201,flag.3D=TRUE,flag.log=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \code{etaclass} object.}
  \item{pdf}{If \code{TRUE}, then 2D plots are sent to a pdf file}
  \item{file}{name of the pdf file}
   \item{ngrid}{number of points for each direction (\code{x}, \code{y}) of a \code{ngrid*ngrid} grid where estimated  intensities must be evaluated. Default value= 201.}
  \item{flag.3D}{If \code{TRUE} a 3D plot is also produced.}
  \item{flag.log}{If \code{TRUE} then a  \code{log} scale is used to plot intensities.}
  \item{...}{other arguments.}

}
\details{

Different plots of the output of an object of class \code{etasclass}.

By default the space-time region is the same used for the estimation of the ETAS model. 
Background, triggered and total space intensities
are also computed and plotted for a grid of values.

If \code{pdf=TRUE} intensities are printed on a pdf file, as spcified by \code{file}; otherwise 
default screen device is used. 
}

\value{This 
\code{plot} method computes, among others, \code{back.grid}, \code{trig.grid}, 
with coordinates \code{x.grid} and \code{y.grid}  used to obtain image plots of background, triggered and total spatial estimated intensities
(see \code{\link{etasclass}} to see the details of the mixed estimation method used).

   \item{x.grid}{x grid values.}
  \item{y.grid}{y grid values.}
 \item{back.grid}{background intensity estimated on a \code{ngrid} x \code{ngrid} grid.}
  \item{trig.grid}{triggered intensities estimated on a grid of \code{ngrid} x \code{ngrid} points.}
  \item{tot.grid}{total intensities estimated on a grid of \code{ngrid} x \code{ngrid} points.}
 \item{back.grid}{background space intensity estimated for observed points.}
  \item{trig.grid}{triggered space intensities estimated for observed points.}
  \item{tot.grid}{total space intensities estimated for observed points.}

}

\note{In this first version the x-y space region, where the point process is defined, by default is  a rectangle embedding the catalog values.
}

\seealso{\code{\link{etasclass}},  \code{\link{eqcat}}, \code{\link{profile.etasclass}}}

\author{
Marcello Chiodi, Giada Adelfio}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\dontrun{
data("italycatalog")
# load a sample catalog of the italian seismicity

class(italycatalog)<-"eqcat"

etas3.1<-etasclass(italycatalog,description="etas flp",magn.threshold=3.1,thinning=FALSE,flp=TRUE,
is.backconstant=FALSE,magn.threshold.back=3.5,sectoday=TRUE,
onlytime=FALSE,declustering=TRUE,epsmax=0.00001,
params.ind=c(1,1,1,1,1,1,1,1),k0=0.005,c=0.005,p=1.01,a=1.05,gamma=0.6,q=1.52,d=1.1,
compsqm=TRUE,usenlm=TRUE,ndeclust=15)

# execution of etasclass for events with minimum magnitude of 3.1. 
# The events with magnitude at least 3.5 are used to build a first approximation
# for the background intensity function
# (magn.threshold.back=3.5)


# plot method

> plot(etas3.1)

# profile likelihood for the 5-th parameter (a), with plot:

prof=profile(etas3.1,nprofile=7,iprofile=5)
plot(prof)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ETAS}
\keyword{plot}
\keyword{earthquake}

% __ONLY ONE__ keyword per line
