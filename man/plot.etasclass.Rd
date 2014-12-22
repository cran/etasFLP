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
\method{plot}{etasclass}(x,pdf=FALSE,file ="etasplot", 
ngrid=201,nclass=10,tfixed=0,flag.3D=FALSE,flag.log=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \code{etaclass} object.}
%  \item{which}{an array with the desired graphs. Default is 1:4}
  \item{pdf}{If \code{TRUE}, then 2D plots are sent to a pdf file}
  \item{file}{name of the pdf file}
   \item{ngrid}{number of points for each direction (\code{x}, \code{y}) of a \code{ngrid*ngrid} grid where estimated  intensities must be evaluated. Default value= 201.}
   \item{nclass}{number of class for each direction (\code{x}, \code{y}) of a grid of \code{nclass*nclass} cells where estimated  intensities must be evaluated.
   Must divide \code{ngrid-1}. Default value= 10.}
  \item{tfixed}{If a positive value is given, then the triggered intensity at time \code{tfixed} is estimated and visualized.}
  \item{flag.3D}{If \code{TRUE} a 3D plot is also produced.}
  \item{flag.log}{If \code{TRUE} then a  \code{log} scale is used to plot intensities.}
  \item{...}{other arguments.}

}
\details{

Different plots of the output of an object of class \code{etasclass}.

By default the space-time region is the same used for the estimation of the ETAS model. 
Background, triggered and total space intensities
are also computed and plotted for a grid of values.

If a positive value is given for \code{tfixed}, then the triggered intensity at time \code{tfixed} is estimated and visualized.
A tipical use can be with \code{tfixed} a day after a big earthquake.

Starting with the package version 1.2.0 different kind of residual analysis are computed and visualized, separately for the space and time dimensions. (8 plot on three windows for the space and 2 plots on one window for the time)

For space dimension, 




Space residuals are computed dividing the observed rectangular space area in
a equally spaced  grid of \code{nclass} intervals for each dimension, so to
divide the observed space area in \code{nclass} x \code{nclass} rectangular cells.
We obtain the classical comparison between observed and
theoretical frequencies. All frequencies are related to the whole time interval (and thus theoretical frequencies are obtained integrating estimated intensities with respect to time).

Fifth graph (image plot)

We define  \code{nclass} x \code{nclass} standardized
residuals:
\deqn{
z_{\ell j} \ = \frac{n_{\ell j}-\hat{\nu}_{\ell j}}{\sqrt{\hat{\nu}_{\ell j}}}
\qquad (\ell =1,2,...,nclass; \ j=1,2,...,nclass)
 }
For each cell \eqn{\ell j}{lj} we have observed (\eqn{n_{\ell j}}) and
theoretical frequency (\eqn{\hat{\nu}_{\ell j}}).

Sixth graph (image plot)

We used a similar technique  to compute residuals for the
background seismicity only, to check if at least the estimation of
the  background component is appropriate. To this purpose the
observed background frequencies (\eqn{{}_b n_{ \ell j}}) are now
computed  by the sum of the estimated weights \code{rho.weights} and the theoretical background frequency \eqn{{}_b \hat{\nu}_{\ell j}} by the  estimated
marginal space background intensity in each cell.

From these quantities we obtain
\code{nclass} x \code{nclass} standardized residuals for the background intensity
only:

\deqn{
{}_b z_{\ell j} \ = \frac{{}_b n_{\ell j} \ - \ {}_b
\hat{\nu}_{\ell j}}{\sqrt{{}_b \hat{\nu}_{\ell j}}} \qquad (\ell =1,2,...,nclass; \
j=1,2,...,nclass)
}

seventh plot: (space intensities (integrated over time))

A 3x2 plot: first column for observed vs.theoretical  , second column  for standardized residuals vs theoretical values. First row for total intensity, second row for background intensity, and third row for their difference, the triggered intensities


eight-th graph:

To check departure of the model for the time dimension, we first
integrated the estimated intensity function with respect to the
observed space region, so to obtain an estimated time process (a
 one dimensional ETAS model):

 \deqn{\hat{\lambda}(t)=
 \int \int_{\Omega_{(x,y)}}\,
\hat{\lambda}(x,y,t)\,d x \, d y }

 As known, a non-homogeneous time
process can be transformed to a homogeneous one through the
integral transformation:


\deqn{\tau_i = 
\int_{t_0}^{t_{i}} 
\hat{\lambda}(t) \ d t
}

Then, a plot of \eqn{\tau_i} versus
\eqn{i} can give information about the departures of the models in the
time dimension. In particular, this plot, together with a plot of
the estimated time intensities, drawn on the same graphic winodw, can inform  on the time at which
departures are more evident



If \code{pdf=TRUE} all graphs are printed on a pdf file, as spcified by \code{file}; otherwise 
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
