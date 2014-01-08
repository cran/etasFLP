\name{summary.etasclass}
\alias{summary.etasclass}
\alias{print.etasclass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary method for etasclass objects}
\description{This is the main method  to  sumarize the output of an object of class \code{etasclass}.

Gives some information on the execution and gives estimates of the ETAS parameters together with the standard errors.
}
\usage{
\method{summary}{etasclass}(object,...)
\method{print}{etasclass}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an \code{etaclass} object to pass to \code{summary}.}
  \item{x}{an \code{etaclass} object to pass to \code{print}.}
  \item{...}{other arguments.}
}
\details{
Displays summary information about an object of class \code{etasclass}.
}

\value{
Display parameters estimates and their standard errors, together with some information on the execution of the \code{etasclass}
estimation process.
}


\seealso{ \code{\link{etasclass}},\code{\link{eqcat}}, \code{\link{profile.etasclass}}}


\author{
Marcello Chiodi, Giada Adelfio}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\dontrun{
data("italycatalog")
# load a sample catalog of the italian seismicity

class(italycatalog)<-"eqcat"
summary(italycatalog)
plot(italycatalog)

# apply summary and plot method to an earthquakes catalog

etas3.1<-etasclass(italycatalog,description="etas flp",magn.threshold=3.1,thinning=FALSE,flp=TRUE,
is.backconstant=FALSE,magn.threshold.back=3.5,sectoday=TRUE,
onlytime=FALSE,declustering=TRUE,epsmax=0.00001,
params.ind=c(1,1,1,1,1,1,1,1),k0=0.005,c=0.005,p=1.01,a=1.05,gamma=0.6,q=1.52,d=1.1,
compsqm=TRUE,usenlm=TRUE,ndeclust=15)

# execution of etasclass for events with minimum magnitude of 3.1. 
# The events with magnitude at least 3.5 are used to build a first approximation
# for the background intensity function
# (magn.threshold.back=3.5)

summary(etas3.1)
> summary(etas3.1)
[1] "etas flp"
Execution started:                  2013-12-17 12:15:51 
Elapsed time of execution (hours)      0.430857 
Number of observations              1700 
Magnitude threshold                    3.1 
Number of declustering iterations     10 
Kind of declustering                weighting 
[1] "sequence of AIC values for each iteration"
 [1] 40444.81 39058.33 39100.61 39101.30 39101.45 39027.12 39025.27 39025.10
 [9] 39025.07 39025.06
[1] "ETAS Parameters"
            Estimates       std.err.
mu           0.299199       0.010178
k0           0.008852       0.002835
c            0.012759       0.002754
p            1.149678       0.020219
a            1.641007       0.070026
gamma        0.932507       0.094307
d            2.010335       0.384201
q            1.927081       0.089653
> 

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ETAS}
\keyword{earthquake}
\keyword{summary}

% __ONLY ONE__ keyword per line
