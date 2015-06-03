\name{summary.etasclass}
\alias{summary.etasclass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary method for etasclass objects}
\description{This is the main method  to  summarize the output of an object of class \code{etasclass}.

Gives some information on the execution and gives estimates of the ETAS parameters together with the standard errors.
}
\usage{
\method{summary}{etasclass}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an \code{etaclass} object to pass to \code{summary}.}
  \item{...}{other arguments.}
}
\details{
Displays summary information about an object of class \code{etasclass}.
}

\value{
Displays AIC values, parameters estimates and their standard errors, together with some information on the execution of the \code{etasclass}
estimation process.  Displays also  the exact call of the function that generated \code{etasclass}
}



\seealso{ \code{\link{etasclass}},\code{\link{eqcat}}, \code{\link{profile.etasclass}}}


\author{
Marcello Chiodi, Giada Adelfio}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\dontrun{
data("italycatalog")
# load a sample catalog of the italian seismicity


etas.flp=etasclass(italycatalog,  
magn.threshold = 3.0,  magn.threshold.back = 3.5,
k0 = 0.005, c = 0.005, p = 1.01, a = 1.05, gamma = 0.6, q = 1.52, d = 1.1,
params.ind = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
declustering = TRUE, thinning = FALSE, flp = TRUE, ndeclust = 15,
onlytime = FALSE, is.backconstant = FALSE,
description = "etas flp",sectoday = TRUE, usenlm = TRUE, epsmax = 0.001)
# execution of etasclass for events with minimum magnitude of 3.0. 
# The events with magnitude at least 3.5 are used to build a first approximation
# for the background intensity function
# (magn.threshold.back=3.5)

# summary method for the etasclass object

summary(etas.flp)
Call: 
 
etasclass(cat.orig = italycatalog, magn.threshold = 3, magn.threshold.back = 3.5, 
    k0 = 0.005, c = 0.005, p = 1.01, a = 1.05, gamma = 0.6, d = 1.1, 
    q = 1.52, params.ind = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
        TRUE, TRUE), declustering = TRUE, thinning = FALSE, flp = TRUE, 
    ndeclust = 15, onlytime = FALSE, is.backconstant = FALSE, 
    description = "etas flp", sectoday = TRUE, usenlm = TRUE, 
    epsmax = 0.001)

 
 
etas flp 
Execution started:                  2015-06-02 13:01:04 
Elapsed time of execution (hours)   0.2473813 
Number of observations             2158 
Magnitude threshold                3 
declustering                        TRUE 
Number of declustering iterations   4 
Kind of declustering                weighting 
flp                                 TRUE 
sequence of AIC values for each iteration 
49620.08 48458.86 48418.2 48415.17 
 
------------------------------------------------------- 
 
ETAS Parameters: 
            Estimates       std.err.
mu           0.355850       0.011294
k0           0.008373       0.002053
c            0.009404       0.001795
p            1.121630       0.016271
a            1.509371       0.064077
gamma        0.857945       0.084688
d            1.915139       0.306384
q            1.836391       0.067067
------------------------------------------------------- 

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ETAS}
\keyword{earthquake}
\keyword{summary}

% __ONLY ONE__ keyword per line
