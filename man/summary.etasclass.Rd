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

etas.flp=etasclass(italycatalog,  magn.threshold = 3.1,  magn.threshold.back = 3.5,
k0 = 0.005,c = 0.005,p = 1.01, a = 1.05, gamma = 0.6, q = 1.52, d = 1.1,
params.ind = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
declustering = TRUE, thinning = FALSE, flp = TRUE, ndeclust = 15,
onlytime = FALSE, is.backconstant = FALSE,
description = "etas flp",sectoday = TRUE, usenlm = TRUE, epsmax = 10e-04)

# execution of etasclass for events with minimum magnitude of 3.1. 
# The events with magnitude at least 3.5 are used to build a first approximation
# for the background intensity function
# (magn.threshold.back=3.5)

summary(etas.flp)
# summary merhod for the etasclass object
>summary(etas.flp)

Call: 
 
etasclass(cat.orig = italycatalog, magn.threshold = 3.1, magn.threshold.back = 3.5, 
    k0 = 0.005, c = 0.005, p = 1.01, a = 1.05, gamma = 0.6, d = 1.1, 
    q = 1.52, params.ind = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
        TRUE, TRUE), declustering = TRUE, thinning = FALSE, flp = TRUE, 
    ndeclust = 15, onlytime = FALSE, is.backconstant = FALSE, 
    description = "etas flp", sectodacatalog-search.htmly = TRUE, usenlm = TRUE, 
    epsmax = 0.001)

 
etas flp 
Execution started:                  2014-02-04 13:17:29 
Elapsed time of execution (hours)   0.358737 
Number of observations             1700 
Magnitude threshold                3.1 
Number of declustering iterations   7 
Kind of declustering                weighting 
sequence of AIC values for each iteration 
40444.81 39058.33 39100.61 39101.3 39101.45 39027.12 39025.27 
 
------------------------------------------------------- 
 
ETAS Parameters: 
            Estimates       std.err.
mu           0.299141       0.010177
k0           0.008847       0.002832
c            0.012747       0.002752
p            1.149504       0.020206
a            1.640902       0.070027
gamma        0.932499       0.094271
d            2.010186       0.384154
q            1.926670       0.089593
------------------------------------------------------- 

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ETAS}
\keyword{earthquake}
\keyword{summary}

% __ONLY ONE__ keyword per line
