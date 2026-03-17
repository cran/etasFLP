\name{summary.etasclass}
\alias{summary.etasclass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary method for etasclass objects}
\description{This is the main method  to  summarize the output of an object of class \code{etasclass}.

It gives some information on the execution and gives estimates of the ETAS parameters together with the standard errors.

More detailed output is avaliable by inspecting \code{str(etasclass.object)}, and printing single objects.
}
\usage{
\method{summary}{etasclass}(object,full=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an \code{etaclass} object to pass to \code{summary}.}
  \item{full}{ logical. New in version 2.2. If \code{TRUE} the full original call will be displayed together with some more output.}
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
# summary method for the etasclass object esecov1 and esecov5 

(see examples in \code{\link{etasclass}})

## only with one covariate, the magnitude, classical ETAS model
>   summary(esecov1)

Execution started:                  2026-02-03 20:34:53 
Elapsed time of execution (hours)   0.072 
Number of observations             2226 
Magnitude threshold                2.5 
declustering                        TRUE 
Number of declustering iterations   7 
Kind of declustering                weighting 
flp                                 TRUE 
sequence of AIC values for each iteration 
43587.95 43252.54 43249.81 43249.27 43249.2 43249.19 43249.19 
 
------------------------------------------------------- 
 
formula for covariates of the triggered components: 
[1] "time ~  magnitude- 1"
ETAS Parameters: 
                Estimates       std.err.
mu               0.667562       0.022621
k0               0.022396       0.005782
c                0.014774       0.002708
p                1.110087       0.015702
gamma            0.000000       0.000000
d                1.905522       0.260380
q                1.947269       0.077635
magnitude        0.740080       0.092559
------------------------------------------------------- 

#### using covariates and option full=TRUE in summary
summary(esecov5,full=TRUE)
Call: 
 
etasclass(cat.orig = catalog.withcov, magn.threshold = 2.5, magn.threshold.back = 3.9, 
    mu = 0.3, k0 = 0.02, c = 0.015, p = 1.01, gamma = 0, d = 1, 
    q = 1.5, params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, 
        TRUE), formula1 = "time ~ z + magnitude +nstaloc_rev +min_distance_rev+distmin- 1", 
    declustering = TRUE, thinning = FALSE, flp = TRUE, ndeclust = 15, 
    onlytime = FALSE, is.backconstant = FALSE, sectoday = FALSE, 
    usenlm = TRUE, compsqm = TRUE, epsmax = 1e-04, iterlim = 100, 
    ntheta = 48)

 
 
Execution started:                  2026-02-03 20:39:13 
Elapsed time of execution (hours)   0.104 
Number of observations             2226 
Magnitude threshold                2.5 
declustering                        TRUE 
Number of declustering iterations   2 
Kind of declustering                weighting 
flp                                 TRUE 
sequence of AIC values for each iteration 
43117.52 42748.26 
 
------------------------------------------------------- 
 
formula for covariates of the triggered components: 
[1] "time ~ z + magnitude +nstaloc_rev +min_distance_rev+distmin- 1"
ETAS Parameters: 
                       Estimates       std.err.
mu                      0.698665       0.022730
k0                      0.068132       0.021121
c                       0.018731       0.003424
p                       1.147265       0.016833
gamma                   0.000000       0.000000
d                       1.900510       0.272133
q                       1.978678       0.084610
z                      -0.039953       0.005776
magnitude               1.150552       0.085389
nstaloc_rev            -0.008982       0.001816
min_distance_rev       -0.011203       0.002803
distmin                -1.831964       0.298728
------------------------------------------------------- 
[1] " summary of background seismicity weights "
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
2.200e-07 5.113e-05 8.611e-02 4.568e-01 9.970e-01 1.000e+00 

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ETAS}
\keyword{earthquake}
\keyword{summary}

% __ONLY ONE__ keyword per line
