\name{summary.etasclass}
\alias{summary.etasclass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary method for etasclass objects}
\description{This is the main method  to  summarize the output of an object of class \code{etasclass}.

It gives some information on the execution and gives estimates of the ETAS parameters together with the standard errors.

More detailed output is avaliable by inspecting \code{str(etasclass.object)}, and printing single objects.
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
# summary method for the etasclass object esecov1 and esecov5 

(see examples in \code{\link{etasclass}})

## only with one covariate, the magnitude, classical ETAS model
>   summary(esecov1)
Call: 
 
etasclass(cat.orig = catalog.withcov, magn.threshold = 2.5, magn.threshold.back = 3.9, 
    mu = 0.3, k0 = 0.02, c = 0.015, p = 0.99, gamma = 0, d = 1, 
    q = 1.5, params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, 
        TRUE), formula1 = "time ~  magn1.work- 1", declustering = TRUE, 
    thinning = FALSE, flp = TRUE, ndeclust = 15, onlytime = FALSE, 
    is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE, 
    compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36)

 
 
 
Execution started:                  2020-05-03 00:24:08 
Elapsed time of execution (hours)   0.2294818 
Number of observations             2226 
Magnitude threshold                2.5 
declustering                        TRUE 
Number of declustering iterations   6 
Kind of declustering                weighting 
flp                                 TRUE 
sequence of AIC values for each iteration 
44887.75 43348.46 43250.77 43249.77 43249.27 43249.19 
 
final AIC value 
44887.75 43348.46 43250.77 43249.77 43249.27 43249.19 
 
------------------------------------------------------- 
 
formula for covariates of the triggered components: 
time ~ magn1.work - 1
<environment: 0x55968d6fd660>
ETAS Parameters: 
                 Estimates       std.err.
mu                0.667509       0.022620
k0                0.022393       0.005781
c                 0.014769       0.002708
p                 1.110059       0.015709
gamma             0.000000       0.000000
d                 1.905461       0.260360
q                 1.947223       0.077627
magn1.work        0.740109       0.092558
------------------------------------------------------- 

#### using covariates
>   summary(esecov5)
Call: 
 
etasclass(cat.orig = catalog.withcov, magn.threshold = 2.5, magn.threshold.back = 3.9, 
    mu = 0.3, k0 = 0.02, c = 0.015, p = 0.99, gamma = 0, d = 1, 
    q = 1.5, params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, 
        TRUE), formula1 = "time ~ z + magn1.work +nstaloc_rev +min_distance_rev+distmin- 1", 
    declustering = TRUE, thinning = FALSE, flp = TRUE, ndeclust = 15, 
    onlytime = FALSE, is.backconstant = FALSE, sectoday = FALSE, 
    usenlm = TRUE, compsqm = TRUE, epsmax = 1e-04, iterlim = 100, 
    ntheta = 36)

 
 
 
Execution started:                  2020-05-03 12:22:31 
Elapsed time of execution (hours)   0.4827933 
Number of observations             2226 
Magnitude threshold                2.5 
declustering                        TRUE 
Number of declustering iterations   3 
Kind of declustering                weighting 
flp                                 TRUE 
sequence of AIC values for each iteration 
44693.04 42884.07 42706.16 
 
final AIC value 
44693.04 42884.07 42706.16 
 
------------------------------------------------------- 
 
formula for covariates of the triggered components: 
time ~ z + magn1.work + nstaloc_rev + min_distance_rev + distmin - 
    1
<environment: 0x55968d5ed118>
ETAS Parameters: 
                       Estimates       std.err.
mu                      0.705351       0.022740
k0                      0.073070       0.021194
c                       0.019396       0.003435
p                       1.154186       0.016874
gamma                   0.000000       0.000000
d                       1.942929       0.272434
q                       2.004915       0.084784
z                      -0.041256       0.005779
magn1.work              1.157698       0.085360
nstaloc_rev            -0.009010       0.001817
min_distance_rev       -0.011020       0.002804
distmin                -1.826717       0.298649
------------------------------------------------------- 


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ETAS}
\keyword{earthquake}
\keyword{summary}

% __ONLY ONE__ keyword per line
