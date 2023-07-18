\name{etasFLP-package}
\alias{etasFLP-package}
\alias{etasFLP}
\docType{package}
\title{Mixed FLP and ML Estimation of ETAS Space-Time Point Processes
}
\description{
\emph{New  version 2.2.0}. Covariates have been introduced to explain the effects of external factors on the induced seismicity. Since the parametrization is changed, the etasclass object created with the previous versions are not compatible with the one obtained with the current version. Estimation of the components of an ETAS (Epidemic Type Aftershock Sequence) model for earthquake description. Non-parametric background seismicity can be estimated through FLP (Forward Likelihood Predictive), while parametric components are estimated through maximum likelihood. The two estimation steps are alternated until convergence is obtained. For each event the probability of being a background event is estimated and used as a weight for declustering steps. Many options to control the estimation process are present, together with some diagnostic tools.
Some descriptive functions for earthquakes catalogs are  included; also \code{plot}, \code{print}, \code{summary}, \code{profile} methods are defined for main output (objects of class \code{etasclass}); update methods are now present.}

\details{
\tabular{ll}{
Package: \tab etasFLP\cr
Type: \tab Package\cr
Version: \tab 2.2.1\cr  
Date: \tab 2023-07-17\cr
License: \tab GPL (>=2) \cr
Imports: fields,maps
Depends: \tab R (>= 3.5.0), mapdata \cr
Suggests: \tab MASS \cr
}
\code{etasclass} is the main function of the package \code{etasFLP}: strongly renewed in version 2.0.
\code{update.etasclass} and \code{timeupdate.etasclass} two new different kind of updating existing \code{etasclass} objects. Very useful for large catalogues 
}
\author{
Marcello Chiodi and Giada Adelfio

Maintainer: Marcello Chiodi<marcello.chiodi@unipa.it>
}

\note{The package is intended for the estimation of the ETAS model for seismicity description (introduced by Ogata (1988), see reference),
but theoretically it can be used for other fields of application.
}

\references{


Adelfio G., Chiodi, M. (2009). Second-Order Diagnostics for Space-Time Point Processes with
Application to Seismic Events. \emph{Environmetrics}, \bold{20}(8), 895-911. doi:10.1002/env.961.

Adelfio, G.,  Chiodi, M. (2013) Mixed estimation technique in semi-parametric space-time point processes for earthquake description. 
\emph{Proceedings of the 28th International Workshop on Statistical Modelling 8-13 July, 2013, Palermo} (Muggeo VMR, Capursi V, Boscaino G, Lovison G, editors). Vol. \bold{1}. pp.65-70.

Adelfio, G., Chiodi, M. (2015) Alternated estimation in semi-parametric space-time branching-type point processes with application to seismic catalogs. 
\emph{Stochastic Environmental Research and Risk Assessment}    \bold{29}(2), pp. 443-450. DOI: 10.1007/s00477-014-0873-8


Adelfio G., Chiodi, M. (2015). FLP Estimation of Semi-Parametric Models for Space-Time
Point Processes and Diagnostic Tools. \emph{Spatial Statistics}, \bold{14}(B), 119-132. doi:10.1016/j.spasta.2015.06.004.

Adelfio G., Chiodi, M. (2020). Including covariates in a space-time point process with application to seismicity. \emph{Statistical Methods and Applications
},  doi:10.1007/s10260-020-00543-5.

Adelfio G., Schoenberg, FP (2009). Point Process Diagnostics Based on Weighted Second-
Order Statistics and Their Asymptotic Properties. \emph{The Annals of the Institute of Statistical
Mathematics}, \bold{61}(4), 929-948. doi:10.1007/s10463-008-0177-1.

Chiodi, M., Adelfio, G., (2011) Forward Likelihood-based predictive approach for space-time processes. \emph{Environmetrics}, vol. \bold{22} (6), pp. 749-757. DOI:10.1002/env.1121.


Chiodi, M., Adelfio, G., (2017) Mixed Non-Parametric and Parametric Estimation Techniques in {R} Package {etasFLP} for Earthquakes' Description. \emph{Journal of Statistical Software}, vol. \bold{76} (3), pp. 1-29.
 DOI: 10.18637/jss.v076.i03 

Console, R.,  Jackson,   D. D. and  Kagan, Y. Y.
Using the ETAS model for Catalog Declustering and Seismic
Background Assessment. \emph{Pure Applied Geophysics}
\bold{167}, 819--830 (2010). DOI:10.1007/s00024-010-0065-5.

Nicolis, O., Chiodi, M. and Adelfio G. (2015) Windowed ETAS models with application to the Chilean seismic catalogs, \emph{Spatial Statistics}, Volume \bold{14}, Part B, November 2015, Pages 151-165, ISSN 2211-6753, http://dx.doi.org/10.1016/j.spasta.2015.05.006.


Ogata, Y.   Statistical models for earthquake
occurrences and residual analysis for point processes.
\emph{Journal of the American Statistical Association},
\bold{83},  9--27 (1988).



Veen, A. , Schoenberg, F.P.  Estimation of space-time
branching process models in seismology using an EM-type algorithm.
\emph{Journal of the American Statistical Association}, \bold{103}(482), 614--624 (2008).


Zhuang, J., Ogata, Y.  and Vere-Jones, D.
Stochastic declustering of space-time earthquake occurrences.
\emph{Journal of the American Statistical Association},
\bold{97},  369--379 (2002). DOI:10.1198/016214502760046925.



}

\keyword{package}
\keyword{earthquakes}
\keyword{ETAS }
\keyword{PointProcess}
\keyword{MLE} 
\keyword{FLP} 
