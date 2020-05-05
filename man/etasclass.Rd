\name{etasclass}
\alias{etasclass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mixed estimation of an ETAS model (renewed in version 2.0)}
\description{\code{etasclass} is the main function of the package \code{etasFLP}.

\code{etassclass} \emph{objects of previous versions
are not compatible with the  current version}

Performs the estimation of the components of the ETAS (Epidemic Type Aftershock Sequence) model for the description of the seismicity in a space-time region. Background seismicity is estimated non-parametrically, while triggered seismicity is estimated by MLE. In particular also the bandwidth for a kernel smoothing can be estimated through the Forward Likelihood Predictive (FLP) approach. For each event the probability of being a background event or a triggered one is estimated. 

\strong{New in version 2.0.0:}
\emph{Covariates have been introduced to explain the effects of external factors on the induced seismicity. Since the parametrization is changed, the etasclass object created with the previous versions are not compatible with the one obtained with the current version.}

An ETAS with up to 7+ncov parameters can be estimated, with several options and different methods.

Returns an \code{etasclass} object, for which \code{plot}, \code{summary}, \code{print} and \code{profile} methods are defined. 
}

%    parallel betacov formula1 offset hvarx hvary


\usage{  etasclass(cat.orig,
         magn.threshold	=2.5,
				 magn.threshold.back=magn.threshold+2,	
         tmax		=max(cat.orig$time),

##### starting values for parameters
				mu		=1,
				k0		=1,
				c		=0.5,
				p		=1.01,
##				a		=1.2,
				gamma		=.5,
				d		=1.,
				q		=1.5,
        betacov         =0.7,
### indicators
				params.ind=replicate(7,TRUE),
### formula for covariates:
                                formula1            ="time~magn1.work-1",
                                offset              =0,
### anisotropic smoothing is not  implemented
###  ellipsoidal spatial function in triggered intensity not yet implemented
				hdef=c(1,1),
				w		=replicate(nrow(cat.orig),1),
                                hvarx            =NULL,
                                hvary            =NULL,
### flags for the kind of declustering and smoothing:
                                declustering	=TRUE,
				thinning	=FALSE,
				flp		=TRUE,
				m1		=NULL,
				ndeclust        =5,
				onlytime	=FALSE,
				is.backconstant	=FALSE,
##### end of  main input arguments. 
##### Control and secondary arguments:
				description	="",
				cat.back   	=NULL,
				back.smooth	=1.0,
				sectoday	=TRUE,
				longlat.to.km   =TRUE,
				usenlm		=TRUE,
				method		="BFGS",
				compsqm 	=TRUE,
				epsmax		= 0.0001,
				iterlim		=50,      
				ntheta		=36) 
				}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cat.orig}{
An  earthquake catalog, possibly an object of class \code{eqcat}, or however a \code{data.frame} with variables of names \code{time},  \code{lat},   \code{long},  \code{z},     \code{magn1}. No missing values are allowed.
}
\item{magn.threshold}{
Threshold magnitude (only events with a magnitude at least \code{magn.threshold} will be used).  Default value = 2.5.
}
  \item{magn.threshold.back}{
Threshold magnitude used to build the  catalog \code{cat.back} for
the first estimation of the background seismicity. Default value = \code{magn.threshold+2}.
}
  \item{tmax}{
Maximum value of time. Only observation before \code{tmax} will be used for estimation. Default value = \code{max(cat.orig$time)}.
}

\emph{Values for the 7 parameters of the ETAS   model} (starting values or fixed values according to  \code{params.ind}):
  \item{mu}{Parameter 1 (\eqn{\mu}) of the ETAS model: background general intensity; see details. Default value = 1.
%%     ~~Describe \code{lambda} here~~
}
  \item{k0}{Parameter 2 (\eqn{\kappa_0}) of the ETAS model: measures the strength
of the aftershock activity; see details. Default value = 1.
%%     ~~Describe \code{k0} here~~
}
  \item{c}{Parameter 3 of the ETAS model; a shift parameter of the Omori law for temporal decay rate of aftershocks; see details. Default value = 0.5.
%%     ~~Describe \code{c} here~~
}
  \item{p}{Parameter 4  of the ETAS model; the exponent of the Omori law for temporal decay rate of aftershocks; see details. Default value = 1.01.
%%     ~~Describe \code{p} here~~
}
  \item{gamma}{Parameter 5 (\eqn{\gamma}) of the ETAS model; together with  \code{a} is related to the efficiency of an event of given magnitude 
  in generating aftershocks; see details. Default value = 0.5.
%%     ~~Describe \code{gamma} here~~
}
  \item{d}{Parameter 6 of the ETAS model; parameter related to the spatial influence of the
mainshock; see details. Default value = 1.
%%     ~~Describe \code{d} here~~
}
  \item{q}{Parameter 7 of the ETAS model; parameter related to the spatial influence of the
mainshock; see details. Default value = 1.5.
%%     ~~Describe \code{q} here~~
}
 \item{betacov}{Parameter of the covariates ETAS model; 
 see details. Default value = 0.7.
%%     ~~Describe \code{q} here~~
}

\emph{End of model pararameter input}
  \item{params.ind}{vector of 7 logical values: \code{params.ind[i] = TRUE} means that the i-th parameter must be estimated.  \code{params.ind[i] = FALSE}  means that the i-th parameter is fixed to its input value  (the order of parametrs is: \code{mu}, \code{k0}, \code{c}, \code{p}, \code{gamma}, \code{d}, \code{q}).
  Default value = \code{replicate(7,TRUE)}, that is, \code{etasclass} estimates all parameters.
}
\item{formula1}{a character variable: Formula which defines the covariates acting on the induced seismicity. In classical etas model   the covariate is the magnitude. The left side (dummy) element must be the time, which is a variable certainly present in the data set. The right part of the formula determines \code{ncov} the number of covariates. Default value=\code{"time~magn1.work-1"}; input \emph{must} be a character value: it is converted in a formula inside the program}
\item{offset}{An offset, for which no parameter will be estimated. Default value=0}

\emph{Flags for the kind of declustering and smoothing}:
%  \item{kern.var}{if \code{TRUE} the background seismicity is estimated through a variable metric anisotropic kernel. yet not used in version 1.1.0. Default value = \code{FALSE}.}
  \item{hdef}{Starting values for the \code{x,y} bandwidths used in the kernel estimator of background seismicity. Default value = \code{1,1}.
}
  \item{w}{Starting values for the weigths used in the kernel estimator of background seismicity. The length must be equal to the number of events of the catalog after event selection (can be less than \code{nrow(cat.orig)}). 
  
Default value = \code{replicate(nrow(cat.orig),1).}
}
  \item{hvarx}{Longitude bandwidths adjustement used in the kernel estimator of background seismicity. The length must be equal to the number of events of the catalog after event selection (can be less than \code{nrow(cat.orig)}). Default value = NULL}
  \item{hvary}{Longitude bandwidths adjustement used in the kernel estimator of background seismicity. The length must be equal to the number of events of the catalog after event selection (can be less than \code{nrow(cat.orig)}). Default value = NULL}
  \item{declustering}{if \code{TRUE} the catalog is iteratively declustered to optimally estimate the background intensity (through thinning, if \code{thinning=TRUE}, or through weighting if \code{thinning=FALSE}). Default value = \code{TRUE}.
}
  \item{thinning}{if \code{thinning=TRUE} a background catalog is obtained sampling from the original catalog with probabilities estimated during the iterations. Default value =\code{FALSE}.
}
  \item{flp}{if \code{flp=TRUE} then background seismicity is estimated through Forward Likelihood Predictive (see details). Otherwise the Silverman rule is used. Default value =\code{TRUE}.
} 
  \item{m1}{Used only if \code{flp=TRUE}. Indicates the range of points used for the FLP steps. See details. If missing it is set  to \code{nrow(cat)/2}.
} 
\item{ndeclust}{maximum number of iterations for the general declustering procedure. Default=5.
}
\item{onlytime}{if \code{TRUE}  then a time process is fitted to data , regardless to space location (in this case \code{is.backconstant} is set to \code{TRUE} and \code{declustering}, \code{flp} are set to \code{FALSE}). Default value = \code{FALSE}.
}
  \item{is.backconstant}{if \code{TRUE}  then background seismicity is assumed to be homogeneous in space  (and \code{declustering, flp} are set to \code{FALSE}).  Default value = \code{FALSE}.
}
\emph{Other control parameters}:
  \item{description}{a description string used for the output. Default value = "".
}
  \item{cat.back}{ external catalog used for the estimation of the background seismicity. 
 Default value = \code{NULL}.
}
\item{back.smooth}{
Controls the level of smoothing for the background seismicity (meaningful only if \code{flp=FALSE}). Default value = 1.
}
  \item{sectoday}{ if \code{TRUE}, then \code{time} variable of \code{cat.orig} is converted from seconds to days.  Default value = \code{TRUE}.
}
  \item{longlat.to.km}{ if \code{TRUE}, then \code{long} and \code{lat} variables of \code{cat.orig} are treated as geographical coordinates and converted to kilometers.  Default value = \code{TRUE}.
}
  \item{usenlm}{if \code{TRUE}, then \code{nlm} function (gauss-newton method) is used in the maximum likelihood steps; if \code{FALSE}, then \code{optim} function is used (with \code{method} \code{=method} ).  Default value = \code{TRUE}.
}
  \item{method}{used if \code{usenlm=FALSE}: method used  by \code{optim}.  Default value = \code{"BFGS"}.
}
  \item{compsqm}{if \code{TRUE}, then standard errors are computed.  Default value = \code{TRUE}.
}
  \item{epsmax}{maximum allowed difference between estimates in subsequent iterations (default = 0.0001).

}
  \item{iterlim}{maximum number of iterations in the maximum likelihood steps (used in \code{nlm} or \code{optim}). Default value = 100.
  }
   \item{ntheta}{number of subdivisions of the round angle, used in the approximation of the integral involved in the likelihood computation of the ETAS model. Default value = 100.
   }



}
\details{
Estimates the components of an ETAS (Epidemic type aftershock sequence)
model for the description of the seismicity of a space-time region. 
Background seismicity is estimated nonparametrically, while triggered seismicity is estimated by MLE.

From version 2.0 of package \code{etasFLP} covariates are allowed to improve the fitting of the triggered part, through the input \code{formula1}, which as a default values of \code{"time ~ magn1.work - 1"}, which corresponds to the previous version of package \code{etasFLP}, that is, magnitude as the only covariate which influence the average number of aftershocks.

The bandwidth of the kernel density estimator is estimated through the 
Forward Likelihood Predictive approach (FLP),
(theoretical reference on Adelfio and Chiodi, 2013) if \code{flp} is set to \code{TRUE}. 
Otherwise the bandwidth is estimated trough the Silverman's rule.
FLP steps for the estimation of nonparametric background component  is alternated with the Maximum Likelihood step for the estimation
of parametric components (only if \code{declustering=TRUE}).
For each event the probability of being a background event or a triggered one is estimated,
according to a declustering procedure in a way similar to the proposal of Zhuang,  Ogata,  and Vere-Jones (2002).


The ETAS model for conditional space time intensity \eqn{\lambda(x,y,t)}{lambda(x,y,t)} is given by:

\deqn{\lambda(x,y,t)=\mu f(x,y)+\kappa_0 \sum_{t_j<t}\frac{ e^{\eta_j}}{(t-t_j +c)^p}
\left\{ \frac{(x-x_j)^2+(y-y_j)^2}{e^{\gamma \ (m_j-m_0)}}+d \right\}^{-q}}{%
lambda(x,y,t)=mu*f(x,y)+
k0\ sum_(t_j<t)(\ e^ eta_j/(t-t_j +c)^p
[ ((x-x_j)^2+(y-y_j)^2)/(e^(gamma \ (m_j-m_0))+d) ]^(-q)}

where \eqn{\eta_j=\sum_{j=1,ncov}betacov_j cov_{ij}}

\eqn{f(x,y)}{f(x,y)} is estimated through a weighted kernel gaussian estimator; if \code{flp} is set to \code{TRUE}
then the bandwidth is estimated through a FLP step.

Weights (computed only if \code{declustering=TRUE}) are given by the estimated probabilities of being  a background event; for the i-th event this is given by
\eqn{\rho_i=\frac{\mu f(x_i,y_i)}{\lambda(x_i,y_i,t_i)}}{rho_i=(mu f(x_i,y_i))/(lambda(x_i,y_i,t_i))}.
The weights \eqn{\rho_i}{rho_i} are updated after a whole iteration. 


 \code{mu} (\eqn{\mu}{mu}) measures the background general intensity (which is assumed temporally homogeneous);

 \code{k0} (\eqn{\kappa_0}{k_0}) is a scale parameter related to the importance of the induced seismicity;

 \code{c} and \code{p} are the characteristic parameters of the seismic
activity of the given region; \code{c} is a shift parameter while \code{p}, which characterizes the pattern of seismicity, is the exponent parameter of the modified Omori law for temporal decay rate of aftershocks; 

\eqn{\eta_j=\sum_{j=1,ncov}betacov_j cov_{ij}} measures the efficiency of an event of a given magnitude in generating aftershock sequences;

\code{d} and \code{q} are two parameters related to the spatial influence of the mainshocks.

Many kinds of ETAS models can be estimated, managing some control input arguments. 
The eight ETAS parameters can be fixed  to some input value, or can be estimated, according to \code{params.ind}:
if \code{params.ind}[i]=FALSE the i-th parameter is kept fixed to its input value, otherwise, if \code{params.ind[i]} \code{ = TRUE},  the i-th parameter is estimated and the input value is used as a starting value. 

By default \code{params.ind=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)}, and so a full 7+\code{ncov} parameters ETAS model will be estimated.

The seven parameters are internally ordered in this way: \code{params} = (\code{mu}, \code{k0}, \code{c}, \code{p},  \code{gamma}, \code{d}, \code{q}); for example
a model with a fixed value \code{p=1} (and  \code{params.ind}[4] = FALSE) can be estimated  and compared with the model where \code{p} is estimated (\code{params.ind}[4]=TRUE);

for example  a 6+\code{ncov} parameters model can be fitted with \code{gamma=0} and \code{params.ind[5]}=\code{FALSE}, 
so that  input must be in this case: \code{params.ind=c(TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE)}; 

if \code{onlytime=TRUE} a time process is fitted to data (with a maximum of 5 parameters), regardless to space location (however the input catalog \code{cat.orig} must contain  three columns named \code{long}, \code{lat, z});

if \code{is.backconstant=TRUE} a process (space-time or time) with a constant background intensity \eqn{\mu}{mu} is fitted;

if \code{mu} is fixed to a very low value a process with very low background intensity is fitted, that is  with only clustered intensity (useful to fit a model to a single  cluster of events).



If \code{flp=TRUE} the bandwidth for the kernel estimation of the background intensity is evaluated maximizing the  
Forward Likelihood Predictive (FLP) quantity, given by
(Chiodi, Adelfio, 2011; Adelfio, Chiodi, 2013):

\deqn{FLP_{k_1,k_2}(\hat{\boldsymbol{\psi}})\equiv\sum_{k=k_1}^{n-1}\delta_{k,{k+1}}(\hat{\boldsymbol{\psi}}(H_{t_k});
H_{t_{k+1}})}

with \eqn{k_1=\frac{n}{2},k_2=n-1} and where \eqn{\delta_{k,k+1}(\hat{\boldsymbol{\psi}}(H_{t_k});
H_{t_{k+1}})}{(d_(k,k+1))(psi(Ht);H(t+1))} is 
the \emph{predictive information}
of the first \eqn{k} observations on the \eqn{k+1}-th observation, and is so defined:

\deqn{\delta_{k,k+1}(\hat{\boldsymbol{\psi}}(H_{t_k});
H_{t_{k+1}})\equiv \log L(\hat{\boldsymbol{\psi}}(H_{t_k}); H_{t_{k+1}} )-\log L(\hat{\boldsymbol{\psi}}(H_{t_k});H_{t_k})}{d_(k,k+1)=log L(psi^(H_k);H_(k+1)-log L(psi^(H_k);H_k)}

where \eqn{H_k} is the  history of the process until time \eqn{t_k} and 
\eqn{\hat{\boldsymbol{\psi}}(H_{t_k})}{psi^(H_{t_k})} is an estimate based only on history until the \eqn{k-th} observation.


In the ML step, the vector of parameter \eqn{\theta=(\mu, \kappa_0, c , p, \alpha, \gamma, d, q)} is estimated maximizing the sample log-likelihood given by:


\deqn{\log L(\boldsymbol{\theta}; H_{t_n}) = \sum_{i=1}^{n}
\log \lambda(x_i,y_i,t_i; \boldsymbol{\theta})-
\int_{T_0}^{T_{max}} \int \int_{\Omega_{(x,y)}}\,
\lambda(x,y,t;\boldsymbol{\theta})\,d x \, d y \,d t
}






}

\value{
returns an object of class \code{etasclass}.

The main items of the output are:
  \item{this.call}{reports the exact call of the function}
  \item{params.ind}{indicates which parameters have been estimated (see details)}
  \item{params}{ML estimates of the ETAS parameters.}
  \item{sqm}{Estimates of standard errors of the ML estimates of the ETAS parameters (\code{sqm[i]}=0 if \code{params.ind[i]}=\code{FALSE} or where the hessian is not computed or near to singularity). }
  \item{AIC.iter}{AIC values at each iteration.}
  \item{hdef}{final bandwidth used for the kernel estimation of background spatial intensity (however estimated, with \code{flp=TRUE} or \code{flp=FALSE}).}
  \item{rho.weights}{Estimated probability for each event to be a background event (\eqn{\rho}).}
  \item{time.res}{rescaled time residuals (for time processes only).}
  \item{params.iter}{A matrix with estimates values at each iteration.}
  \item{sqm.iter}{A matrix with the estimates of the standard errors at each iteration.}
  \item{rho.weights.iter}{A matrix with the values of \code{rho.weights} at each iteration.}
  \item{l}{A vector with estimated intensities, corresponding to observed points}

\code{summary},
\code{print} and
\code{plot} methods are defined for an object of class \code{etasclass} to obtain main output. 

A \code{profile} method (\code{\link{profile.etasclass}}
) is also defined  to make approximate inference on a single parameter   
}
\note{In this version the x-y space region, where the point process is defined, is  a rectangle embedding the catalog values.

The optimization algorithm  depends on the choice of initial values. Some default guess choice is performed 
inside the function for parameters without input starting values; the function \code{etas.starting} gives rough first guess for initial values. If convergence problem are experienced, a useful strategy can be starting  with an higher magnitude threshold value \eqn{m_0}{m0} (that is, with a smaller catalog with bigger earthquakes), and then using this first output as starting  guess for a running with a lower magnitude threshold value \eqn{m_0}{m0}. 
In this trial executions avoid declustering (\code{declustering=FALSE}) or at least use a small value of \code{ndeclust}; small values of \code{iterlim} and \code{ntheta} can speed first executions.

Quicker executions are obtained using smaller values of \code{iterlim} and \code{ntheta} in the input.

Also a first execution with \code{is.backconstant = TRUE}, to fit a first approximation model with constant background, can be useful.

Some other useful information can be obtained estimating a pure time process, that can give a good guess at least for some  parameters, like  \eqn{\mu, \kappa_0, \alpha,c,p}{mu, k0, a, c, p}.


Input times are expected in days, and so final intensities are expected number of events per day. If input values are in seconds, then set \code{sectoday=TRUE}
}

\seealso{ \code{\link{eqcat}}, \code{\link{plot.etasclass}},\code{\link{print.etasclass}}, \code{\link{summary.etasclass}}, \code{\link{profile.etasclass}}, \code{\link{etas.starting}}}

\references{

Adelfio G, Chiodi M (2015). Alternated Estimation in Semi-Parametric Space-Time
Branching-Type Point Processes with Application to Seismic Catalogs. \emph{Stochastic Environmental Research and Risk Assessment}, \bold{29}(2), 443-450. doi:10.1007/s00477-014-0873-8.

Adelfio G, Chiodi M (2015). FLP Estimation of Semi-Parametric Models for Space-Time
Point Processes and Diagnostic Tools. \emph{Spatial Statistics}, \bold{14}(B), 119-132. doi:10.1016/j.spasta.2015.06.004.

Adelfio, G. and Chiodi, M. (2013) Mixed estimation technique in semi-parametric space-time point processes for earthquake description. 
\emph{Proceedings of the 28th International Workshop on Statistical Modelling 8-13 July, 2013, Palermo} (Muggeo V.M.R., Capursi V., Boscaino G., Lovison G., editors). Vol. \bold{1} pp.65-70. 

Chiodi, M. and Adelfio, G., (2011) Forward Likelihood-based predictive approach for space-time processes. \emph{Environmetrics}, vol. \bold{22} (6), pp. 749-757. DOI:10.1002/env.1121.

Chiodi, M. and Adelfio, G., (2017) Mixed Non-Parametric and Parametric Estimation Techniques in {R} Package {etasFLP} for Earthquakes' Description. \emph{Journal of Statistical Software}, vol. \bold{76} (3), pp. 1-29.
 DOI: 10.18637/jss.v076.i03.


Zhuang, J., Ogata, Y.  and Vere-Jones, D.
Stochastic declustering of space-time earthquake occurrences.
\emph{Journal of the American Statistical Association},
\bold{97},  369--379 (2002). DOI:10.1198/016214502760046925.
}

\author{
Marcello Chiodi, Giada Adelfio}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\dontrun{
data("italycatalog")
# load a sample catalog of the italian seismicity
esecov1<-etasclass(cat.orig = catalog.withcov, magn.threshold = 2.5, magn.threshold.back = 3.9, 
mu = 0.3, k0 = 0.02, c = 0.015, p = 0.99, gamma = 0, d = 1, 
q = 1.5, params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, 
TRUE), formula1 = "time ~  magn1.work- 1", declustering = TRUE, 
thinning = FALSE, flp = TRUE, ndeclust = 15, onlytime = FALSE, 
is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE, 
compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36)

# execution of etasclass for events with minimum magnitude of 3.0. 
# The events with magnitude at least 3.5 are used to build a first approximation
# for the background intensity function
# (magn.threshold.back=3.5)
# The magnitude effect is given by the covariate magn1.work 
# in the formula "time ~  magn1.work- 1"
# magn1.work is the internal name for magn1-magn.threshold
# print method for the etasclass object

print(esecov1)

> print(esecov1)
Call: 
 
etasclass(cat.orig = catalog.withcov, magn.threshold = 2.5, magn.threshold.back = 3.9, 
    mu = 0.3, k0 = 0.02, c = 0.015, p = 0.99, gamma = 0, d = 1, 
    q = 1.5, params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, 
        TRUE), formula1 = "time ~  magn1.work- 1", declustering = TRUE, 
    thinning = FALSE, flp = TRUE, ndeclust = 15, onlytime = FALSE, 
    is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE, 
    compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36)

 
 
 
Number of observations             2226 
ETAS Parameters: 
        mu         k0          c          p      gamma          d          q 
  0.667509   0.022393   0.014769   1.110059   0.000000   1.905461   1.947223 
magn1.work 
  0.740109 
# summary method for more informative output etasclass object

summary(esecov1)
# plot results with maps of intensities and diagnostic tools 

plot(esecov1)

## an application with 5 covariates
esecov5<-etasclass(cat.orig = catalog.withcov, magn.threshold = 2.5, magn.threshold.back = 3.9, 
mu = 0.3, k0 = 0.02, c = 0.015, p = 0.99, gamma = 0, d = 1, 
q = 1.5, params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), 
formula1 = "time ~ z + magn1.work +nstaloc_rev +min_distance_rev+distmin- 1", 
declustering = TRUE, thinning = FALSE, flp = TRUE, ndeclust = 15, 
onlytime = FALSE, is.backconstant = FALSE, sectoday = FALSE, 
usenlm = TRUE, compsqm = TRUE, epsmax = 1e-04, iterlim = 100, 
ntheta = 36)

## print results, more out put with summary
 print(esecov5)
Call: 
 
etasclass(cat.orig = catalog.withcov, magn.threshold = 2.5, magn.threshold.back = 3.9, 
    mu = 0.3, k0 = 0.02, c = 0.015, p = 0.99, gamma = 0, d = 1, 
    q = 1.5, params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, 
        TRUE), formula1 = "time ~ z + magn1.work +nstaloc_rev +min_distance_rev+distmin- 1", 
    declustering = TRUE, thinning = FALSE, flp = TRUE, ndeclust = 15, 
    onlytime = FALSE, is.backconstant = FALSE, sectoday = FALSE, 
    usenlm = TRUE, compsqm = TRUE, epsmax = 1e-04, iterlim = 100, 
    ntheta = 36)

 
 
 
Number of observations             2226 
ETAS Parameters: 
              mu               k0                c                p 
        0.705351         0.073070         0.019396         1.154186 
           gamma                d                q                z 
        0.000000         1.942929         2.004915        -0.041256 
      magn1.work      nstaloc_rev min_distance_rev          distmin 
        1.157698        -0.009010        -0.011020        -1.826717 

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ETAS}
\keyword{earthquake}
\keyword{kernel}
\keyword{covariates}
\keyword{flp}

% __ONLY ONE__ keyword per line
