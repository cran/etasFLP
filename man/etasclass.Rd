\name{etasclass}
\alias{etasclass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mixed estimation of an ETAS model}
\description{\code{etasclass} is the main function of the package \code{etasFLP}.

Performs the estimation of the components of the ETAS (Epidemic type aftershock sequences) model for the description of the seismicity in a space-time region. Background seismicity is estimated non-parametrically, while triggered seismicity is estimated by MLE. In particular also the bandwidth for a kernel smoothing can be estimated through the Forward Likelihood Predictive approach (FLP). For each event the probability of being a background event or a triggered one is estimated. 

An ETAS with up to 8 parameters can be estimated, with several options and different methods.

Returns an \code{etasclass} object, for which \code{plot}, \code{summary}, \code{print} and \code{profile} methods are defined. 
}
\usage{etasclass(cat.orig, 
    magn.threshold=2.5, magn.threshold.back=magn.threshold+2,
    mu=1,k0=1,c=0.5,p=1.01,	a=1.2,gamma=.5,d=1.,q=1.5, params.ind=replicate(8,TRUE),
%    kern.var=FALSE, alpha=0.5,
    hdef=c(1,1),
declustering=TRUE,thinning=FALSE,
    flp=TRUE, m1=NULL, ndeclust=5, onlytime=FALSE,is.backconstant=FALSE,
    w=replicate(nrow(cat.orig),1),
##### end of  main input arguments. 
##### Control and secondary arguments:
    description="", cat.back=NULL, back.smooth=1.0,
    sectoday=TRUE,longlat.to.km=TRUE,
    usenlm=TRUE, method	="BFGS", compsqm=TRUE,
    epsmax=0.0001, iterlim=100, ntheta=100)
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

\emph{Values for the 8 parameters of the ETAS   model} (starting values or fixed values according to  \code{params.ind}):
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
  \item{a}{Parameter 5 (\eqn{\alpha}) of the ETAS model;
efficiency of an event of given magnitude in generating aftershocks;  see details. Default value = 1.2.
%%     ~~Describe \code{a} here~~
}
  \item{gamma}{Parameter 6 (\eqn{\gamma}) of the ETAS model; together with  \code{a} is related to the efficiency of an event of given magnitude 
  in generating aftershocks; see details. Default value = 0.5.
%%     ~~Describe \code{gamma} here~~
}
  \item{d}{Parameter 7 of the ETAS model; parameter related to the spatial influence of the
mainshock; see details. Default value = 1.
%%     ~~Describe \code{d} here~~
}
  \item{q}{Parameter 8 of the ETAS model; parameter related to the spatial influence of the
mainshock; see details. Default value = 1.5.
%%     ~~Describe \code{q} here~~
}
\emph{End of model pararameter input}
  \item{params.ind}{vector of 8 logical values: \code{params.ind[i] = TRUE} means that the i-th parameter must be estimated.  \code{params.ind[i] = FALSE}  means that the i-th parameter is fixed to its input value  (the order of parametrs is: \code{mu}, \code{k0}, \code{c}, \code{p}, \code{a}, \code{gamma}, \code{d}, \code{q}).
  Default value = \code{replicate(8,TRUE)}, that is, \code{etasclass} estimates all parameters.
}

\emph{Flags for the kind of declustering and smoothing}:
%  \item{kern.var}{if \code{TRUE} the background seismicity is estimated through a variable metric anisotropic kernel. yet not used in version 1.1.0. Default value = \code{FALSE}.}
  \item{hdef}{Starting values for the \code{x,y} bandwidths used in the kernel estimator of background seismicity. Default value = \code{1,1}.
}
%  \item{alpha}{if \code{kern.var=TRUE} indicates the mixing proportion in the background seismicity between  variable metric and fixed metric kernel.yet not used in version 1.1.0. Default value = \code{0.5}.}

  \item{declustering}{if \code{TRUE} the catalog is iteratively declustered to optimally estimate the background intensity (through thinning, if \code{thinning=TRUE}, or through weighting if \code{thinning=FALSE}). Default value = \code{TRUE}.
}
  \item{thinning}{if \code{thinning=TRUE} a background catalog is obtained sampling from the original catalog with probabilities estimated during the iterations. Default value =\code{FALSE}.
}
  \item{flp}{if \code{flp=TRUE} then background seismicity is estimated through Forward Likelihood Predictive (see details). Otherwise the Silverman rule is used. Default value =\code{TRUE}.
%%     ~~Describe \code{flp} here~~
} 
  \item{m1}{Used only if \code{flp=TRUE}. Indicates the range of points used for the FLP steps. See details. If missing it is set  to \code{nrow(cat)/2}.
%%     ~~Describe \code{flp} here~~
} 
\item{ndeclust}{maximum number of iterations for the general declustering procedure. Default=5.
%%     ~~Describe \code{ndeclust} here~~
}

\item{onlytime}{if \code{TRUE}  then a time process is fitted to data , regardless to space location (in this case \code{is.backconstant} is set to \code{TRUE} and \code{declustering}, \code{flp} are set to \code{FALSE}). Default value = \code{FALSE}.
%%     ~~Describe \code{onlytime} here~~
}
  \item{is.backconstant}{if \code{TRUE}  then background seismicity is assumed to be homogeneous in space  (and \code{declustering, flp} are set to \code{FALSE}).  Default value = \code{FALSE}.
%%     ~~Describe \code{is.backconstant} here~~
}
  \item{w}{initial weights}
\emph{Other control parameters}:
  \item{description}{a description string used for the output. Default value = "".
%%     ~~Describe \code{description} here~~
}
  \item{cat.back}{ external catalog used for the estimation of the background seismicity. 
 Default value = \code{NULL}.
}
%If \code{declustering=TRUE} it is used only as a first approximation.
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
Estimates the components of an ETAS (Epidemic type aftershock sequences)
model for the description of the seismicity of a space-time region. 
Background seismicity is estimated nonparametrically, while triggered seismicity is estimated by MLE.

The bandwidth of the kernel density estimator is estimated through the 
Forward Likelihood Predictive approach (FLP),
(theoretical reference on Adelfio and Chiodi, 2013) if \code{flp} is set to \code{TRUE}. 
Otherwise the bandwidth is estimated trough  Silverman's rule.
FLP steps for the estimation of nonparametric background component  is alternated with the Maximum Likelihood step for the estimation
of parametric components (only if \code{declustering=TRUE}).
For each event the probability of being a background event or a triggered one is estimated,
according to a declustering procedure in a way similar to the proposal of Zhuang,  Ogata,  and Vere-Jones (2002).

The ETAS model for conditional space time intensity \eqn{\lambda(x,y,t)}{lambda(x,y,t)} is given by:

\deqn{\lambda(x,y,t)=\mu f(x,y)+\sum_{t_j<t}\frac{\kappa_0\ e^{(\alpha-\gamma) \ (m_j-m_0)}}{(t-t_j +c)^p}
\left\{ \frac{(x-x_j)^2+(y-y_j)^2}{e^{\gamma \ (m_j-m_0)}}+d \right\}^{-q}}{%
lambda(x,y,t)=mu*f(x,y)+
sum_(t_j<t)(k0\ e^(a-gamma) \ (m_j-m_0))/(t-t_j +c)^p
[ ((x-x_j)^2+(y-y_j)^2)/(e^(gamma \ (m_j-m_0))+d) ]^(-q)}


\eqn{f(x,y)}{f(x,y)} is estimated through a weighted kernel gaussian estimator; if \code{flp} is set to \code{TRUE}
then the bandwidth is estimated through a FLP step.

Weights (computed only if \code{declustering=TRUE}) are given by the estimated probabilities of being  a background event; for the i-th event this is given by
\eqn{\rho_i=\frac{\mu f(x_i,y_i)}{\lambda(x_i,y_i,t_i)}}{rho_i=(mu f(x_i,y_i))/(lambda(x_i,y_i,t_i))}.
The weights \eqn{\rho_i}{rho_i} are updated after a whole iteration. 


 \code{mu} (\eqn{\mu}{mu}) measures the background general intensity (which is assumed temporally homogeneous);

 \code{k0} (\eqn{\kappa_0}{k_0}) is a scale parameter related to the importance of the induced seismicity;

 \code{c} and \code{p} are the characteristic parameters of the seismic
activity of the given region; \code{c} is a shift parameter while \code{p}, which characterizes the pattern of seismicity, is the exponent parameter of the modified Omori law for temporal decay rate of aftershocks; 

\code{a} (\eqn{\alpha}{alpha}) and \code{gamma} (\eqn{\gamma}{gamma}) measure the efficiency of an event of given magnitude in generating aftershock sequences;

\code{d} and \code{q} are two parameters related to the spatial influence of the mainshocks.

Many kinds of ETAS models can be estimated, managing some control input arguments. 
The eight ETAS parameters can be fixed  to some input value, or can be estimated, according to \code{params.ind}:
if \code{params.ind}[i]=FALSE the i-th parameter is kept fixed to its input value, otherwise, if \code{params.ind[i]} \code{ = TRUE},  the i-th parameter is estimated and the input value is used as a starting value. 

By default \code{params.ind=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)}, and so a full 8 parameters ETAS model will be estimated.

The eight parameters are internally ordered in this way: \code{params} = (\code{mu}, \code{k0}, \code{c}, \code{p}, \code{a}, \code{gamma}, \code{d}, \code{q}); for example
a model with a fixed value \code{p=1} (and  \code{params.ind}[4] = FALSE) can be estimated  and compared with the model where \code{p} is estimated (\code{params.ind}[4]=TRUE);

for example  a 7 parameters model can be fitted with \code{gamma=0} and \code{params.ind[6]}=\code{FALSE}, 
so that  input must be in this case: \code{params.ind=c(TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE)}; 

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
  \item{sqm}{Estimates of standard errors of the ML estimates of the ETAS parameters (\code{sqm[i]}=0 if \code{params.ind[i]}=\code{FALSE} and for the situation where hessian is not computed or near to singularity). }
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
\note{In this first version the x-y space region, where the point process is defined, is  a rectangle embedding the catalog values.

The optimization algorithm  depends on the choice of initial values. Some default guess choice is performed 
inside the function for parameters without input starting values. If convergence problem are experienced, a useful strategy can be to start  with an high magnitude threshold value \eqn{m_0}{m0} (that is, with a smaller catalog with bigger earthquakes), and then using this first output as starting  guess for a running with a lower magnitude threshold value \eqn{m_0}{m0}. 
In this trial executions avoid declustering (\code{declustering=FALSE}) or at least use a small value of \code{ndeclust}; small values of \code{iterlim} and \code{ntheta} can speed first executions.

Quicker executions are obtained using smaller values of \code{iterlim} and \code{ntheta} in the input.

Also a first execution with \code{is.backconstant = TRUE}, to fit a first approximation model with constant background, can be useful.

Some other useful information can be obtained estimating a pure time process, that can give a good guess at least for some  parameters, like  \eqn{\mu, \kappa_0, \alpha,c,p}{mu, k0, a, c, p}.


Input times are expected in days, and so final intensities are expected number of events per day. If input values are in seconds, then set \code{sectoday=TRUE}
}

\seealso{ \code{\link{eqcat}}, \code{\link{plot.etasclass}}, \code{\link{summary.etasclass}}, \code{\link{profile.etasclass}}}

\references{


Adelfio, G. and Chiodi, M. (2014) Alternated estimation in semi-parametric space-time branching-type point processes with application to seismic catalogs. \emph{Stochastic Environmental Research and Risk Assessment},
doi={10.1007/s00477-014-0873-8}

Adelfio, G. and Chiodi, M. (2013) Mixed estimation technique in semi-parametric space-time point processes for earthquake description. 
\emph{Proceedings of the 28th International Workshop on Statistical Modelling 8-13 July, 2013, Palermo} (Muggeo V.M.R., Capursi V., Boscaino G., Lovison G., editors). Vol. 1. pp.65-70. 

Chiodi, M. and Adelfio, G., (2011) Forward Likelihood-based predictive approach for space-time processes. \emph{Environmetrics}, vol. 22 (6), pp. 749-757.

Zhuang, J., Ogata, Y.  and Vere-Jones, D.
Stochastic declustering of space-time earthquake occurrences.
\emph{Journal of the American Statistical Association},
\bold{97},  369--379 (2002).
}

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
description = "etas flp",sectoday = TRUE, usenlm = TRUE,epsmax = 0.001)


# execution of etasclass for events with minimum magnitude of 3.1. 
# The events with magnitude at least 3.5 are used to build a first approximation
# for the background intensity function
# (magn.threshold.back=3.5)

# summary method for the etasclass object

summary(etas.flp)
Call: 
 
etasclass(cat.orig = italycatalog, magn.threshold = 3.1, magn.threshold.back = 3.5, 
    k0 = 0.005, c = 0.005, p = 1.01, a = 1.05, gamma = 0.6, d = 1.1, 
    q = 1.52, params.ind = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
        TRUE, TRUE), declustering = TRUE, thinning = FALSE, flp = TRUE, 
    ndeclust = 15, onlytime = FALSE, is.backconstant = FALSE, 
    description = "etas flp", sectoday = TRUE, usenlm = TRUE, 
    epsmax = 0.001)

 
 
etas flp 
Execution started:                  2015-04-14 17:08:50 
Elapsed time of execution (hours)   0.1735942 
Number of observations             1700 
Magnitude threshold                3.1 
declustering                        TRUE 
Number of declustering iterations   4 
Kind of declustering                weighting 
flp                                 TRUE 
sequence of AIC values for each iteration 
40444.81 39042.21 39025.69 39025.1 
 
------------------------------------------------------- 
 
ETAS Parameters: 
            Estimates       std.err.
mu           0.299153       0.010178
k0           0.008847       0.002832
c            0.012754       0.002753
p            1.149591       0.020213
a            1.640813       0.070028
gamma        0.932248       0.094336
d            2.010304       0.384207
q            1.926604       0.089574
------------------------------------------------------- 

# plot results with maps of intensities and diagnostic tools 

plot(etas.flp)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ETAS}
\keyword{earthquake}
\keyword{kernel}
\keyword{flp}

% __ONLY ONE__ keyword per line
