% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update.etasclass.R
\name{update.etasclass}
\alias{update.etasclass}
\title{update.etasclass}
\usage{
\method{update}{etasclass}(object, ...)
}
\arguments{
\item{object}{an etasclass object  obtained from etasflp versions 2.2 or newer that will be updated.}

\item{...}{optional arguments that will override the corresponding arguments in \code{object}}
}
\value{
un updated \code{etasclass} object
}
\description{
New in version 2.2. A method update for etasclass objects: a very experimental version that can be used only on etasclass objects obtained from etasflp versions 2.2 or newer.
}
\details{
It is a beta version. The catalog must be the same, and options in "..." must leave unchanged the number of observations used for estimation. Arguments given in "..." will override arguments already present in \code{object}. Not all arguments are suitable for updating: among them \code{formula} and \code{params.ind} should not be included in "..." list (to update such parameters it is better to assign them to a variable and then pass the variable name) . A new etasclass execution will start, using as arguments  values of input \code{object}, eventually integrated with the list in "...". Tipically a first execution can be given with low values of iterlim, ndeclust, ntheta and high values of epsmax (e.g. iterlim=5, ndeclust=1, ntheta=24, epsmax=0.01), to obtain good starting values for parameters and for weights. Then an update can be run with better values such as iterlim=50, ndeclust=10, ntheta=60, epsmax=0.0001.
}
\seealso{
\code{\link{timeupdate.etasclass}}
}
