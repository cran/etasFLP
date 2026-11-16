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
  \item{object}{an \code{etasclass} object to pass to \code{summary}.}
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
data("italycatalog")
## a small sample of the catalog keeps the execution time short:
## please, use larger values of n to obtain reasonable datasets and estimates
sub <- italycatalog[1:150, ]
fit <- etasclass(sub, magn.threshold = 3.0, magn.threshold.back = 3.5,
    k0 = 0.005, c = 0.005, p = 1.01, gamma = 0.6, q = 1.52, d = 1.1,
    ndeclust = 1, sectoday = TRUE, epsmax = 0.001, trace = FALSE)

summary(fit)
## more output, with the summary of the weights of background seismicity
summary(fit, full = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ETAS}
\keyword{earthquake}
\keyword{summary}

% __ONLY ONE__ keyword per line
