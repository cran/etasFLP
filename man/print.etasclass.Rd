\name{print.etasclass}
\alias{print.etasclass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print method for etasclass objects}
\description{Print method for an object of class \code{etasclass}.

Gives some information on the execution and gives estimates of the ETAS parameters.
}
\usage{
\method{print}{etasclass}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \code{etaclass} object.}
  \item{...}{other arguments.}
}
\details{
Print brief information about an object of class \code{etasclass}. More output is obtained with \code{summary}.
}

\value{
Displays parameters estimates and information on the execution of the \code{etasclass}
estimation process. Displays also  the exact call of the function that generated \code{etasclass}
}


\seealso{ \code{\link{etasclass}},\code{\link{eqcat}}, \code{\link{profile.etasclass}}}


\author{
Marcello Chiodi, Giada Adelfio}

%% ~Make other sections like Warning with \section{Warning }{....} ~


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ETAS}
\keyword{earthquake}
\keyword{print}

% __ONLY ONE__ keyword per line
