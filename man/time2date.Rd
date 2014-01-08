\name{time2date}
\alias{time2date}
\alias{timecharunique2seq}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Date time conversion tools}
\description{
Date time conversion tools, useful in connection with package etasFLP for earthquake description. Base data is Jan. 1st 1900}
\usage{
time2date(t)

timecharunique2seq(timestring)
}
%- maybe also 'usage' for other objects documented here.
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{A time:  seconds elapsed from 1900-1-1.}
%%     ~~Describe \code{t} here~~

  \item{timestring}{A time string.
%%     ~~Describe \code{t} here~~
}
}


\details{\code{time2date} converts sequential time into character string;
\code{timecharunique2seq}  converts character times of catalogs into sequential time (secs elapsed from the base date): the input is a single string.}

\value{\code{time2date} returns a character string; \code{timecharunique2seq} returns a list:
\item{char}{the input string.}
\item{sec}{secs elapsed from the base date.}
\item{day}{days elapsed from the base date.}
}
%\references{}
\author{Marcello Chiodi}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{}
%\examples{
  

\examples{
\dontrun{
tchar="1960-11-06 11:09:35.000"
tsec =timecharunique2seq(tchar)[["sec"]]
time2date(tsec)
}}
\keyword{date}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
