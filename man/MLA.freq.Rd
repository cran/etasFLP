\name{MLA.freq}
\alias{MLA.freq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Display a pretty frequency table
}
\description{Display a pretty frequency table. It is only a wrapper to the function \code{table}
but with a richer  output, at least for numerical variables.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MLA.freq(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector.
%%     ~~Describe \code{x} here~~
}
}
\details{The output gives the different kinds of frequencies and cumulated frequencies:
single frequencies, cumulated and back cumulated (absolute and relatives).
%%  ~~ If necessary, more details than the description above ~~
}
\value{return a matrix with 7 columns: the modaldistinct values of \code{x}, frequencies, relative frequencies, cumulated frequencies,
cumulated relative frequencies,  back cumulated frequencies
and back cumulated relative frequencies.
}
\examples{
%\dontrun{

x=trunc(runif(1000)*10)
MLA.freq(x)

data(italycatalog)
MLA.freq(italycatalog$magn1)

%}
}
\author{
Marcello Chiodi}



