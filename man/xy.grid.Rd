\name{xy.grid}
\alias{xy.grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a 2-d grid}
\description{
Creates a 2-d grid.}
\usage{
xy.grid(rangex, rangey, nx, ny = nx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rangex}{A length 2 numeric vector: the range of the x-variable.}
  \item{rangey}{A length 2 numeric vector: the range of the y-variable.}
  \item{nx}{The number of points of the grid in the x-direction.}
  \item{ny}{The number of points of the grid in the y-direction.}
}
\value{A grid of the coordinates of  \code{nx*ny} points on the x-y plane, expanded in a matrix of \code{nx*ny} rows and 2 columns: 
a row gives the (x,y) coordinates of a point.}
\examples{
xy.grid(c(3,7),c(11,17),nx=5,ny=4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
