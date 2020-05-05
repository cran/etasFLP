\name{catalog.withcov}
\docType{data}
\alias{catalog.withcov}
\title{Small sample catalog of italian earthquakes with covariates}
\description{
A small sample catalog of  italian earthquakes of magnitude at least 2.5 from May 2012 to May 2016 with extra information.}
\usage{catalog.withcov}
\format{a data matrix with 2,226 observations and 11 variables:
\code{Date}, \code{id_ev},  \code{time}, \code{lat}, \code{long}, \code{z}, \code{magn1}, \code{err_h_rev}, \code{min_distance_rev}, \code{rms_rev}, \code{nstaloc_rev}           .
}
\details{
A small sample catalog of  italian earthquakes of magnitude at least 2.5 from May 2012 to May 2016 with extra information.
\code{Date} is the date, \code{id_ev} is an identifier of the events,  \code{time}, \code{lat}, \code{long}, \code{z}, \code{magn1} are time (days from 1900-01-01), latitude, longitude, depth, magnitude, \code{err_h_rev} is an estimate of hypocentral uncertainty,  \code{min_distance_rev} is the distance from the nearest station, \code{rms_rev} is a measure of the quality of the location, \code{nstaloc_rev} is the number of stations that registered the event, \code{distmin} is the
distance from the nearest fault}    


\examples{
data(catalog.withcov)
str(catalog.withcov)
}
\source{INGV (Istituto Nazionale di Geofisica e Vulcanologia) ISIDE Data Base.}
\references{INGV  home page: \url{http://www.ingv.it}.}
\keyword{datasets}
\keyword{earthquake}
