summary.eqcat <-
function(object,extended=TRUE,...){
x=object
# summary for earthquake.catalogs
cat("Earthquake Catalog of ",length(x$long),"  events. ", "\n")
cat("------------------------------------------------------------------", "\n")
cat("Catalog extension: ", "\n")
cat("Longitude range  ",round(range(x$long),2),"\n")
cat("Latitude  range  ",round(range(x$lat),2),"\n")
cat("Depth     range  ",round(range(x$z),2),"\n")
cat("Time      range  ",round(range(x$time),2),"\n")
cat("------------------------------------------------------------------", "\n")
cat("Magnitude range  ",round(range(x$magn1),2),"\n")
if (extended){
cat("------------------------------------------------------------------", "\n")
print("Summary of time differences")
print(summary(diff(x$t)))
cat("------------------------------------------------------------------", "\n")
print("Magnitude Distribution")
print(summary(diff(x$magn1)))
MLA.freq(x$magn1)
}
}
