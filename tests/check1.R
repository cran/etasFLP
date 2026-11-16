library(etasFLP)
data(italycatalog)

fit <- etasclass(italycatalog[1:150, ], magn.threshold = 3.0,
                 magn.threshold.back = 3.5, k0 = 0.005, c = 0.005, p = 1.01,
                 gamma = 0.6, q = 1.52, d = 1.1, ndeclust = 1,
                 sectoday = TRUE, epsmax = 0.001, trace = FALSE)

stopifnot(inherits(fit, "etasclass"))
ref.params <- c(mu = 0.511418, k0 = 0.005861596, c = 0.04996627,
                p = 1.348313, gamma = 0.882482, d = 1.124963, q = 1.768607)
stopifnot(isTRUE(all.equal(fit$params, ref.params, tolerance = 1e-3)))
stopifnot(isTRUE(all.equal(fit$hdef, c(91.97022, 76.10657), tolerance = 1e-3)))
stopifnot(is.null(fit$parallel), is.null(fit$fastML.))

out <- capture.output(print(fit))
stopifnot(length(out) > 0)
out <- capture.output(summary(fit))
stopifnot(length(out) > 0)

fit2 <- update(fit, iterlim = 3)
stopifnot(inherits(fit2, "etasclass"), length(fit2$params) == 7)

cmp <- compare.etasclass(fit, fit2)
stopifnot(is.finite(cmp$AIC))

n0 <- nrow(fit$cat)
fit$cat.orig <- italycatalog[1:160, ]
fit3 <- suppressWarnings(timeupdate.etasclass(fit))
stopifnot(inherits(fit3, "etasclass"), nrow(fit3$cat) > n0)

d1 <- daily.etasclass(fit2, ngrid = 41, nclass = 10, tfixed = 1)
stopifnot(!is.null(d1))

mla <- MLA.freq(italycatalog$magn1)
stopifnot(is.data.frame(mla), sum(mla$freq) == nrow(italycatalog))
