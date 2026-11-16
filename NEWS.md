# etasFLP 2.4.0

* Faster computations: the Fortran routines for the likelihood and for the
  intensities used in the plots have been rewritten, moving out of the loops
  the terms that depend on a single event and using `exp(-p*log(.))` in place
  of powers with real exponents. In our benchmarks the computing time is about
  45-50% shorter. Estimates can differ from those of version 2.3.0 in the
  seventh digit or so, because of the different rounding.
* Removed the unused, never activated options `parallel` and `fastML` (with
  the corresponding R code and Fortran routines, including the OpenMP
  directives). Objects of class `etasclass` no longer contain the components
  `parallel` and `fastML.`.
* Removed the unused internal functions `etas.mod2NEW_withoutfastML` and
  `fastML.init`.
* Fixed `timeupdate.etasclass`, which stopped with an error.
* Fixed `compare.etasclass`: the difference of AIC values was not computed
  with the objects of the current version.
* `plot.etasclass` and `plot.eqcat`: new windows are opened only in interactive
  sessions (otherwise all graphs go to the current device); graphical
  parameters changed by `plot.etasclass` are restored; the pdf device opened
  with `pdf = TRUE` is closed also in case of error.
  The list returned by `plot.etasclass` is now returned invisibly.
* Removed a warning in `etasclass` due to the recycling of a vector of
  weights.
* Documentation revised: wrong default values (`iterlim`, `ntheta`), the
  argument `cat.back`, the description of the returned values and several
  typos.
* Examples are now run (the longest ones in `\donttest`); new examples have
  been added.
* `NAMESPACE`: imports restricted to the functions actually used;
  `DESCRIPTION` updated.
* Added tests.

# etasFLP 2.3.0

* See the file `ChangeLog` for the changes of this and earlier versions.
