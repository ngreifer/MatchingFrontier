#Needed to pass CRAN checks for NSE calls with ggplot2
utils::globalVariables(c("Val", "Covariate", "Est", "Var"))

#Used to load backports functions. No need to touch, but must always be included somewhere.
.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
}

.onAttach <- function(libname, pkgname) {
  pkgLib <- dirname(system.file(package = pkgname))
  version <- packageDescription(pkgname, lib.loc = pkgLib)$Version
  BuildDate <- packageDescription(pkgname, lib.loc = pkgLib)$Date

  foo <- paste0(" ", pkgname, " (Version ", version, ", Build Date: ", format(BuildDate, "%F"), ")")
  packageStartupMessage(foo)
}
