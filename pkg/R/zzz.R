# (c) William R. Engels, 2014

#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname){
	version <- packageDescription("HWxtest", fields = "Version")
	packageStartupMessage("------------------------------------------------")
	packageStartupMessage("                  HWxtest")
	packageStartupMessage(paste("               version", version))
	packageStartupMessage("Please see the package vignette for instructions")
	packageStartupMessage("------------------------------------------------")
}