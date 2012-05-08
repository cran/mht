## .onLoad <- function(lib, pkg)
##   {
##     require(methods)
##     require(emulator)
##   }

##.onAttach <- function(libname, pkgname){
##  cat("-----------------",
##      "Loaded mht",as.character(packageDescription("mht")[["Version"]]),
##      "-----------------",
##      sep = "\n")
##}

.onAttach <- function(libname, pkgname){ packageStartupMessage("Loaded mht ",as.character(packageDescription("mht")[["Version"]]))}