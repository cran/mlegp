.First.lib <- function(lib, pkg) {
  library.dynam("mlegp", pkg, lib)
  
  #  load package adapt, or display warning if it does not exist
#if (FALSE) {
#  success = suppressWarnings(require(adapt, quietly=TRUE, warn.conflicts=FALSE))
#        if (success) {
#                suppressWarnings(detach(package:adapt))
#                library(adapt)   ## will show warnings, if there are any
#                cat("loading recommended package: adapt\n")
#        }
#        else {
#                cat("warning: package adapt could not be loaded\n")
#                cat("the following functions will not be available: FANOVADecomposition, plotInteractionEffects\n")
#       }
#  }
}

