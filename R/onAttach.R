.onAttach <- function(libname, pkgname) {
    # Runs when attached to search() path such as by library() or require()
    if (interactive()) {
        packageStartupMessage('v', as.character(packageVersion("gcount")), 
        ', type vignette("gcount-vignette", package="gcount") to get started.')
    }
}
