.onLoad <- function(libname, pkgname) {
    gcount_verbose = "FALSE"
    # global options
    opts = c(
              "gcount_verbose" = gcount_verbose
            )
    for (i in setdiff(names(opts), names(options())) ) {
        text = paste('options(', i, '=', opts[i], ')', sep="")
        eval(parse(text=text))
    }
    invisible()
}
