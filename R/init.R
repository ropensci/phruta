.onAttach <- function(libname, pkgname) {
  message <- c(
    "       ", "\U0002705", "  Thank you for using the {phruta} R package!   ",
    "\U0002705", "\n", "\n                        ", "\U0001F642",
    "Happy coding!!", "\U0001F642", "\n                       ",
    "\U0001F34A", "\U0001F34D", "\U0001F350", "\U0001F95D", "\U0001F965",
    "\U0001F336", "\U0001F951", "\U0001F352", "\U0001F96D"
  )
  packageStartupMessage(paste(message, collapse = ""))
}

pkg.env <- NULL
.onLoad <- function(libname, pkgname) {
  pkg.env <<- new.env()
  assign(".testMode", FALSE, envir = pkg.env)
}

