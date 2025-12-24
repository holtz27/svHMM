#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib holtz, .registration = TRUE
NULL

.onLoad <- function(libname, pkgname) {
  ncores <- parallel::detectCores(logical = FALSE)

  if (is.na(ncores) || ncores < 2) {
    ncores <- 1
  } else {
    ncores <- ncores - 1
  }

  RcppParallel::setThreadOptions(numThreads = ncores)
}

#' Set number of threads used by holtz
#'
#' @param n Number of threads to be used by RcppParallel
#' @export
set_threads <- function(n) {
  stopifnot(is.numeric(n), length(n) == 1, n >= 1)

  max_cores <- tryCatch(
    parallel::detectCores(logical = FALSE),
    error = function(e) NA_integer_
  )

  if (is.na(max_cores)) {
    stop("Unable to detect the number of CPU cores on this system.", call. = FALSE)
  }

  if (n > max_cores) {
    stop(
      sprintf(
        "Invalid number of threads: n = %d. It must satisfy n <= %d (number of physical cores).",
        n, max_cores
      ),
      call. = FALSE
    )
  }

  RcppParallel::setThreadOptions(numThreads = as.integer(n))
  invisible(NULL)
}

