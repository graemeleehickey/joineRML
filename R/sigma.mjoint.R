#' Print \code{mjoint} object
#'
#' @keywords internal
#' @export
sigma.mjoint <- function(object, ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  sig <- sqrt(object$coef$sigma2)
  sig.names <- names(object$formLongFixed)
  if (is.null(sig.names)) {
    names(sig) <- paste0("sigma_", 1:length(sig))
  } else {
    names(sig) <- sig.names
  }

  sig

}
