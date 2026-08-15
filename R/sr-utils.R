#' Convert Logit-Transformed Steepness to Original Scale
#'
#' @param logit_h Numeric vector of logit-transformed steepness values
#' @return Numeric vector of steepness values on original scale (0.2 to 1.0)
#' @export
from_logits <- function(logit_h) {
  if (!is.numeric(logit_h)) stop("logit_h must be numeric")
  0.2001 + 0.7998 * 1 / (1 + exp(-logit_h))
}

#' Convert Steepness to Logit Scale
#'
#' @param h Numeric vector of steepness values (0.2 to 1.0)
#' @return Numeric vector of logit-transformed steepness values
#' @export
to_logits <- function(h) {
  if (!is.numeric(h)) stop("h must be numeric")
  if (any(h <= 0.2001 | h >= 1.0)) stop("h must be between 0.2001 and 1.0")
  -log(0.7998 / (h - 0.2001) - 1)
}
