#' Von Bertalanffy parameter estimation from stock weights
#'
#' Fits a simple VB curve to length-at-age derived from stock weights via
#' length-weight parameters. Returns named coefficients c(linf, k, t0).
#'
#' @param stk An FLStock object
#' @return A named numeric vector with elements linf, k, t0
#' @export
vonB <- function(stk) {
  ln <- FLCore::wt2len(FLCore::stock.wt(stk), FLCore::FLPar(a = 0.1, b = 3))
  ln <- as.data.frame(ln, drop = TRUE)

  vbFn <- function(age, linf, k, t0) linf * (1 - exp(-k * (age - t0)))

  linf_start <- max(ln$data, na.rm = TRUE)
  k_start    <- 0.2

  fit <- stats::nls(
    data ~ vbFn(age, linf, k, t0),
    data  = ln,
    start = list(linf = linf_start, k = k_start, t0 = -0.2)
  )
  stats::coef(fit)
}

#' Age at which an age-specific vector reaches a target (linear interpolation)
#'
#' Canonical interpolation used by \code{\link{mat50}}, \code{\link{sel50}},
#' and \code{\link{findAge}}. Sorts by \code{age}, finds the first bracket where
#' \code{value} crosses \code{target}, and linearly interpolates. Returns
#' \code{NA} if the target is not bracketed (no extrapolation).
#'
#' @param value Numeric vector (e.g. maturity or selectivity).
#' @param age Numeric vector of ages (same length as \code{value}).
#' @param target Level to interpolate (default 0.5).
#' @return Numeric age, or \code{NA_real_}.
#' @export
ageAtTarget <- function(value, age, target = 0.5) {
  ok <- is.finite(value) & is.finite(age)
  value <- value[ok]
  age <- age[ok]
  if (!length(value)) {
    return(NA_real_)
  }
  ord <- order(age)
  age <- age[ord]
  value <- value[ord]
  if (any(abs(value - target) < 1e-8, na.rm = TRUE)) {
    return(age[which.min(abs(value - target))])
  }
  i <- which(value > target)[1]
  if (is.na(i) || i == 1L) {
    return(NA_real_)
  }
  x1 <- age[i - 1L]
  y1 <- value[i - 1L]
  x2 <- age[i]
  y2 <- value[i]
  x1 + (target - y1) * (x2 - x1) / (y2 - y1)
}

#' Age at 50% maturity by year
#'
#' Applies \code{\link{ageAtTarget}} to the \code{mat} slot, one value per year.
#'
#' @param stk An FLStock object.
#' @param target Maturity level (default 0.5).
#' @return Numeric vector named by year.
#' @seealso \code{\link{sel50}}, \code{\link{findAge}}, \code{\link{ageAtTarget}}
#' @export
mat50 <- function(stk, target = 0.5) {
  md <- as.data.frame(FLCore::mat(stk))
  md <- md[is.finite(md$data), , drop = FALSE]
  md$age <- as.numeric(as.character(md$age))
  split_y <- split(md, md$year)
  vapply(
    split_y,
    function(df) ageAtTarget(df$data, df$age, target = target),
    numeric(1)
  )
}

#' Age at 50% selectivity by year
#'
#' Normalises \code{catch.sel} to its maximum each year, then applies
#' \code{\link{ageAtTarget}}.
#'
#' @param stk An FLStock object.
#' @param target Selectivity level (default 0.5).
#' @return Numeric vector named by year.
#' @seealso \code{\link{mat50}}, \code{\link{ageAtTarget}}
#' @export
sel50 <- function(stk, target = 0.5) {
  sel <- as.data.frame(FLCore::catch.sel(stk))
  sel <- sel[is.finite(sel$data), , drop = FALSE]
  sel$age <- as.numeric(as.character(sel$age))
  split_y <- split(sel, sel$year)
  vapply(
    split_y,
    function(df) {
      v <- df$data / max(df$data, na.rm = TRUE)
      ageAtTarget(v, df$age, target = target)
    },
    numeric(1)
  )
}

#' Log-linear mortality ~ weight relationship
#'
#' Regresses log(M) on log(weight) across ages to return intercept and slope.
#'
#' @param stk An FLStock object
#' @return Named numeric vector with elements m1 (intercept) and m2 (slope)
#' @export
lorFn <- function(stk) {
  m  <- as.data.frame(FLCore::m(stk))
  wt <- as.data.frame(FLCore::stock.wt(stk))

  df <- merge(
    m, wt,
    by = c("year", "age", "unit", "season", "area", "iter"),
    suffixes = c(".m", ".w")
  )
  df <- df[is.finite(df$data.m) & is.finite(df$data.w) & df$data.w > 0, ]
  df <- transform(df, logM = log(data.m), logW = log(data.w))

  co <- stats::coef(stats::glm(logM ~ logW, data = df, family = gaussian()))
  names(co) <- c("m1", "m2")
  co
}

#' Transform log(M)~log(W) linear relation to length-inverse mortality params
#'
#' Uses constants c=1.44 and d=1 to convert the linear relation parameters,
#' together with VB parameters linf and k, into the parameterization used
#' by inverse-length mortality models.
#'
#' @param par An FLPar with at least params m1, linf, k
#' @return An FLPar with elements m1, m2, m3, m4
#' @export
lor2gis <- function(par) {
  dimnames(par)$params <- tolower(dimnames(par)$params)
  c <- 1.44
  d <- 1
  b <- par["m1", drop = TRUE]
  a <- par["m1", drop = TRUE] - (c * log(par["linf", drop = TRUE])) - (d * log(par["k", drop = TRUE]))
  FLCore::FLPar(m1 = a, m2 = b, m3 = c, m4 = d)
}
