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

#' Age at 50% maturity by year via linear interpolation
#'
#' Computes, per year, the age where maturity crosses 0.5 using linear
#' interpolation between adjacent ages.
#'
#' @param stk An FLStock object
#' @return A numeric vector named by year
#' @export
mat50 <- function(stk) {
  md <- as.data.frame(FLCore::mat(stk))
  md <- md[is.finite(md$data), ]
  md$age <- as.numeric(as.character(md$age))

  fn <- function(df) {
    df <- df[order(df$age), ]
    if (any(abs(df$data - 0.5) < 1e-6))
      return(df$age[which.min(abs(df$data - 0.5))])
    i <- which(df$data > 0.5)[1]
    if (is.na(i) || i == 1) return(NA_real_)
    x1 <- df$age[i - 1]; y1 <- df$data[i - 1]
    x2 <- df$age[i];     y2 <- df$data[i]
    x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)
  }

  split_y <- split(md, md$year)
  sapply(split_y, fn)
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
