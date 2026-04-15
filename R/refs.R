#' Build Stock Reference Summary Table
#'
#' Returns a wide reference table for a list/FLStocks object, including
#' benchmark, eqsim, fishlife-derived fields and depletion summaries.
#'
#' @param object A list/FLStocks collection of stock objects.
#'
#' @return A data.frame with one row per stock and columns including
#'   `shape`, `psi`, `ssb.minyear`, `ssb.maxyear`, `ssb.minyear.year`,
#'   and `ssb.maxyear.year`.
#'
#' @export
refs <- function(object) {
  eqsm <- eqsim(object)
  benchm <- benchmark(object)

  initial <- plyr::ldply(object, function(x) {
    ts <- tseries(x)
    ts[ts$year == min(ts$year), ]
  })

  current <- plyr::ldply(object, function(x) {
    ts <- tseries(x)
    ts[ts$year == max(ts$year), ]
  })

  fl <- fishlife(object)

  priors <- merge(benchm[, -8], eqsm[, c(".id", "bmsy", "b0")], by=".id")
  priors <- merge(priors, fl[, c(".id", "r")], by=".id")
  priors <- merge(
    priors,
    initial[, c(".id", "ssb", "year")],
    by=".id"
  )
  names(priors)[names(priors) == "ssb"] <- "initial"
  names(priors)[names(priors) == "year"] <- "ssb.minyear.year"

  priors <- merge(
    priors,
    current[, c(".id", "ssb", "year")],
    by=".id"
  )
  names(priors)[names(priors) == "ssb"] <- "current"
  names(priors)[names(priors) == "year"] <- "ssb.maxyear.year"

  transform(
    priors,
    shape = bmsy / b0,
    ssb.maxyear = current / bmsy,
    ssb.minyear = initial / b0,
    psi = initial / b0
  )
}
