#' Non-stationary reference points
#'
#' Sets M, wt, sel and mat to vary by iter based on annual values to examine
#' non-stationarity in reference points.
#'
#' @param object An \code{FLStock} object
#' @param sr An \code{FLSR} or \code{FLBRP} providing stock-recruit model and params
#' @param slots Character vector of biology/selectivity slots to vary by year
#' @param matSel Logical; if TRUE set selectivity equal to maturity
#' @param abi Logical; if TRUE return only the ABI time series
#'
#' @return A data.frame of reference points by year, or an FLQuant if \code{abi=TRUE}
#' @export
nonStationarity <- function(object, sr,
                            slots = c("m", "mat", "stock.wt", "catch.wt", "catch.sel"),
                            matSel = FALSE, abi = FALSE) {

  eq = FLBRP(object)
  eq = propagate(eq, dim(object)[2])

  year2iter <- function(x) {
    tmp = as.data.frame(x, drop = TRUE)
    names(tmp)[names(tmp) == "year"] = "iter"
    as.FLQuant(tmp)
  }

  if ("m" %in% slots)     m(eq) = year2iter(m(object))
  if ("mat" %in% slots) mat(eq) = year2iter(mat(object))

  if ("stock.wt" %in% slots) stock.wt(eq) = year2iter(stock.wt(object))
  if ("catch.wt" %in% slots | "landings.wt" %in% slots) landings.wt(eq) = year2iter(landings.wt(object))
  if ("catch.wt" %in% slots | "discards.wt" %in% slots) discards.wt(eq) = year2iter(discards.wt(object))

  if ("catch.sel" %in% slots | "landings.sel" %in% slots) {
    sel = catch.sel(object) %*% landings.n(object) %/% catch.n(object)
    sel[is.na(sel)] = 0
    sel[!is.finite(sel)] = 0
    landings.sel(eq) = year2iter(sel)
  }

  if ("catch.sel" %in% slots | "discards.sel" %in% slots) {
    sel = catch.sel(object) %*% discards.n(object) %/% catch.n(object)
    sel[is.na(sel)] = 0
    sel[!is.finite(sel)] = 0
    discards.sel(eq) = year2iter(sel)
  }

  if (matSel) {
    landings.sel(eq) = mat(eq)
    discards.sel(eq)[] = 0
  }

  nms = dimnames(refpts(eq))
  nms[[1]] = c(nms[[1]], "spr.100")
  refpts(eq) = FLPar(array(NA, lapply(nms, length), nms))

  model(eq)  = model(sr)
  params(eq) = params(sr)

  if (abi)
    return(FLQuant(c(abiMsy(eq)), dimnames = dimnames(fbar(object))))

  refpts(eq) = computeRefpts(eq)
  refpts(eq) = rbind(refpts(eq), refpts(eq)[1, ])
  dimnames(refpts(eq))$refpt[dim(refpts(eq))[1]] = "current"
  refpts(eq)["current"] = NA
  refpts(eq)["current", "ssb"] = c(FLCore::iter(ssb.obs(eq), 1))

  rtn = rbind(refpts(eq), FLBRP::properties(eq))
  rtn = rtn[!duplicated(dimnames(rtn)[[1]])]
  rtn = rtn[, apply(rtn, 2, function(x) all(is.na(x))) == 0]
  names(dimnames(rtn)) = names(refpts(eq))

  ebio = computeRefpts(eq)[, "harvest", drop = TRUE]
  ebio = FLPar(c(ebio), dimnames = dimnames(rtn[, 1]))
  dimnames(ebio)[[2]] = "eb"

  abiVal = FLQuant(c(abiMsy(eq)), dimnames = dimnames(fbar(object)))
  abiVal = FLPar(c(abiVal), dimnames = dimnames(rtn[, 1]))
  dimnames(abiVal)[[2]] = "abi"

  spr0 = refpts(eq)["virgin", "ssb"] %/% refpts(eq)["virgin", "rec"]
  dimnames(spr0)[[2]] = "spr0"

  yrs = as.numeric(dimnames(object)$year[seq_len(dim(object)[2])])
  rtn = rbind(
    transform(as.data.frame(rtn),    year = yrs[iter])[, -3],
    transform(as.data.frame(spr0),   year = yrs[iter])[, -3],
    transform(as.data.frame(ebio),   year = yrs[iter])[, -3],
    transform(as.data.frame(abiVal), year = yrs[iter])[, -3]
  )

  rtn[do.call("order", rtn[, c(1, 2, 4)]), ]
}
