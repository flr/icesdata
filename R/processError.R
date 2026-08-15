#' Surplus production along the equilibrium production curve
#'
#' @param stk An \code{FLStock}
#' @param eq An \code{FLBRP}
#' @param stock Function returning the stock metric (default \code{FLCore::ssb})
#' @return An \code{FLQuant} of surplus production
#' @export
sp <- function(stk, eq, stock = FLCore::ssb) {
  fbar(eq) = FLQuant(seq(0, 1, length.out = 201)) * computeRefpts(eq)["crash", "harvest"]
  mf = model.frame(FLQuants(eq,
                            "stock" = function(x) stock(x),
                            "catch" = function(x) catch(x)), drop = TRUE)
  dat = with(mf, approx(stock, catch, xout = c(stock(stk))))
  FLQuant(dat$y, dimnames = dimnames(stock(stk)))
}

#' Process error from stock and equilibrium production
#'
#' @param stk An \code{FLStock}
#' @param eq An \code{FLBRP}
#' @param stock Function returning the stock metric (default \code{FLCore::ssb})
#' @return An \code{FLQuant} of process error
#' @export
pe <- function(stk, eq, stock = FLCore::ssb) {
  (stock(stk) %-%
     window(stock(stk)[, -1], end = dims(stk)$maxyear + 1) -
     catch(stk) %+% sp(stk, eq, stock)) %/% stock(stk)
}

#' @rdname processError
#' @export
setMethod("processError", signature(object = "FLBRP"),
          function(object) {
            rfs = FLPar(c(ssb.obs(object)),
                        dimnames = list(refpts = "ssb",
                                        quant = dimnames(refpts(object))$quant,
                                        iter = seq(dim(ssb.obs(object))[2])))
            rfs[, -4] = NA
            refpts(object) = rfs

            rtn = data.frame(model.frame(FLQuants(object,
                                                  ssb = ssb.obs,
                                                  catch = catch.obs),
                                         drop = TRUE),
                             sp = c(computeRefpts(object)[, "yield"]))

            rtn$pe = (c(rtn$ssb[-1] - rtn$ssb[-dim(rtn)[1]] +
                          rtn$catch[-dim(rtn)[1]] -
                          rtn$sp[-dim(rtn)[1]], NA)) / rtn$ssb

            FLQuants(
              ssb   = as.FLQuant(data.frame(year = rtn$year, data = rtn$ssb)),
              catch = as.FLQuant(data.frame(year = rtn$year, data = rtn$catch)),
              pe    = as.FLQuant(data.frame(year = rtn$year, data = rtn$pe)),
              sp    = as.FLQuant(data.frame(year = rtn$year, data = rtn$sp))
            )
          })

#' Process error with year-varying biology from an FLStock and SR
#'
#' @param object An \code{FLStock}
#' @param sr An \code{FLSR} or \code{FLBRP}
#' @param slots Character vector of slots to vary by year
#' @param log Unused; retained for compatibility
#' @return An \code{FLQuants} with ssb, catch, production and error
#' @export
processErrorStock <- function(object, sr,
                              slots = c("m", "mat", "stock.wt", "catch.wt", "catch.sel"),
                              log = FALSE) {

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

  if ("catch.sel" %in% slots | "landings.sel" %in% slots)
    landings.sel(eq) = year2iter(catch.sel(object) %*% landings.n(object) %/% catch.n(object))
  if ("catch.sel" %in% slots | "discards.sel" %in% slots)
    discards.sel(eq) = year2iter(catch.sel(object) %*% discards.n(object) %/% catch.n(object))

  nms = dimnames(refpts(eq))
  nms[[1]] = c(nms[[1]], "spr.100")
  refpts(eq) = FLPar(array(NA, lapply(nms, length), nms))

  model(eq)  = model(sr)
  params(eq) = params(sr)

  dmns = list(refpt = c("production"),
              quant = c("harvest", "yield", "rec", "ssb", "biomass", "revenue", "cost", "profit"),
              iter = dimnames(eq)$iter)
  prd = FLPar(array(NA, dim = vapply(dmns, length, 1L), dimnames = dmns))
  prd[, "ssb"] = c(ssb(object))
  refpts(eq) = prd
  prd = computeRefpts(eq)

  production = as.FLQuant(c(prd[, "yield"]), dimnames = dimnames(ssb(object)))
  ssb.t = window(ssb(object),
                 start = dims(ssb(object))$minyear + 1,
                 end = dims(ssb(object))$maxyear + 1)

  FLQuants(ssb = ssb(object),
           catch = catch(object),
           production = production,
           error = (1 / ssb(object)) * (ssb.t - catch(object) + production))
}
