#' Leslie matrix demographic properties
#'
#' @param x An \code{FLBRP}
#' @param fbar Fishing mortality; default is MSY harvest
#' @return Named numeric vector: damp, r, gt, r0
#' @export
leslieFn <- function(x, fbar = FLQuant(c(refpts(x)["msy", "harvest"]))) {
  L = leslie(x, fbar = fbar)
  if (length(dim(L)) == 3)
    L = L[, , 1, drop = TRUE]
  else
    L = L[drop = TRUE]

  egnMod = Mod(eigen(L)$values)
  sorted = sort(egnMod, decreasing = TRUE)
  dampingRatio = if (length(sorted) > 1) sorted[1] / sorted[2] else NA_real_

  lam = max(egnMod)
  r = log(lam)

  fecundity = L[1, ]
  nages = length(fecundity)
  lx = numeric(nages)
  lx[1] = 1
  for (i in 2:nages)
    lx[i] = lx[i - 1] * L[i, i - 1]

  R0 = sum(lx * fecundity)
  genTime = sum((1:nages) * lx * fecundity) / R0

  c(damp = dampingRatio, r = r, gt = genTime, r0 = R0)
}

ssbN <- function(object) {
  f = harvest(object) %*% harvest.spwn(object)
  m = m(object) %*% m.spwn(object)
  expZ = exp(-(f %+% m))
  stock.n(object) %*% expZ %*% stock.wt(object) %*% mat(object)
}

#' @rdname covarFn
#' @param fbar An \code{FLQuant} fishing mortality (default MSY harvest)
#' @param model SR model name for steepness (\code{"bevholt"} or \code{"ricker"})
#' @export
setMethod("covarFn", signature(object = "FLBRP"),
          function(object,
                   fbar = as.FLQuant(refpts(object)["msy", "harvest", drop = TRUE]),
                   model = "bevholt") {

            smry <- function(x) {
              sel = catch.sel(x) / max(catch.sel(x))
              ecdf1 = ecdf(c(sel))
              ecdf2 = ecdf(c(mat(x)))
              rangeVals = seq(min(c(sel, mat(x))), max(c(sel, mat(x))), length.out = 1000)
              matsel = sum(ecdf1(rangeVals) - ecdf2(rangeVals)) * diff(rangeVals)[1]

              abiSel = as.numeric(c(ages(catch.sel(x))) >= c(abiAge(x)))
              ssbSel = c(mat(x) %*% stock.wt(x))
              ssbSel = ssbSel / max(ssbSel)
              rangeVals = seq(0, 1, length.out = 100)
              abiSSB = sum(ecdf(abiSel)(rangeVals) - ecdf(ssbSel)(rangeVals)) *
                diff(rangeVals)[1]

              peVal = tryIt({
                peQ = processError(x)$pe
                as.numeric(sqrt(var(c(peQ), na.rm = TRUE)))
              })
              if (is.null(peVal)) peVal = NA_real_

              sVal = tryIt(c(sv(params(x), spr0 = spr0(x), model = model)["s"]))
              if (is.null(sVal)) sVal = NA_real_

              rtn = c(
                m.spwn       = any(m.spwn(x) != 0),
                harvest.spwn = any(harvest.spwn(x) != 0),
                wts          = mean(catch.sel(x) %*% (catch.wt(x) - stock.wt(x)) / stock.wt(x)),
                ks           = unname(ks.test(c(sel), c(mat(x)))[[1]]),
                matsel       = matsel,
                abiSSB       = abiSSB,
                spr0         = c(spr0(x)),
                shape        = c(refpts(x)["msy", "ssb"] / refpts(x)["virgin", "ssb"]),
                s            = sVal,
                pe           = peVal
              )
              c(rtn, leslieFn(x))
            }

            fbar(object) = fbar

            apexCtc = catch.n(object) %*% catch.wt(object)
            apexCtc[] = FLCore::fapex(apexCtc)
            apexCtc = apexCtc == catch.n(object) %*% catch.wt(object)
            apexCtc = subset(as.data.frame(as.FLQuant(apexCtc)), data == 1)[, -7]

            apexSSB = ssbN(object)
            apexSSB[] = FLCore::fapex(ssbN(object))
            apexSSB = apexSSB == ssbN(object)
            apexSSB = subset(as.data.frame(as.FLQuant(apexSSB)), data == 1)[, -7]

            names(apexCtc)[1] = "ctcAge"
            names(apexSSB)[1] = "ssbAge"

            mnCwt = as.data.frame(
              quantSums(catch.n(object) %*% catch.wt(object)) %/%
                quantSums(catch.n(object)))[, -1]
            mnSwt = as.data.frame(
              quantSums(stock.n(object) %*% stock.wt(object) %*% mat(object)) %/%
                quantSums(stock.n(object) %*% mat(object)))[, -1]

            names(mnSwt)[6] = "swt"
            names(mnCwt)[6] = "cwt"

            rtn = unlist(merge(
              merge(mnSwt, mnCwt, by = names(mnSwt)[1:5]),
              merge(apexSSB, apexCtc, by = names(apexSSB)[2:6]),
              by = names(mnSwt)[1:5])[, 6:9])

            c(rtn, smry(object))
          })

#' @rdname covarFn
#' @export
setMethod("covarFn", signature(object = "FLBRPs"),
          function(object, ...) {
            plyr::ldply(object, function(x) {
              res = tryIt(covarFn(x, ...))
              if (is.null(res)) return(NULL)
              as.data.frame(as.list(res))
            })
          })

#' @rdname leslie
#' @export
setMethod("leslie", signature(object = "FLBRP", fec = "missing"),
          function(object, fbar = refpts(object)["crash", "harvest"]) {
            fbar(object) = as.FLQuant(fbar[drop = TRUE])

            survivors = exp(-(m(object) %+% harvest(object)))
            zspawn = (harvest(object) %*% harvest.spwn(object)) %+%
              (m(object) %*% m.spwn(object))
            fec = stock.n(object) %*% exp(-zspawn) %*%
              stock.wt(object) %*% mat(object)

            mkL = getFromNamespace(".leslie", "FLCore")
            L = array(0, c(dim(fec)[1], dim(fec)[1], dim(fec)[6]))
            for (i in seq(dims(object)$iter))
              L[, , i] = mkL(iter(survivors, i)[drop = TRUE],
                             iter(fec, i)[drop = TRUE])
            L
          })
