#' Equilibrium model fitting
#'
#' Fits a stock-recruitment model with \code{ftmb2} (scalar SPR0) and returns an
#' \code{FLBRP} with reference points. Year-varying SPR0 uses \code{eqlFn}
#' (\code{ftmb}). Dynamic-B0 TMB fits are \code{FLRebuild::ftmb_b0dyn}.
#'
#' @param object An \code{FLStock}
#' @param model Stock-recruitment model (default \code{"bevholtSV"})
#' @param ... Passed to methods (\code{nyears}, priors, etc.)
#'
#' @return An \code{FLBRP}
#' @export
setGeneric("eql", function(object, model, ...) standardGeneric("eql"))

#' @rdname eql
#' @param nyears Number of years used for mean SPR0 (default all years)
#' @param prior_s Optional steepness prior mean
#' @param cv_s Optional steepness prior CV
#' @param prior_r0 Optional R0 prior mean
#' @param cv_r0 Optional R0 prior CV
#' @export
setMethod("eql", signature(object = "FLStock"),
          function(object, model = "bevholtSV", nyears = dim(object)[2],
                   prior_s = NULL, cv_s = NULL, prior_r0 = NULL, cv_r0 = NULL) {

            spFn <- function(x) {
              rfs = FLPar(c(ssb.obs(x)),
                          dimnames = list(refpts = "ssb",
                                          quant = dimnames(refpts(x))$quant,
                                          iter = seq(dim(ssb.obs(x))[2])))
              rfs[, -4] = NA
              refpts(x) = rfs
              rtn = data.frame(model.frame(FLQuants(x, ssb = ssb.obs, catch = catch.obs),
                                           drop = TRUE),
                               sp = c(computeRefpts(x)[, "yield"]))
              rtn$pe = (c(rtn$ssb[-1] - rtn$ssb[-dim(rtn)[1]] + rtn$catch[-dim(rtn)[1]] -
                            rtn$sp[-dim(rtn)[1]], NA)) / rtn$ssb
              rtn
            }

            sr   = as.FLSR(object, model = model)
            spr0 = mean(spr0Yr(object)[, dim(object)[2] - (1:nyears) + 1])

            if (model != "segreg") {
              sr = ftmb2(sr, s.est = TRUE,
                         s = 0.7,
                         s.logitsd = 0.4,
                         spr0 = spr0,
                         prior_s = prior_s, cv_s = cv_s,
                         prior_r0 = prior_r0, cv_r0 = cv_r0)
            } else {
              inflect = if ("benchmark" %in% names(attributes(object)))
                benchmark(object)["blim", drop = TRUE] else NA
              sr = ftmb2(sr, s.est = TRUE, inflect = inflect, spr0 = spr0)
            }

            if (is(sr, "FLPar")) {
              sr_params = FLPar(apply(sr, 1, median))
              sr_obj = as.FLSR(object, model = model)
              params(sr_obj) = sr_params
            } else {
              sr_params = FLPar(apply(params(sr), 1, median))
              sr_obj = sr
            }

            rtn = brp(FLBRP(object, nyears = nyears,
                            sr = list(model = do.call(gsub("SV", "", model), list())$model,
                                      params = sr_params)))

            attributes(rtn)[["sr"]]      = sr_obj
            attributes(rtn)[["logLik"]]  = logLik(sr_obj)
            attributes(rtn)[["prod"]]    = spFn(rtn)
            attributes(rtn)[["tseries"]] = tseries(object)
            attributes(rtn)[["eb.obs"]]  = ebiomass(object)

            if ("benchmark" %in% names(attributes(object)))
              attributes(rtn)[["benchmark"]] = benchmark(object)

            rtn
          })

#' Equilibrium FLBRP with median SR parameters from \code{ftmb}
#'
#' @param object An \code{FLStock}
#' @param model Stock-recruitment model (default \code{"bevholtSV"})
#' @return An \code{FLBRP}
#' @export
eqlFn <- function(object, model = "bevholtSV") {
  spr0 = spr0Yr(object)
  sr   = as.FLSR(object, model = model)
  sr   = ftmb(sr, s.est = TRUE,
              s = 0.7,
              s.logitsd = 0.4,
              spr0 = spr0)

  rtn = brp(FLBRP(object, nyears = dim(object)[2],
                  sr = list(model = do.call(gsub("SV", "", model), list())$model,
                            params = FLPar(apply(params(sr), 1, median)))))

  attributes(rtn)[["logLik"]]        = logLik(sr)
  attributes(rtn)[["rec.residuals"]] = residuals(sr)
  attributes(rtn)[["eb.obs"]]        = ebiomass(object)
  rtn
}
