#' Dynamic B0 approximation using FLStock and FLSR
#'
#' Approximate Stock Synthesis (SS3)'s dynamic \eqn{B_0} within the FLR
#' framework using \code{FLStock} and \code{FLSR} objects. The function
#' constructs a time series of unfished spawning stock biomass by replaying
#' the stock history under an unfished scenario while preserving recruitment
#' variability and time-varying biology.
#'
#' Specifically, the procedure is:
#' \enumerate{
#' \item Fit a Beverton--Holt stock--recruitment relationship (SRR) to the
#'   assessment \code{FLStock}, using the \code{FLSR} class, with the model
#'   parameterised to constrain \eqn{R_0}. Recruitment deviations are taken
#'   from \code{residuals(FLSR)} and treated as lognormal multiplicative
#'   effects on expected recruitment.
#' \item Specify a recruitment time series as the SR expectation multiplied
#'   by the recruitment deviations. This recruitment series is then used in
#'   forward projections of the \code{FLStock}.
#' \item Run a forward projection of the \code{FLStock} with fishing mortality
#'   set to zero (\eqn{F = 0} in all projection years), while maintaining the
#'   same time-varying biology as in the assessment: weight-at-age,
#'   maturity-at-age, natural mortality, and other biological schedules are
#'   held at their estimated values through time.
#' \item Derive the unfished spawning stock biomass (SSB) at each year from
#'   this \eqn{F = 0} projection using the standard FLR calculations
#'   (i.e., SSB = sum over ages of N-at-age \eqn{\times} weight
#'   \eqn{\times} maturity \eqn{\times} spawning fraction). The resulting
#'   time series is FLR's approximation to SS3's dynamic \eqn{B_0}, i.e. the
#'   SSB the stock would have had in each year if it had experienced the same
#'   recruitment and biology, but without fishing.
#' \item Define depletion as the ratio of the fitted SSB (from the assessment
#'   model) to the dynamic unfished SSB in the corresponding year. The dynamic
#'   trajectory can also be linked to a notion of virgin biomass obtained from
#'   a per-recruit analysis and the Beverton--Holt production function
#'   (following Sissenwine & Shepherd), which provides a reference \eqn{B_0}
#'   under constant recruitment and unfished conditions.
#' }
#'
#' The main differences from SS3's internal implementation of dynamic
#' \eqn{B_0} are structural and relate to the treatment of the initial
#' population:
#' \itemize{
#' \item FLR projections operate on age-by-year matrices for a single
#'   \code{FLStock}, with optional unit/season/area dimensions, but typically
#'   without the full fleet-by-area-by-season structure that SS3 uses. SS3's
#'   dynamic \eqn{B_0} calculations account explicitly for fleet-specific
#'   fishing mortality, spatial structure, and seasonal timing of biology and
#'   fishing.
#' \item In SS3 the initial state is the same for the fished and unfished
#'   populations, and the unfished trajectory is constructed internally from
#'   that state. In the FLR implementation described here, the population
#'   vector (numbers at age in the first modelled year) is taken from the
#'   fitted assessment and then adjusted to represent an unfished state by
#'   inflating numbers at age cohort-by-cohort using the historical fishing
#'   mortality schedule (approximately via \eqn{\exp(\mathrm{cumsum}(F))}
#'   along each cohort's F history).
#' }
#'
#' As a result, dynamic \eqn{B_0} in FLR is conceptually aligned with SS3
#' (unfished SSB through time given the same SRR, recruitment deviations, and
#' time-varying biology) but uses a simpler representation of fleets and areas,
#' and an explicit reconstruction of the initial cohort structure.
#'
#' @param object An \code{FLStock} object representing the assessed stock.
#' @param sr An \code{FLSR} object representing the fitted
#'   Beverton--Holt stock--recruitment relationship.
#' @param ... Additional arguments passed to the forward projection
#'   function used internally (e.g. control objects).
#'
#' @return An \code{FLStock} object containing the dynamic \eqn{B_0}
#'   trajectory (unfished SSB through time) under the specified SRR and
#'   recruitment deviations, with the same biological schedules as
#'   \code{object}. Depletion and related reference points can be computed
#'   by comparing \code{ssb(object)} and \code{ssb(b0dyn(object, sr))}.
#'
#' @seealso \code{\link{FLSR}}, \code{\link{ssb}}, and SS3 documentation
#'   on dynamic \eqn{B_0}.
#'
#' @examples
#' \dontrun{
#'   # Fit Beverton-Holt SRR
#'   sr0 <- fmle(as.FLSR(stk, model = "bevholtSV"))
#'
#'   # Compute dynamic B0 trajectory
#'   stk.b0 <- b0dyn(stk, sr0)
#'
#'   # Depletion relative to dynamic B0
#'   dep <- ssb(stk) / ssb(stk.b0)
#' }
#'
#' @export
#' 
#' 
# ab  {{{
setMethod('ab', signature(x='FLSR', model='missing'),
          function(x)
          {
            res <- x
            model(res) <- sub('SV', '', SRModelName(model(x)))
            params(res) <- ab(params(x), SRModelName(model(x)))
            residuals(res) <- residuals(x)
            fitted(res)    <- fitted(x)
            ssb(res)    <- ssb(x)
            
            return(res)
          }
) # }}}

dims.=function(x) unlist(dims(x))

setGeneric("inRec", function(object, nyr, ...) standardGeneric("inRec"))

setMethod("inRec", signature(object="FLStock", nyr="missing"), function(object, mort=FLCore::z) {
            
      # Cumulative survival from age 0 to current age in first year
      cumZ=as.data.frame(exp(-apply(mort(object)[, 1], 2:6, cumsum)))
            
      # Shift ages by +1 and invert survival to get per‑age recruitment Z
      cumZ=rbind(
              cumZ[1, ],
              transform(head(cumZ, -1),
                        age  = age + 1,
                        data = 1 / data)
            )
            cumZ[1, "data"]=1
            names(cumZ)[7]="Z"
            
            # Numbers at age in first year
            intlN=as.data.frame(stock.n(object)[, 1])
            names(intlN)[7]="N"
            
            # Reconstructed recruitment R = Z * N
            rtn=transform(merge(cumZ, intlN), R = Z * N)
            
            # Put recruitment at year = year ‑ age, age = min(age)
            rtn=as.FLQuant(
              transform(rtn,
                        year = year - age,
                        age  = min(age),
                        data = R)[, -(7:9)]
            )
            
            # Expand rec window backwards to minyear - max age
            rec.ext=window(rec(object),
                              start = an(dims(object)["minyear"]) -
                                an(dims(object)["max"]))
            
            rec.ext[, ac((an(dims(object)["minyear"]) - an(dims(object)["max"])):
                           (an(dims(object)["minyear"])))]=rtn
            
            rec.ext
          })

setMethod("inRec", signature(object = "FLStock", nyr = "numeric"),
          function(object, nyr) {
            
            # Extended rec window from minyear - maxage
            rec.ext <- window(rec(object),
                              start = an(dims(object)["minyear"]) -
                                an(dims(object)["max"]))
            
            # Typical recruitment (geometric mean) over first nyr years
            # apply over iterations etc, dimensions 3:6
            typical_rec <- apply(rec(object)[, seq(nyr)], 3:6,
                                 function(x) exp(mean(log(x))))
            
            # Fill pre‑assessment years (minyear - maxage : minyear - 1)
            rec.ext[, ac((an(dims(object)["minyear"]) - an(dims(object)["max"])):
                           (an(dims(object)["minyear"]) - 1))] <- typical_rec
            
            rec.ext
          }
)

setGeneric("survivors", function(object,sr,...) standardGeneric("survivors"))

#' Survivorship from total mortality
#'
#' @param object \code{FLQuant} of total mortality, \code{FLStock}, or \code{FLBRP}.
#' @param cohort Logical; cohort-based survivorship for stock/BRP objects.
#' @param ... Additional arguments.
#' @export

## FLQuant method: at end of the year
setMethod("survivors", signature(object="FLQuant", sr="missing"), 
          function(object, ...) {
  
  minyr =an(dims(object)["minyear"])
  maxage=an(dims(object)["max"])
            
  mrt    =object
  mrt[1] =0
  mrt[-1]=object[-dim(object)[1]]

  mrt    =window(mrt,start=minyr-maxage)
  mrt[,seq(maxage)]=mrt[,maxage+1]

  mrt=FLCohort(mrt)
  
  rtn=exp(-apply(mrt,c(2,3,5,6),function(object) {
            flag=is.na(object)
            
            if (all(flag)) 
              return(object)
            
            object[!flag]=cumsum(object[!flag])
            object}))
  
  rtn=rtn[, ac((an(dims(object)["minyear"])-an(dims(object)["max"])):an(dims(object)["maxyear"]))]
  
  return(rtn)})


setGeneric("b0dyn", function(object,sr,...) standardGeneric("b0dyn"))

#' Survivorship from total mortality
#'
#' @param object \code{FLStock}.
#' @param object \code{FLSR}.
#' @param ... Additional arguments.
#' @export

setMethod("b0dyn", signature(object="FLStock", sr="FLSR"), function(object,sr,residuals=1,...) {
            
            yrs       =seq(max(1,dims(object)$min))
            srvAdj    =exp(apply(apply(harvest(object)[,yrs],c(1,3:6),mean),3:6,cumsum))
            srvAdj[-1]=srvAdj[-dim(srvAdj)[1]]
            srvAdj[1] =1
            stock.n(object)[,yrs]=stock.n(object)[,yrs]+stock.n(object)[,yrs]%*%srvAdj
            
            catch.n(   object)[,yrs]=0
            landings.n(object)[,yrs]=0
            discards.n(object)[,yrs]=0
            harvest(   object)[,yrs]=0
            
            rsdl=window(exp(residuals(sr)*residuals),start=dims(object)$minyear)
            
            # Fbar set to 0 for projection years
            stk0=FLCore::ffwd(object,fbar=fbar(object)[,-yrs] %=% 0,
                                     sr=sr, deviances=rsdl)
            
            stk0})

