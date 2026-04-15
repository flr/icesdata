#' Generate Length Samples from Operating Model
#' 
#' @title lenSamples
#' 
#' @description Generates length frequency samples from an FLStock operating model.
#' Converts catch-at-age to length frequencies using inverse age-length keys.
#' 
#' @param object \code{FLStock} object with catch data
#' @param par \code{FLPar} object with life history parameters (linf, k, l50)
#' @param n Numeric. Sample size for length frequencies
#' @param bin Numeric. Length bin width (default: 0.5)
#' @param extendToLinf Logical. Extend length bins to Linf * 1.1 (default: TRUE)
#' @param ... any additional arguments
#' 
#' @return \code{FLQuant} with length frequency data
#' 
#' @author Based on code from generic.Rmd and case.Rmd
#' 
#' @aliases lenSamples lenSamples-method lenSamples,FLStock-method
#' 
#' @docType methods
#' 
#' @rdname lenSamples
#' 
#' @export lenSamples
#' @exportMethod lenSamples
#' 
#' @examples
#' \dontrun{
#' lsmps = lenSamples(om, par, n = 5000)
#' }
setGeneric('lenSamples', function(object, par, n, ...) standardGeneric('lenSamples'))

setMethod('lenSamples', signature(object = 'FLStock', par = 'FLPar', n = 'numeric'),
  function(object, par, n, bin = 0.5, extendToLinf = TRUE, ...) {
    
    # Convert catch weight to length
    cln = wt2len(catch.wt(object), par)
    cln[is.na(cln)] = 0
    
    # Create inverse age-length key
    iAlk = FLCandy:::invALK(cln, bin = bin)
    iAlk[is.na(iAlk)] = 0
    
    # Ensure catch numbers are not NA
    catch.n(object)[is.na(catch.n(object))] = 0
    
    # Generate length samples
    lsmps = FLCandy:::lenSamp2(catch.n(object), iAlk, n = n)
    
    # Extend length bins to Linf * 1.1 if needed
    if (extendToLinf && max(an(dimnames(lsmps)[[1]])) <= par["linf", drop = TRUE]) {
      newDims = dimnames(lsmps)
      binWidth = an(newDims[[1]][2]) - an(newDims[[1]][1])
      newDims[[1]] = seq(an(newDims[[1]][1]), par["linf", drop = TRUE] * 1.1, binWidth)
      newLsmps = FLQuant(1e-6, dimnames = newDims, units = units(lsmps))
      newLsmps[dimnames(lsmps)[[1]], , , , ] = lsmps
      lsmps = newLsmps
    }
    
    return(lsmps)
  })

