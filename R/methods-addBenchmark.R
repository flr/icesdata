#' Convert a SAG reference-point row to an FLPar benchmark
#'
#' Maps SAG / advice names onto the lowercase params used by
#' \code{\link{benchmark}}: \code{fmsy}, \code{flim}, \code{fpa},
#' \code{blim}, \code{bpa}, \code{btrigger}.
#'
#' @param x A one-row \code{data.frame}, named numeric vector, or \code{FLPar}
#' @return An \code{FLPar} with params \code{fmsy}, \code{flim}, \code{fpa},
#'   \code{blim}, \code{bpa}, \code{btrigger}
#' @export
rfptsToBenchmark <- function(x) {
  if (is(x, "FLPar")) {
    nms = tolower(dimnames(x)$params)
    dimnames(x)$params = nms
    if ("msybtrigger" %in% nms && !"btrigger" %in% nms) {
      x = rbind(x, FLPar(btrigger = c(x["msybtrigger"])))
      nms = dimnames(x)$params
    }
    want = c("fmsy", "flim", "fpa", "blim", "bpa", "btrigger")
    have = intersect(want, nms)
    out = FLPar(NA, dimnames = list(params = want, iter = 1))
    if (length(have)) out[have] = x[have]
    return(out)
  }

  if (is.data.frame(x)) {
    if (!nrow(x)) stop("rfpts data.frame has no rows.", call. = FALSE)
    if (nrow(x) > 1L) x = x[1, , drop = FALSE]
    nms = names(x)
    vals = unlist(x[1, ], use.names = TRUE)
  } else if (is.numeric(x) || is.list(x)) {
    vals = unlist(x)
    nms = names(vals)
  } else {
    stop("x must be a data.frame, named numeric, or FLPar.", call. = FALSE)
  }

  if (is.null(nms) || !any(nzchar(nms)))
    stop("Reference points must be named.", call. = FALSE)

  key = tolower(nms)
  # drop id / year columns
  drop = key %in% c("sid", ".id", "assyear", "assyr", "year", "assessmentyear",
                    "fishstock", "stock")
  vals = vals[!drop]
  key = key[!drop]

  rename = c(
    fmsy = "fmsy", f = "fmsy",
    flim = "flim",
    fpa = "fpa",
    blim = "blim",
    bpa = "bpa",
    btrigger = "btrigger", msybtrigger = "btrigger", btrig = "btrigger",
    bmsy = "bmsy")
  std = rename[key]
  keep = !is.na(std)
  vals = suppressWarnings(as.numeric(vals[keep]))
  std = unname(std[keep])

  out = FLPar(NA, dimnames = list(
    params = c("fmsy", "flim", "fpa", "blim", "bpa", "btrigger"),
    iter = 1))
  for (i in seq_along(std)) {
    p = std[[i]]
    if (p %in% dimnames(out)$params && !is.na(vals[[i]]))
      out[p] = vals[[i]]
  }
  out
}

.pickRfptsRow <- function(rfpts, sid) {
  rfpts = .sagStdKeys(as.data.frame(rfpts, stringsAsFactors = FALSE))
  idCol = if ("sid" %in% names(rfpts)) "sid" else if (".id" %in% names(rfpts)) ".id"
  if (is.null(idCol))
    stop("rfpts must contain a sid column.", call. = FALSE)
  row = rfpts[as.character(rfpts[[idCol]]) == as.character(sid), , drop = FALSE]
  if (!nrow(row))
    stop("No SAG reference points for sid '", sid, "'.", call. = FALSE)
  if ("assYear" %in% names(row) || "assYr" %in% names(row)) {
    yrCol = if ("assYear" %in% names(row)) "assYear" else "assYr"
    y = suppressWarnings(as.integer(row[[yrCol]]))
    row = row[which.max(y), , drop = FALSE]
  } else {
    row = row[1, , drop = FALSE]
  }
  row
}

.fetchBenchmarkRfpts <- function(sid, year = NULL, yearBack = 0L, root = NULL,
                                 quiet = TRUE) {
  sid = as.character(sid)[1]
  if (!nzchar(sid) || is.na(sid))
    stop("sid is missing; set name(object) or pass sid=.", call. = FALSE)

  if (is.null(year)) {
    yearBack = as.integer(yearBack)[1]
    if (is.na(yearBack) || yearBack < 0L)
      stop("yearBack must be a non-negative integer.", call. = FALSE)
    years = as.integer(format(Sys.Date(), "%Y")) - seq.int(0L, yearBack)
  } else {
    years = as.integer(year)[1]
  }

  lastErr = NULL
  for (y in years) {
    rp = tryCatch(
      loadSagRefpts(sid = sid, year = y, root = root),
      error = function(e) {
        lastErr <<- conditionMessage(e)
        NULL
      })
    if (is.null(rp) || !nrow(rp)) next
    row = tryCatch(.pickRfptsRow(rp, sid), error = function(e) NULL)
    if (!is.null(row) && nrow(row)) {
      if (!quiet)
        message("SAG benchmarks for ", sid, " from assessment year ",
                if ("assYear" %in% names(row)) row$assYear else y)
      return(row)
    }
  }
  stop("Could not load SAG benchmarks for '", sid, "' (tried ",
       paste(years, collapse = ", "), ")",
       if (!is.null(lastErr)) paste0(": ", lastErr) else ".",
       call. = FALSE)
}

#' @rdname addBenchmark
#' @param sid Stock id (default \code{name(object)})
#' @param year Assessment year, or \code{NULL} to use current year (and
#'   \code{yearBack} fall-backs)
#' @param yearBack Extra prior years to try when \code{year} is \code{NULL}
#' @param root Optional local \code{sdGraphs} path for \code{\link{loadSagRefpts}}
#' @param rfpts Optional SAG refpts \code{data.frame} (e.g. \code{getSAG()$rfpts});
#'   skips the download when supplied
#' @param quiet Suppress progress messages
#' @export
setMethod("addBenchmark", signature(object = "FLStock"),
          function(object, sid = name(object), year = NULL, yearBack = 0L,
                   root = NULL, rfpts = NULL, quiet = TRUE, ...) {
            if (is.null(rfpts)) {
              row = .fetchBenchmarkRfpts(sid, year = year, yearBack = yearBack,
                                        root = root, quiet = quiet)
            } else {
              row = .pickRfptsRow(rfpts, sid)
            }
            bm = rfptsToBenchmark(row)
            # Store as one-row data.frame — same form as packaged icesdata stocks
            attributes(object)$benchmark = data.frame(
              Fmsy     = c(bm["fmsy"]),
              Flim     = c(bm["flim"]),
              Fpa      = c(bm["fpa"]),
              Blim     = c(bm["blim"]),
              Bpa      = c(bm["bpa"]),
              Btrigger = c(bm["btrigger"]),
              stringsAsFactors = FALSE)
            object
          })

#' @rdname addBenchmark
#' @export
setMethod("addBenchmark", signature(object = "FLStocks"),
          function(object, year = NULL, yearBack = 0L, root = NULL,
                   rfpts = NULL, quiet = TRUE, ...) {
            sids = names(object)
            if (is.null(sids) || !all(nzchar(sids)))
              sids = vapply(object, function(x) as.character(name(x))[1],
                            character(1))

            if (is.null(rfpts)) {
              sag = getSAG(sid = unique(sids), year = year, yearBack = yearBack,
                           root = root, quiet = quiet)
              rfpts = sag$rfpts
            }

            for (i in seq_along(object)) {
              object[[i]] = tryCatch(
                addBenchmark(object[[i]], sid = sids[[i]], rfpts = rfpts,
                             quiet = quiet),
                error = function(e) {
                  warning(sids[[i]], ": ", conditionMessage(e), call. = FALSE)
                  object[[i]]
                })
            }
            object
          })

#' @rdname addBenchmark
#' @param updateRefs If \code{TRUE}, also place Blim / MSYBtrigger / FMSY on
#'   \code{refpts(object)} and recompute
#' @export
setMethod("addBenchmark", signature(object = "FLBRP"),
          function(object, sid = name(object), year = NULL, yearBack = 0L,
                   root = NULL, rfpts = NULL, quiet = TRUE,
                   updateRefs = FALSE, ...) {
            if (is.null(rfpts)) {
              row = .fetchBenchmarkRfpts(sid, year = year, yearBack = yearBack,
                                        root = root, quiet = quiet)
            } else {
              row = .pickRfptsRow(rfpts, sid)
            }
            bm = rfptsToBenchmark(row)
            attributes(object)$benchmark = data.frame(
              Fmsy     = c(bm["fmsy"]),
              Flim     = c(bm["flim"]),
              Fpa      = c(bm["fpa"]),
              Blim     = c(bm["blim"]),
              Bpa      = c(bm["bpa"]),
              Btrigger = c(bm["btrigger"]),
              stringsAsFactors = FALSE)

            if (isTRUE(updateRefs)) {
              # Attach advice points on refpts (same idea as ARrefs)
              rfs = refpts(object)
              need = c("Blim", "MSYBtrigger", "FMSY")
              have = dimnames(rfs)$refpt
              for (nm in setdiff(need, have)) {
                add = rfs[1, ]
                add[] = NA
                dimnames(add)$refpt = nm
                rfs = rbind(rfs, add)
              }
              if (!is.na(c(bm["blim"])))
                rfs["Blim", "ssb"] = c(bm["blim"])
              if (!is.na(c(bm["btrigger"])))
                rfs["MSYBtrigger", "ssb"] = c(bm["btrigger"])
              if (!is.na(c(bm["fmsy"])))
                rfs["FMSY", "harvest"] = c(bm["fmsy"])
              refpts(object) = rfs
              refpts(object) = computeRefpts(object)
            }
            object
          })
