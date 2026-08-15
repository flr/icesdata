#' SAG time series, reference points, and advice helpers for stock-group Rmds.

defaultProjectRoot <- function() {
  if (exists("projectRoot", inherits = TRUE) &&
      is.character(projectRoot) && length(projectRoot) == 1 && nzchar(projectRoot))
    return(projectRoot)
  if (exists("findProjectRoot", mode = "function"))
    return(findProjectRoot())
  envRoot = Sys.getenv("CREST_ROOT", Sys.getenv("RESILIENCE_ROOT", ""))
  if (nzchar(envRoot))
    return(normalizePath(envRoot, winslash = "/", mustWork = TRUE))
  stop("Cannot resolve project root. Set CREST_ROOT or attach CREST.",
       call. = FALSE)
}

#' Locate a StockAssessmentGraphs CSV for an ICES year folder.
#'
#' Picks the lexicographically last \code{StockAssessmentGraphs_<digits>...csv}
#' (excludes Yieldrecruit companions). Keep one Graphs CSV per year folder.
#'
#' @param year Folder under \code{sdGraphs} (e.g. 2025).
#' @param root Path to the \code{sdGraphs} directory. Required; if \code{NULL}
#'   or empty, returns \code{NULL} (no search).
#' @return Full path to one CSV file, or \code{NULL} if none exists.
#' @export
findSagGraphsCsv <- function(year = 2025, root = NULL) {
  if (is.null(root) || !nzchar(root)) return(NULL)
  d = file.path(root, as.character(year))
  files = list.files(d, pattern = "^StockAssessmentGraphs_[0-9].*\\.csv$")
  if (!length(files)) return(NULL)
  file.path(d, sort(files)[length(files)])
}

#' Fetch ICES SAG reference points from the web API.
#'
#' For each stock, tries \code{year} first, then \code{year - 1} if no
#' assessment key or reference points are available.
#'
#' @param sid Character vector of FishStock codes.
#' @param year Preferred assessment year (default: current calendar year).
#' @return Data frame with \code{sid}, \code{assYear}, \code{Blim},
#'   \code{MSYBtrigger}, \code{FMSY}, \code{Fmanagement}, \code{Fpa}.
#' @export
fetchSagRefpts <- function(sid, year = as.integer(format(Sys.Date(), "%Y"))) {
  years = unique(as.integer(c(year, year - 1L)))
  rows = lapply(sid, function(stock) {
    for (y in years) {
      row = tryCatch({
        keys = icesSAG::findAssessmentKey(stock, year = y)
        if (!length(keys) || all(is.na(keys))) {
          NULL
        } else {
          rp = icesSAG::getFishStockReferencePoints(keys[[1]])
          if (!is.data.frame(rp) || !nrow(rp)) {
            NULL
          } else {
            data.frame(
              assYear     = rp$AssessmentYear[[1]],
              sid         = stock,
              Blim        = rp$Blim[[1]],
              Bpa         = if ("Bpa" %in% names(rp)) rp$Bpa[[1]] else NA_real_,
              MSYBtrigger = rp$MSYBtrigger[[1]],
              FMSY        = rp$FMSY[[1]],
              Flim        = if ("Flim" %in% names(rp)) rp$Flim[[1]] else NA_real_,
              Fmanagement = rp$Fmanagement[[1]],
              Fpa         = rp$Fpa[[1]],
              stringsAsFactors = FALSE)
          }
        }
      }, error = function(e) NULL)
      if (!is.null(row)) return(row)
    }
    warning("SAG refpts skipped for ", stock, call. = FALSE)
    NULL
  })
  out = do.call(rbind, rows)
  if (is.null(out) || !nrow(out))
    stop("No SAG reference points returned from the web API.", call. = FALSE)
  rownames(out) = NULL
  out
}

#' Load ICES SAG reference points for a set of stock ids.
#'
#' By default fetches from the ICES SAG web API (\code{\link{fetchSagRefpts}}).
#' Pass \code{root} to read a local StockAssessmentGraphs CSV first (optional
#' offline mirror); \code{root} must be the full path to the \code{sdGraphs}
#' folder — nothing is inferred from the project layout.
#'
#' @param sid Character vector of FishStock codes.
#' @param year Preferred sdGraphs / assessment year (default: current calendar
#'   year).
#' @param root Optional path to a local \code{sdGraphs} directory. \code{NULL}
#'   (default) uses the web API only.
#' @return Data frame with \code{sid} and ref-point columns.
#' @export
loadSagRefpts <- function(sid, year = as.integer(format(Sys.Date(), "%Y")),
                          root = NULL) {
  if (!is.null(root) && nzchar(root)) {
    years = unique(as.integer(c(year, year - 1L)))
    for (y in years) {
      csv = findSagGraphsCsv(year = y, root = root)
      if (!is.null(csv) && file.exists(csv)) {
        rfpts = read.csv(csv, stringsAsFactors = FALSE)
        # AssessmentYear / FishStock are external SAG CSV column names
        keep = c("AssessmentYear", "FishStock", "Blim", "Bpa", "MSYBtrigger",
                 "FMSY", "Flim", "Fmanagement", "Fpa")
        missing = setdiff(c("AssessmentYear", "FishStock", "Blim", "MSYBtrigger",
                            "FMSY", "Fpa"), names(rfpts))
        if (length(missing))
          stop("SAG graphs CSV missing columns: ", paste(missing, collapse = ", "),
               call. = FALSE)
        keep = intersect(keep, names(rfpts))
        rfpts = rfpts[, keep, drop = FALSE]
        names(rfpts)[names(rfpts) == "FishStock"] = "sid"
        names(rfpts)[names(rfpts) == "AssessmentYear"] = "assYear"
        rfpts = subset(rfpts, sid %in% sid)
        rfpts = rfpts[!duplicated(rfpts$sid), ]
        if (nrow(rfpts)) return(rfpts)
      }
    }
  }
  fetchSagRefpts(sid = sid, year = year)
}

#' Read advice.csv into an FLQuants keyed by sid.
#'
#' Prefers local \code{data/advice/advice.csv}, else the packaged copy under
#' \code{inst/extdata/advice.csv}. When \code{sid} is supplied, only matching
#' rows are loaded (never silently returns unmatched stocks).
#'
#' @param adviceFile Path to advice CSV (wide year columns).
#' @param sid Optional character vector of stock ids to keep. \code{NULL}
#'   loads every row in the CSV.
#' @return \code{FLQuants} of advice catch by stock (empty if none match).
#' @export
loadAdviceFlqs <- function(adviceFile = NULL, sid = NULL) {
  if (is.null(adviceFile))
    adviceFile = shippedAdvicePath(defaultProjectRoot())
  wide = read.csv(adviceFile, stringsAsFactors = FALSE)
  if (!is.null(sid)) {
    sid = unique(as.character(sid))
    wide = wide[wide$sid %in% sid, , drop = FALSE]
  }
  if (!nrow(wide))
    return(FLCore::FLQuants())
  advice = reshape::melt(wide)
  FLCore::FLQuants(plyr::dlply(advice, "sid", function(d) {
    FLCore::as.FLQuant(data.frame(
      year = FLCore::an(substr(as.character(d$variable), 2, 5)),
      data = d$value))
  }))
}

#' Pull ICES SAG summary tables and standardise column names.
#'
#' Uses named \code{icesSAG::getSAG} fields (not positional indices).
#' When multiple \code{sagYears} are supplied, each stock keeps rows from its
#' newest assessment year that still has usable data (after catch filtering).
#' Prefer \code{\link{getSAG}} when a per-sid year cascade with early exit is
#' needed: that tries years one stock at a time rather than fetching every
#' stock\times year combination up front.
#'
#' @param sid Stock ids.
#' @param sagYears Assessment year(s). With more than one value, newest usable
#'   year is kept per stock.
#' @param includeTB Include total biomass (\code{TBiomass} as \code{stock}).
#' @param dropAssYearRows Drop rows where data year equals assessment year.
#' @param requireCatches Keep only assessment pulls with non-missing catches
#'   (pelagic/demersal). Set \code{FALSE} for Nephrops.
#' @return Data frame with FLR-style columns: \code{sid}, \code{assYear},
#'   \code{year}, \code{rec}, \code{ssb}, \code{fbar}, \code{catch},
#'   \code{landings}, \code{discards}, and \code{stock} when \code{includeTB}.
#' @export
fetchSagTs <- function(sid,
                       sagYears = 2025,
                       includeTB = FALSE,
                       dropAssYearRows = FALSE,
                       requireCatches = TRUE) {
  grid = expand.grid(stock = sid, year = sagYears, stringsAsFactors = FALSE)
  stdG = plyr::mdply(grid, function(stock, year) {
    tryCatch({
      res = icesSAG::getSAG(stock = stock, year = year)
      # Non-data.frame means "not available" for that year; skip.
      if (!is.data.frame(res)) return(NULL)
      res
    }, error = function(e) {
      warning("SAG series skipped for ", stock, " (", year, "): ",
              conditionMessage(e), call. = FALSE)
      NULL
    })
  })

  if (!is.data.frame(stdG) || !nrow(stdG))
    stop("No SAG time series returned for the requested sid(s).", call. = FALSE)

  # AssessmentYear is the external icesSAG field; mapped to assYear below
  need = c("Year", "recruitment", "SSB", "F", "catches", "landings",
           "discards", "AssessmentYear")
  missing = setdiff(need, names(stdG))
  if (length(missing))
    stop("icesSAG::getSAG result missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)

  stdG = plyr::ddply(stdG, "stock", function(d) {
    if (requireCatches) d = subset(d, !is.na(catches))
    if (!nrow(d)) return(d)
    subset(d, year == max(year))
  })
  if (!nrow(stdG))
    stop("No SAG time series left after filtering.", call. = FALSE)

  # Names follow FLCore accessors: rec, ssb, fbar, catch, landings, discards, stock
  out = data.frame(
    sid      = stdG$stock,
    assYear  = stdG$AssessmentYear,
    year     = stdG$Year,
    rec      = stdG$recruitment,
    ssb      = stdG$SSB,
    fbar     = stdG$F,
    catch    = stdG$catches,
    landings = stdG$landings,
    discards = stdG$discards,
    stringsAsFactors = FALSE)
  if (includeTB && "TBiomass" %in% names(stdG))
    out$stock = stdG$TBiomass
  if (dropAssYearRows)
    out = subset(out, year != assYear)
  out
}

# Fetch one stock for one assessment year; NULL if unavailable / unusable.
fetchSagTsOne <- function(stock, sagYear,
                          includeTB = FALSE,
                          dropAssYearRows = FALSE,
                          requireCatches = TRUE) {
  tryCatch(
    fetchSagTs(
      sid = stock,
      sagYears = sagYear,
      includeTB = includeTB,
      dropAssYearRows = dropAssYearRows,
      requireCatches = requireCatches),
    error = function(e) NULL)
}

# Row-bind data frames that may differ in columns (e.g. TB for Nephrops).
rbindFill <- function(dfs) {
  dfs = Filter(function(d) is.data.frame(d) && nrow(d), dfs)
  if (!length(dfs)) return(NULL)
  cols = unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs = lapply(dfs, function(d) {
    for (m in setdiff(cols, names(d))) d[[m]] = NA
    d[cols]
  })
  out = do.call(rbind, dfs)
  rownames(out) = NULL
  out
}

#' Load SAG time series and reference points for given stock ids.
#'
#' Downloads ICES Stock Assessment Graphs (SAG) summaries and reference points
#' via \code{icesSAG}. Duplicate \code{sid} values are ignored
#' (\code{unique}). Distinct from \code{icesSAG::getSAG}, which returns a
#' single raw assessment table.
#'
#' Does \strong{not} load the local catch-bridge \code{advice.csv}; use
#' \code{\link{loadAdviceFlqs}} (or the case-study script
#' \code{case-studies/bim/02_advice.R}) for that. Does \strong{not} read a
#' cached \code{bim-sag.RData}; that file is written by
#' \code{case-studies/bim/01_retrieve_sag.R} for offline reuse via
#' \code{load()}, separately from this function.
#'
#' \strong{Assessment year (per sid).} When \code{year = NULL} (the default),
#' each requested \code{sid} independently tries the current calendar year,
#' then current year \eqn{- 1}, \ldots, down to current year
#' \eqn{-} \code{yearBack}, and keeps the first year that yields usable time
#' series for that stock (same success criteria as \code{\link{fetchSagTs}}).
#' Default \code{yearBack = 0} therefore tries the current calendar year only;
#' set \code{yearBack = 1} to also try the previous year, or \code{2} for a
#' further step back. Different stocks in the returned \code{ts}/\code{rfpts} may
#' therefore come from different assessment years; the year is recorded in
#' \code{ts$assYear} (and \code{rfpts$assYear}). A short
#' per-sid year summary is messaged when not \code{quiet}. An explicitly
#' supplied \code{year} is used for all sids with no per-sid fallback
#' (\code{yearBack} is ignored).
#'
#' Time series and reference points are fetched from the ICES SAG web API
#' (\code{icesSAG}). No project paths are resolved. For an optional offline
#' ref-point CSV, call \code{\link{loadSagRefpts}(..., root = "path/to/sdGraphs")}
#' directly with an explicit folder path.
#'
#' \strong{Missing stocks.} Codes that are not on ICES SAG (e.g. ICCAT
#' \code{alb-n}), or that fail for every tried assessment year, are skipped
#' with a warning; remaining stocks are returned. Derive returned ids with
#' \code{unique(ts$sid)}; skipped ids with
#' \code{setdiff(requested, unique(ts$sid))}. BIM fishery metadata is
#' separate: \code{bimSids()}. If every requested sid fails, stops.
#'
#' \strong{Nephrops.} For each id that starts with \code{nep.}, total biomass
#' is requested, incomplete assessment-year rows are dropped, and missing
#' catches are allowed. Early horse-mackerel years (before 1985) are dropped
#' when that stock is included.
#'
#' @param sid Character vector of ICES FishStock codes. Default
#'   \code{NULL} uses \code{bimSids()$sid}.
#' @param year Assessment / SAG year. Default \code{NULL} cascades
#'   independently per sid from the current calendar year back
#'   \code{yearBack} years. An explicit year is used for all sids with no
#'   cascade.
#' @param yearBack Non-negative integer. When \code{year = NULL}, maximum
#'   steps back from the current calendar year to try per sid (inclusive).
#'   Default \code{0} (current calendar year only). Ignored when
#'   \code{year} is set.
#' @param quiet If \code{TRUE}, suppress SAG download messages/warnings
#'   (for knitted reports).
#' @param root Optional path to a local \code{sdGraphs} directory. Used for
#'   reference points, and as a fallback for time series when the SAG API
#'   returns nothing (e.g. HTTP 503).
#' @return Named list with only:
#' \describe{
#'   \item{\code{ts}}{Standardised SAG time series (requested ids with data).}
#'   \item{\code{rfpts}}{Reference points (\code{\link{loadSagRefpts}}) for
#'     the same set, using each sid's chosen assessment year.}
#' }
#' @seealso \code{\link{fetchSagTs}}, \code{\link{fetchSagRefpts}},
#'   \code{\link{readSag}}, \code{\link{loadAdviceFlqs}}, \code{\link{bimSids}}
#' @export
#' @examples
#' \dontrun{
#' # year = NULL: each sid tries current calendar year only (yearBack = 0)
#' mac = getSAG("mac.27.nea")
#' names(mac)  # ts, rfpts
#' unique(mac$ts$sid)
#' utils::tail(mac$ts, 5)
#'
#' # catch bridge is separate
#' advice = loadAdviceFlqs(sid = unique(mac$ts$sid))
#'
#' requested = bimSids()$sid
#' all = getSAG(requested)
#' # mixed assessment years are possible when year = NULL and yearBack > 0
#' unique(all$ts[, c("sid", "assYear")])
#' setdiff(requested, unique(all$ts$sid))  # typically includes alb-n
#'
#' # allow previous year / further steps back (current, -1, -2)
#' all2 = getSAG(requested, yearBack = 2)
#'
#' pel = getSAG(subset(bimSids(), fishery == "Pelagics")$sid)
#' # explicit year: no per-sid cascade (yearBack ignored)
#' mac25 = getSAG("mac.27.nea", year = 2025)
#'
#' # local CSVs if the API is down
#' sag = getSAG(sids, yearBack = 3, root = "data/sdGraphs")
#' }
getSAG <- function(sid = NULL,
                   year = NULL,
                   yearBack = 0L,
                   quiet = FALSE,
                   root = NULL) {
  if (is.null(sid))
    sid = bimSids()$sid
  sid = unique(as.character(sid))
  if (!length(sid) || any(!nzchar(sid)))
    stop("sid must be a non-empty character vector.", call. = FALSE)

  currentYear = as.integer(format(Sys.Date(), "%Y"))
  if (is.null(year)) {
    yearBack = as.integer(yearBack)
    if (length(yearBack) != 1L || is.na(yearBack) || yearBack < 0L)
      stop("yearBack must be a single non-negative integer.", call. = FALSE)
    yearsToTry = as.integer(currentYear - seq.int(0L, yearBack))
  } else {
    yearsToTry = as.integer(year)
    if (length(yearsToTry) != 1L || is.na(yearsToTry))
      stop("year must be a single integer or NULL.", call. = FALSE)
  }

  run <- function() {
    tsParts = vector("list", length(sid))
    yearUsed = setNames(rep(NA_integer_, length(sid)), sid)

    for (i in seq_along(sid)) {
      stock = sid[[i]]
      isNep = grepl("^nep\\.", stock)
      for (y in yearsToTry) {
        out = fetchSagTsOne(
          stock = stock,
          sagYear = y,
          includeTB = isNep,
          dropAssYearRows = isNep,
          requireCatches = !isNep)
        if (is.null(out) || !nrow(out)) next
        # Horse mackerel early years are not used in BIM workflows
        out = subset(out, !(sid == "hom.27.2a3a4a5b6a7a-ce-k8" & year < 1985))
        if (!nrow(out)) next
        tsParts[[i]] = out
        yearUsed[[i]] = y
        break
      }
    }

    missed = names(yearUsed)[is.na(yearUsed)]
    if (length(missed) && !is.null(root) && nzchar(root)) {
      if (!quiet)
        message("SAG API missed ", length(missed),
                " sid(s); reading local CSVs under ", root)
      local = tryCatch(
        readSag(root, sid = missed, assYear = yearsToTry, verbose = !quiet),
        error = function(e) {
          warning("Local SAG read failed: ", conditionMessage(e), call. = FALSE)
          NULL
        })
      if (!is.null(local) && is.data.frame(local$ts) && nrow(local$ts)) {
        idCol = if ("sid" %in% names(local$ts)) "sid" else ".id"
        yrCol = if ("assYear" %in% names(local$ts)) "assYear" else "assYr"
        for (stock in missed) {
          d = local$ts[local$ts[[idCol]] == stock, , drop = FALSE]
          if (!nrow(d)) next
          y = max(as.integer(d[[yrCol]]), na.rm = TRUE)
          tsParts[[match(stock, sid)]] = d[as.integer(d[[yrCol]]) == y, , drop = FALSE]
          yearUsed[[stock]] = y
        }
      }
    }

    ts = rbindFill(tsParts)
    if (is.null(ts) || !nrow(ts)) {
      hint = if (is.null(root) || !nzchar(root))
        " Pass root = 'path/to/sdGraphs' to use local CSVs when the API is down."
      else
        ""
      stop("No usable SAG data for any requested sid after trying year(s): ",
           paste(yearsToTry, collapse = ", "), ".", hint, call. = FALSE)
    }

    got = names(yearUsed)[!is.na(yearUsed)]
    skipped = setdiff(sid, got)
    if (length(skipped) && !quiet)
      warning("No SAG series for: ", paste(skipped, collapse = ", "),
              call. = FALSE)

    if (!quiet) {
      summ = data.frame(
        sid = got,
        assYear = unname(yearUsed[got]),
        stringsAsFactors = FALSE)
      message("SAG assessment year by sid:")
      message(paste(utils::capture.output(print(summ, row.names = FALSE)),
                    collapse = "\n"))
    }

    # Refpts for each sid at the same assessment year used for ts
    rfParts = lapply(got, function(stock) {
      y = yearUsed[[stock]]
      tryCatch(
        loadSagRefpts(sid = stock, year = y, root = root),
        error = function(e) {
          warning("SAG refpts skipped for ", stock, " (", y, "): ",
                  conditionMessage(e), call. = FALSE)
          NULL
        })
    })
    rfpts = rbindFill(rfParts)
    if (is.null(rfpts) || !nrow(rfpts))
      stop("No SAG reference points returned for stocks with time series.",
           call. = FALSE)

    list(ts = ts, rfpts = rfpts)
  }

  if (isTRUE(quiet))
    suppressMessages(suppressWarnings(run()))
  else
    run()
}

# --- SAG S4 helpers (readSag / getSag / meltFl) ------------------------------

.sagAsNum <- function(x) suppressWarnings(as.numeric(as.character(x)))

.sagRefptMap <- c(
  Flim = "flim", Fpa = "fpa", Bpa = "bpa", Blim = "blim",
  FMSY = "fmsy", BMSY = "bmsy", MSYBtrigger = "msybtrigger",
  Fmanagement = "fmanagement", Bmanagement = "bmanagement",
  TerminalYear = "terminalyear")

.sagIdKeys <- c(".id", "sid", "assYr", "assYear", "year", "age", "unit", "iter")

.sagStdKeys <- function(df) {
  if (is.null(df) || !is.data.frame(df) || !ncol(df)) return(df)
  if ("FishStock" %in% names(df) && !"sid" %in% names(df))
    df$sid <- df$FishStock
  if (".id" %in% names(df) && !"sid" %in% names(df))
    df$sid <- df$.id
  if ("sid" %in% names(df) && !".id" %in% names(df))
    df$.id <- df$sid
  if ("AssessmentYear" %in% names(df) && !"assYear" %in% names(df))
    df$assYear <- df$AssessmentYear
  if ("assYr" %in% names(df) && !"assYear" %in% names(df))
    df$assYear <- df$assYr
  if ("assYear" %in% names(df) && !"assYr" %in% names(df))
    df$assYr <- df$assYear
  df
}

.rbindFill <- function(dfs) {
  dfs = Filter(function(d) is.data.frame(d) && nrow(d), dfs)
  if (!length(dfs)) return(NULL)
  cols = unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs = lapply(dfs, function(d) {
    for (m in setdiff(cols, names(d))) d[[m]] = NA
    d[cols]
  })
  out = do.call(rbind, dfs)
  rownames(out) = NULL
  out
}

#' @rdname readSagRefpts
#' @export
readSagRefpts <- function(file, sid = NULL) {
  if (is.character(file) && length(file) == 1L)
    file = utils::read.csv(file, stringsAsFactors = FALSE)
  if (!is.data.frame(file))
    stop("'file' must be a path or data.frame.", call. = FALSE)
  df = .sagStdKeys(file)
  idCol = if ("sid" %in% names(df)) "sid" else if (".id" %in% names(df)) ".id"
  if (is.null(idCol))
    stop("SAG CSV missing stock id column.", call. = FALSE)
  assCol = if ("assYear" %in% names(df)) "assYear" else "assYr"
  if (!assCol %in% names(df))
    stop("SAG CSV missing assessment year column.", call. = FALSE)
  if (!is.null(sid))
    df = df[df[[idCol]] %in% sid, , drop = FALSE]
  key = df[[idCol]]
  ass = df[[assCol]]
  keep = !duplicated(paste(key, ass))
  out = df[keep, , drop = FALSE]
  for (src in names(.sagRefptMap)) {
    dst = .sagRefptMap[[src]]
    if (src %in% names(out))
      out[[dst]] = .sagAsNum(out[[src]])
  }
  keep = unique(c(".id", "sid", "assYr", "assYear", intersect(unname(.sagRefptMap), names(out))))
  out = .sagStdKeys(out[, keep, drop = FALSE])
  rownames(out) = NULL
  out
}

#' @rdname readSagTs
#' @export
readSagTs <- function(file, sid = NULL) {
  if (is.character(file) && length(file) == 1L)
    file = utils::read.csv(file, stringsAsFactors = FALSE)
  if (!is.data.frame(file))
    stop("'file' must be a path or data.frame.", call. = FALSE)
  df = .sagStdKeys(file)
  idCol = if ("sid" %in% names(df)) "sid" else ".id"
  assCol = if ("assYear" %in% names(df)) "assYear" else "assYr"
  yearCol = if ("year" %in% names(df)) "year" else if ("Year" %in% names(df)) "Year"
  if (is.null(yearCol))
    stop("SAG CSV missing year column.", call. = FALSE)
  if (!is.null(sid))
    df = df[df[[idCol]] %in% sid, , drop = FALSE]
  pick <- function(cands) {
    hit = cands[cands %in% names(df)]
    if (length(hit)) hit[1L] else NA_character_
  }
  ssbCol = pick(c("ssb", "StockSize", "TBiomass"))
  fCol = pick(c("f", "fbar", "FishingPressure", "F"))
  recCol = pick(c("rec", "Recruitment", "recruitment"))
  catchCol = pick(c("catch", "Catches", "catches"))
  out = data.frame(
    sid = as.character(df[[idCol]]),
    assYear = as.integer(df[[assCol]]),
    year = .sagAsNum(df[[yearCol]]),
    ssb = if (!is.na(ssbCol)) .sagAsNum(df[[ssbCol]]) else NA_real_,
    f = if (!is.na(fCol)) .sagAsNum(df[[fCol]]) else NA_real_,
    rec = if (!is.na(recCol)) .sagAsNum(df[[recCol]]) else NA_real_,
    catch = if (!is.na(catchCol)) .sagAsNum(df[[catchCol]]) else NA_real_,
    stringsAsFactors = FALSE)
  out = .sagStdKeys(out)
  out = out[is.finite(out$year), , drop = FALSE]
  rownames(out) = NULL
  out
}

#' @rdname readSagDir
#' @export
readSagDir <- function(root, sid = NULL, assYear = NULL, verbose = TRUE) {
  if (is.null(assYear)) {
    yrs = list.dirs(root, full.names = FALSE, recursive = FALSE)
    assYear = sort(as.integer(yrs[grepl("^[0-9]{4}$", yrs)]))
  }
  assYear = as.integer(assYear)
  tsParts = rfParts = vector("list", length(assYear))
  for (i in seq_along(assYear)) {
    y = assYear[[i]]
    csv = findSagGraphsCsv(year = y, root = root)
    if (is.null(csv)) {
      if (verbose)
        message("No StockAssessmentGraphs CSV for ", y, " under ", root)
      next
    }
    if (verbose) message("Reading ", csv)
    raw = utils::read.csv(csv, stringsAsFactors = FALSE)
    if (!"AssessmentYear" %in% names(raw) || all(is.na(raw$AssessmentYear)))
      raw$AssessmentYear = y
    tsParts[[i]] = readSagTs(raw, sid = sid)
    rfParts[[i]] = readSagRefpts(raw, sid = sid)
  }
  list(ts = .rbindFill(tsParts), rfpts = .rbindFill(rfParts))
}

#' @export
setGeneric("readSag", function(object, ...) standardGeneric("readSag"))

#' @rdname readSag
#' @export
setMethod("readSag", signature(object = "character"),
          function(object, sid = NULL, assYear = NULL, verbose = TRUE, ...) {
            if (length(object) != 1L)
              stop("'object' must be a single file or directory path.", call. = FALSE)
            if (dir.exists(object))
              return(readSagDir(object, sid = sid, assYear = assYear, verbose = verbose))
            list(
              ts = readSagTs(object, sid = sid),
              rfpts = readSagRefpts(object, sid = sid))
          })

#' @rdname readSag
#' @export
setMethod("readSag", signature(object = "data.frame"),
          function(object, sid = NULL, ...) {
            list(
              ts = readSagTs(object, sid = sid),
              rfpts = readSagRefpts(object, sid = sid))
          })

.sagPull <- function(stock, assYear, data = "summary", purpose = "Advice") {
  if (!requireNamespace("icesSAG", quietly = TRUE))
    stop("Package 'icesSAG' is required.", call. = FALSE)
  tryCatch(
    icesSAG::getSAG(stock = stock, year = assYear, data = data, purpose = purpose),
    error = function(e) {
      warning(stock, " ", assYear, " (", data, "): ", conditionMessage(e),
              call. = FALSE)
      NULL
    })
}

.sagTidySummary <- function(dat, sid, assYear) {
  if (is.null(dat) || !is.data.frame(dat) || !nrow(dat)) return(NULL)
  pick <- function(cands) {
    hit = cands[tolower(cands) %in% tolower(names(dat))]
    if (length(hit)) names(dat)[match(tolower(hit[1L]), tolower(names(dat)))] else NA_character_
  }
  yearCol = pick(c("Year", "year"))
  if (is.na(yearCol)) return(NULL)
  ssbCol = pick(c("SSB", "ssb", "StockSize", "TBiomass"))
  fCol = pick(c("F", "f", "fbar", "FishingPressure"))
  recCol = pick(c("recruitment", "Recruitment", "rec", "recruits"))
  catchCol = pick(c("catches", "Catches", "catch"))
  landCol = pick(c("landings", "Landings"))
  discCol = pick(c("discards", "Discards"))
  catches = if (!is.na(catchCol)) .sagAsNum(dat[[catchCol]]) else NA_real_
  landings = if (!is.na(landCol)) .sagAsNum(dat[[landCol]]) else NA_real_
  discards = if (!is.na(discCol)) .sagAsNum(dat[[discCol]]) else NA_real_
  catch = ifelse(is.finite(catches), catches,
                 landings + ifelse(is.finite(discards), discards, 0))
  out = data.frame(
    sid = as.character(sid),
    assYear = as.integer(assYear),
    year = .sagAsNum(dat[[yearCol]]),
    ssb = if (!is.na(ssbCol)) .sagAsNum(dat[[ssbCol]]) else NA_real_,
    f = if (!is.na(fCol)) .sagAsNum(dat[[fCol]]) else NA_real_,
    rec = if (!is.na(recCol)) .sagAsNum(dat[[recCol]]) else NA_real_,
    catch = catch,
    stringsAsFactors = FALSE)
  out = .sagStdKeys(out)
  out = out[is.finite(out$year), , drop = FALSE]
  out[order(out$year), , drop = FALSE]
}

.sagTidyRefpts <- function(dat, sid, assYear) {
  if (is.null(dat) || !is.data.frame(dat) || !nrow(dat)) return(NULL)
  r = dat[nrow(dat), , drop = FALSE]
  out = data.frame(sid = as.character(sid), assYear = as.integer(assYear),
                   stringsAsFactors = FALSE)
  for (src in names(.sagRefptMap)) {
    dst = .sagRefptMap[[src]]
    out[[dst]] = if (src %in% names(r)) .sagAsNum(r[[src]]) else NA_real_
  }
  .sagStdKeys(out)
}

.sagGrid <- function(sid, assYear) {
  sid = unique(as.character(sid))
  assYear = unique(as.integer(assYear))
  expand.grid(sid = sid, assYear = assYear, stringsAsFactors = FALSE)
}

#' @rdname getSagTs
#' @export
getSagTs <- function(sid, assYear, purpose = "Advice", dropAssYearRows = FALSE) {
  grid = .sagGrid(sid, assYear)
  parts = vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    parts[[i]] = .sagTidySummary(
      .sagPull(grid$sid[[i]], grid$assYear[[i]], "summary", purpose),
      grid$sid[[i]], grid$assYear[[i]])
  }
  out = .rbindFill(parts)
  if (is.null(out)) {
    out = data.frame(sid = character(), assYear = integer(), year = integer(),
                    ssb = double(), f = double(), rec = double(), catch = double(),
                    stringsAsFactors = FALSE)
  }
  if (dropAssYearRows && nrow(out))
    out = subset(out, year != assYear)
  .sagStdKeys(out)
}

#' @rdname getSagRefpts
#' @export
getSagRefpts <- function(sid, assYear, purpose = "Advice") {
  grid = .sagGrid(sid, assYear)
  parts = vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    parts[[i]] = .sagTidyRefpts(
      .sagPull(grid$sid[[i]], grid$assYear[[i]], "refpts", purpose),
      grid$sid[[i]], grid$assYear[[i]])
  }
  out = .rbindFill(parts)
  if (is.null(out)) {
    out = data.frame(sid = character(), assYear = integer(),
                     stringsAsFactors = FALSE)
    for (nm in unname(.sagRefptMap)) out[[nm]] = double()
  }
  .sagStdKeys(out)
}

#' @export
setGeneric("getSag", function(sid, assYear, ...) standardGeneric("getSag"))

.getSagList <- function(sid, assYear, purpose = "Advice", dropAssYearRows = FALSE) {
  list(
    ts = getSagTs(sid, assYear, purpose = purpose, dropAssYearRows = dropAssYearRows),
    rfpts = getSagRefpts(sid, assYear, purpose = purpose))
}

#' @rdname getSag
#' @export
setMethod("getSag", signature(sid = "character", assYear = "numeric"),
          function(sid, assYear, purpose = "Advice", dropAssYearRows = FALSE, ...) {
            .getSagList(sid, assYear, purpose, dropAssYearRows)
          })

#' @rdname getSag
#' @export
setMethod("getSag", signature(sid = "character", assYear = "integer"),
          function(sid, assYear, purpose = "Advice", dropAssYearRows = FALSE, ...) {
            .getSagList(sid, assYear, purpose, dropAssYearRows)
          })

#' @rdname getSag
#' @export
setMethod("getSag", signature(sid = "data.frame", assYear = "missing"),
          function(sid, assYear, purpose = "Advice", dropAssYearRows = FALSE, ...) {
            sid = .sagStdKeys(sid)
            idCol = if ("sid" %in% names(sid)) "sid" else ".id"
            assCol = if ("assYear" %in% names(sid)) "assYear" else "assYr"
            .getSagList(sid[[idCol]], sid[[assCol]], purpose, dropAssYearRows)
          })

#' @export
setGeneric("meltFl", function(object, ...) standardGeneric("meltFl"))

.meltFlCore <- function(object, id.vars = NULL, measure.vars = NULL, na.rm = FALSE) {
  if (is.null(id.vars))
    id.vars = intersect(.sagIdKeys, names(object))
  if (is.null(measure.vars))
    measure.vars = setdiff(names(object), id.vars)
  if (!length(measure.vars))
    stop("No measure columns to melt.", call. = FALSE)
  n = nrow(object)
  k = length(measure.vars)
  out = object[rep(seq_len(n), times = k), id.vars, drop = FALSE]
  out$qname = factor(rep(measure.vars, each = n), levels = measure.vars)
  out$data = unlist(object[, measure.vars, drop = FALSE], use.names = FALSE)
  if (na.rm) out = out[is.finite(out$data), , drop = FALSE]
  rownames(out) = NULL
  out
}

#' @rdname meltFl
#' @export
setMethod("meltFl", signature(object = "data.frame"),
          function(object, id.vars = NULL, measure.vars = NULL, na.rm = FALSE, ...) {
            .meltFlCore(object, id.vars, measure.vars, na.rm)
          })

#' @rdname meltFl
#' @export
setMethod("meltFl", signature(object = "matrix"),
          function(object, id.vars = NULL, measure.vars = NULL, na.rm = FALSE, ...) {
            .meltFlCore(as.data.frame(object, stringsAsFactors = FALSE),
                        id.vars, measure.vars, na.rm)
          })

#' @rdname meltFl
#' @export
setMethod("meltFl", signature(object = "FLQuant"),
          function(object, ...) {
            .meltFlCore(as.data.frame(object))
          })

#' @rdname meltFl
#' @export
setMethod("meltFl", signature(object = "FLQuants"),
          function(object, ...) {
            .meltFlCore(as.data.frame(object))
          })
