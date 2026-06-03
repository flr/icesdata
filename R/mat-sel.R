#' Interpolate age at a target value along an age–vector curve
#'
#' Wrapper for \code{\link{ageAtTarget}} (same linear bracketing rule as
#' \code{\link{mat50}} / \code{\link{sel50}}).
#'
#' @param data Numeric vector (e.g. maturity or selectivity by age).
#' @param age Numeric vector of ages (same length as \code{data}).
#' @param target Value to interpolate (default 0.5).
#' @return Numeric age, or \code{NA} if not bracketed.
#' @export
findAge <- function(data, age = seq_along(data) - 1L, target = 0.5) {
  ageAtTarget(value = data, age = age, target = target)
}

#' Normalized overlap between two age-specific curves
#'
#' @param curve1,curve2 Numeric vectors of the same length.
#' @return Proportion in \eqn{[0, 1]}.
#' @export
calcOverlap <- function(curve1, curve2) {
  norm1 <- curve1 / max(curve1, na.rm = TRUE)
  norm2 <- curve2 / max(curve2, na.rm = TRUE)
  sum(pmin(norm1, norm2), na.rm = TRUE) / sum(pmax(norm1, norm2), na.rm = TRUE)
}

#' Maturity vs selectivity metrics by age for one \code{FLStock}
#'
#' Per-year \code{mat50} and \code{sel50} use \code{\link{mat50}} and
#' \code{\link{sel50}}; cumulative schedules and lag metrics are by age row.
#'
#' @param x An \code{FLStock} with \code{catch.sel} and \code{mat} slots.
#' @return A data frame (from \code{model.frame} on \code{FLQuants}) with extra
#'   columns \code{cumSel}, \code{cumMat}, \code{dif}, \code{mat50}, \code{sel50},
#'   \code{maxDifAge}, \code{maxDifVal}.
#' @export
matSel <- function(x) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for matSel().", call. = FALSE)
  }
  metrics <- model.frame(
    FLCore::FLQuants(
      x,
      sel = function(x) FLCore::catch.sel(x) / max(FLCore::catch.sel(x)),
      mat = FLCore::mat
    ),
    drop = TRUE
  )
  m50 <- mat50(x)
  s50 <- sel50(x)
  metrics$mat50 <- unname(m50[as.character(metrics$year)])
  metrics$sel50 <- unname(s50[as.character(metrics$year)])
  suppressWarnings(
    metrics <- dplyr::mutate(
      metrics,
      cumSel = cumsum(sel) / max(cumsum(sel)),
      cumMat = cumsum(mat) / max(cumsum(mat)),
      dif = cumSel - cumMat,
      maxDifAge = age[which.max(dif)],
      maxDifVal = max(dif)
    )
  )
  metrics
}

#' Summarise maturity–selectivity metrics for a list of stocks
#'
#' @param x A named list of \code{FLStock} objects (e.g. from \code{icesdata}).
#' @return A data frame with one row per stock (\code{.id}) and columns
#'   \code{ageMaxDif}, \code{mat50}, \code{sel50}, \code{selMatLag},
#'   \code{vuln}, \code{overlap}.
#' @export
selMetrics <- function(x) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for selMetrics().", call. = FALSE)
  }
  if (!requireNamespace("plyr", quietly = TRUE)) {
    stop("Package 'plyr' is required for selMetrics().", call. = FALSE)
  }
  dat <- plyr::ldply(x, matSel)
  dplyr::summarise(
    dplyr::group_by(dat, .data$.id),
    ageMaxDif = age[which.max(dif)][1],
    mat50 = mean(mat50, na.rm = TRUE),
    sel50 = mean(sel50, na.rm = TRUE),
    selMatLag = mean(sel50 - mat50, na.rm = TRUE),
    vuln = sum(sel * mat, na.rm = TRUE) / sum(mat, na.rm = TRUE),
    overlap = calcOverlap(sel, mat),
    .groups = "drop"
  )
}

#' Scatterplot of age at 50\% selectivity vs maturity
#'
#' @param x Data frame with columns \code{mat50} and \code{sel50} (e.g. from
#'   \code{\link{selMetrics}}).
#' @return A \code{ggplot} object.
#' @export
plotMatSel <- function(x) {
  ggplot2::ggplot(x, ggplot2::aes(x = .data$mat50, y = .data$sel50)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::labs(
      x = "Age at 50% Maturity",
      y = "Age at 50% Selectivity",
      title = "Comparison of Maturity and Selectivity Patterns"
    ) +
    ggplot2::theme_minimal()
}
