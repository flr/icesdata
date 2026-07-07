#' Coerce JABBA fit to a standard list
#' Minimal assumption: we can build chains, residuals, ts from the raw fit.
as_jabba_tidy <- function(raw_fit, stock="", scenario = "base") {
  # TODO: adapt this to your actual JABBA object
  fit <- list(
    chains    = raw_fit$mc.sims,     # mcmc samples as data.frame/matrix
    residuals = raw_fit$cpue.res,    # long data frame
    ts        = raw_fit$timeseries   # yearly summaries
  )
  attr(fit, "stock")    <- stock
  attr(fit, "scenario") <- scenario
  class(fit) <- c("jabba_tidy", class(raw_fit))
  fit
}

library(dplyr)
library(tidyr)
library(ggplot2)

#' Tidy MCMC chains
tidy_chains <- function(fit, params) {
  ch <- as.data.frame(fit$chains)
  ch <- ch[, intersect(colnames(ch), params), drop = FALSE]
  ch$iter <- seq_len(nrow(ch))
  ch |>
    pivot_longer(-iter, names_to = "param", values_to = "value")
}

#' Trace + density plots
gg_mcmc_trace_density <- function(fit, params = c("r", "K", "sigma.proc", "sigma.obs")) {
  df <- tidy_chains(fit, params)
  stock <- attr(fit, "stock")
  
  p_trace <- ggplot(df, aes(iter, value)) +
    geom_line(alpha = 0.6) +
    facet_wrap(~ param, scales = "free_y") +
    labs(x = "Iteration", y = "Value",
         title = paste0(stock, " – MCMC traces"))
  
  p_dens <- ggplot(df, aes(value)) +
    geom_histogram(aes(y = ..density..), bins = 40,
                   fill = "grey80", colour = "grey40") +
    facet_wrap(~ param, scales = "free") +
    labs(x = "Parameter value", y = "Density",
         title = paste0(stock, " – posterior densities"))
  
  list(trace = p_trace, density = p_dens)
}

# expects residuals: year, fleet, resid, obs, pred
gg_residuals_time <- function(fit) {
  res <- fit$residuals
  stock <- attr(fit, "stock")
  
  ggplot(res, aes(year, resid)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_line(colour = "grey40") +
    geom_point(aes(colour = fleet)) +
    facet_wrap(~ fleet, scales = "free_x") +
    labs(x = "Year", y = "Standardised residual",
         title = paste0(stock, " – CPUE residuals over time"))
}

gg_obs_pred <- function(fit) {
  res <- fit$residuals
  stock <- attr(fit, "stock")
  
  ggplot(res, aes(pred, obs, colour = fleet)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    geom_point(alpha = 0.7) +
    facet_wrap(~ fleet, scales = "free") +
    labs(x = "Predicted index", y = "Observed index",
         title = paste0(stock, " – observed vs predicted"))
}

gg_residuals_qq <- function(fit) {
  res <- fit$residuals
  stock <- attr(fit, "stock")
  
  ggplot(res, aes(sample = resid)) +
    stat_qq() + stat_qq_line() +
    facet_wrap(~ fleet) +
    labs(x = "Theoretical quantiles", y = "Sample quantiles",
         title = paste0(stock, " – residual QQ‑plots"))
}

# expects ts: year, Bmed, Blow, Bhigh, Umed, Ulow, Uhigh, B_Bmsy, U_Umsy
gg_traj_B_U <- function(fit) {
  ts <- fit$ts
  stock <- attr(fit, "stock")
  
  pB <- ggplot(ts, aes(year, Bmed)) +
    geom_ribbon(aes(ymin = Blow, ymax = Bhigh), alpha = 0.25) +
    geom_line() +
    labs(x = "Year", y = "Biomass",
         title = paste0(stock, " – biomass with uncertainty"))
  
  pU <- ggplot(ts, aes(year, Umed)) +
    geom_ribbon(aes(ymin = Ulow, ymax = Uhigh), alpha = 0.25) +
    geom_line() +
    labs(x = "Year", y = "Harvest rate",
         title = paste0(stock, " – harvest rate with uncertainty"))
  
  list(B = pB, U = pU)
}

gg_traj_rel <- function(fit) {
  ts <- fit$ts
  stock <- attr(fit, "stock")
  
  pBB <- ggplot(ts, aes(year, B_Bmsy)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
    geom_line() +
    labs(x = "Year", y = "B / Bmsy",
         title = paste0(stock, " – relative biomass"))
  
  pUU <- ggplot(ts, aes(year, U_Umsy)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
    geom_line() +
    labs(x = "Year", y = "U / Umsy",
         title = paste0(stock, " – relative harvest"))
  
  list(BB = pBB, UU = pUU)
}

#' Prior–posterior overlay using a supplied prior density
#' priors: named list of density functions, e.g. list(r = function(x) dlnorm(x, m, s), ...)
gg_prior_posterior <- function(fit, params, priors, n_grid = 200) {
  ch <- as.data.frame(fit$chains)
  stock <- attr(fit, "stock")
  
  lapply(params, function(p) {
    stopifnot(p %in% names(priors))
    
    post <- ch[[p]]
    rng  <- range(post, na.rm = TRUE)
    xseq <- seq(rng[1], rng[2], length.out = n_grid)
    
    df_prior <- data.frame(
      x   = xseq,
      pdf = priors[[p]](xseq),
      type = "prior"
    )
    df_post <- data.frame(
      x   = density(post, na.rm = TRUE)$x,
      pdf = density(post, na.rm = TRUE)$y,
      type = "posterior"
    )
    df <- rbind(df_prior, df_post)
    
    ggplot(df, aes(x, pdf, colour = type)) +
      geom_line() +
      labs(x = p, y = "Density",
           colour = "",
           title = paste0(stock, " – prior vs posterior: ", p))
  })
}

#' JABBA B vs ICES SSB
#' ices_ts: year, SSB (or TB)
gg_external_ices <- function(fit, ices_ts) {
  ts <- fit$ts
  stock <- attr(fit, "stock")
  
  df <- full_join(
    ts |> select(year, JABBA_B = Bmed),
    ices_ts |> select(year, ICES_SSB = SSB),
    by = "year"
  )
  
  ggplot(df, aes(year)) +
    geom_line(aes(y = JABBA_B, colour = "JABBA B")) +
    geom_line(aes(y = ICES_SSB, colour = "ICES SSB")) +
    labs(x = "Year", y = "Biomass / SSB", colour = "",
         title = paste0(stock, " – JABBA vs ICES biomass"))
}

#' Sum of substocks vs ICES
gg_external_sum_substocks <- function(fits, ices_ts) {
  ts_all <- bind_rows(lapply(fits, function(f) {
    ts <- f$ts
    ts$stock <- attr(f, "stock")
    ts
  }))
  
  sum_ts <- ts_all |>
    group_by(year) |>
    summarise(JABBA_sumB = sum(Bmed, na.rm = TRUE), .groups = "drop")
  
  df <- full_join(
    sum_ts,
    ices_ts |> select(year, ICES_SSB = SSB),
    by = "year"
  )
  
  ggplot(df, aes(year)) +
    geom_line(aes(y = JABBA_sumB, colour = "Sum JABBA B")) +
    geom_line(aes(y = ICES_SSB, colour = "ICES SSB")) +
    labs(x = "Year", y = "Biomass / SSB", colour = "",
         title = "Sum of substocks vs ICES SSB")
}

#' External survey vs B
#' external_idx: year, index
gg_external_survey <- function(fit, external_idx) {
  ts <- fit$ts
  stock <- attr(fit, "stock")
  
  df <- full_join(
    ts |> select(year, B = Bmed),
    external_idx |> select(year, index),
    by = "year"
  ) |>
    mutate(
      B_std   = as.numeric(scale(B)),
      idx_std = as.numeric(scale(index))
    )
  
  ggplot(df, aes(year)) +
    geom_line(aes(y = B_std, colour = "JABBA B")) +
    geom_line(aes(y = idx_std, colour = "External index")) +
    labs(x = "Year", y = "Standardised units", colour = "",
         title = paste0(stock, " – biomass vs external index"))
}

summarise_status <- function(fit, ref_year = max(fit$ts$year)) {
  ts   <- fit$ts
  row  <- ts[ts$year == ref_year, ]
  
  data.frame(
    stock     = attr(fit, "stock"),
    scenario  = attr(fit, "scenario"),
    year      = ref_year,
    B_Bmsy    = row$B_Bmsy,
    U_Umsy    = row$U_Umsy,
    # if you have per‑year posterior arrays B_Bmsy_post, etc.:
    # P_BltBmsy = mean(fit$post$B_Bmsy[fit$post$year == ref_year] < 1),
    # P_UgtUmsy = mean(fit$post$U_Umsy[fit$post$year == ref_year] > 1),
    stringsAsFactors = FALSE
  )
}

run_core_diags <- function(fits, ices_ts, external_idx_list = NULL,
                           params = c("r", "K", "sigma.proc", "sigma.obs")) {
  lapply(fits, function(fit) {
    list(
      chains   = gg_mcmc_trace_density(fit, params),
      resid_ts = gg_residuals_time(fit),
      obs_pred = gg_obs_pred(fit),
      qq       = gg_residuals_qq(fit),
      traj     = gg_traj_B_U(fit),
      rel      = gg_traj_rel(fit),
      ext_ices = gg_external_ices(fit, ices_ts),
      ext_survey = if (!is.null(external_idx_list)) {
        gg_external_survey(fit, external_idx_list[[attr(fit, "stock")]])
      } else NULL
    )
  })
}