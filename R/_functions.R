bf_01 <- function(x, prior_sigma) {
  dens <- density(x = as.numeric(x))
  dens_at_0 <- approx(dens$x, dens$y, xout = 0, n = 1000, rule = 2)$y
  dnorm(x = 0, mean = 0, sd = prior_sigma) / dens_at_0
}

summary_posterior <- function(x, prior_sigma = 3) {
  estimate = round(x = mean(x), digits = 2)
  ci = round(x = tidybayes::hdi(x), digits = 2)
  ci = paste0("[", paste0(ci, collapse = ", "), "]")
  pd_side = if (estimate < 0) {x < 0} else {x > 0}
  pd = round(x = sum(pd_side)/length(x), digits = 3)
  ps_side = if (estimate < 0) {x < -0.1} else {x > 0.1}
  ps = round(x = sum(ps_side)/length(x), digits = 3)
  bf = round(x = bf_01(x, prior_sigma), digits = 2)
  ess = round(x = posterior::ess_basic(x), digits = 1)
  rhat = round(x = posterior::rhat(x), digits = 3)
  
  list(estimate = estimate, ci = ci,
       pd = pd, ps = ps, bf = bf,
       ess = ess, rhat = rhat)
}

report_posterior <- function(x, prior_sigma = 3) {
  m_report <- summary_posterior(x, prior_sigma) |>
    as.data.table()
  m_report[, list(
    effect = paste0("$\beta$ = ", estimate, ", CI~95%~", ci),
    significance = paste0("pd = ", pd*100, "%, ps = ", ps*100, "%, BF~10~ = ", bf),
    convergence = paste0("ESS = ", ess, ", R-hat = ", rhat)
  )]
}

summary_model <- function(model, variable = NULL) {
  brms::as_draws_df(model, variable = variable, regex = TRUE) |>
    lapply(summary_posterior) |>
    rbindlist(idcol = "var")
}

report_model <- function(model, variable = NULL) {
  m_report <- summary_model(model, variable)
  m_report[!var %in% c(".chain", ".iteration", ".draw"), list(
    effect = paste0("$\beta$ = ", estimate, ", CI~95%~", ci),
    significance = paste0("pd = ", pd*100, "%, ps = ", ps*100, "%, BF~10~ = ", bf),
    convergence = paste0("ESS = ", ess, ", R-hat = ", rhat)
  ), by = "var"][, var := gsub("_", " ", var)][]
}
