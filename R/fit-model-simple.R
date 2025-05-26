
# Prepare workspace -------------------------------------------------------

## Load libraries
library(data.table)
library(brms)
library(ggplot2)

## Load data
data("senescence_d")

# Prepare data objects ----------------------------------------------------

# Response variables
var_response <- c("alpha","beta","c","lambda","phi","tau","delta")

# Predictors
var_predictors <- c("abc_linfocitos_total", "abc_linfocitosb_total",
                    "cd21_m_cd11c_p","cd21_p_cd11c_p",
                    "cd21_p_cd11c_m","cd21_m_cd11c_m")

var_all <- c(var_response, var_predictors)

## Get model data, only one observation per subject
m_data <- senescence_d[, .SD[1], .SDcols = var_all, id]

## Standardize response variables
m_data_std <- datawizard::standardize(m_data)

## Define the multivariate model with residual correlation
m_formula <- 
  bf(alpha ~ abc_linfocitos_total + abc_linfocitosb_total + cd21_m_cd11c_p + cd21_m_cd11c_m + cd21_p_cd11c_m + cd21_p_cd11c_p) +
  bf(beta ~ abc_linfocitos_total + abc_linfocitosb_total + cd21_m_cd11c_p + cd21_m_cd11c_m + cd21_p_cd11c_m + cd21_p_cd11c_p) +
  bf(c ~ abc_linfocitos_total + abc_linfocitosb_total + cd21_m_cd11c_p + cd21_m_cd11c_m + cd21_p_cd11c_m + cd21_p_cd11c_p) +
  bf(lambda ~ abc_linfocitos_total + abc_linfocitosb_total + cd21_m_cd11c_p + cd21_m_cd11c_m + cd21_p_cd11c_m + cd21_p_cd11c_p) +
  bf(phi ~ abc_linfocitos_total + abc_linfocitosb_total + cd21_m_cd11c_p + cd21_m_cd11c_m + cd21_p_cd11c_m + cd21_p_cd11c_p) +
  bf(tau ~ abc_linfocitos_total + abc_linfocitosb_total + cd21_m_cd11c_p + cd21_m_cd11c_m + cd21_p_cd11c_m + cd21_p_cd11c_p) +
  bf(delta ~ abc_linfocitos_total + abc_linfocitosb_total + cd21_m_cd11c_p + cd21_m_cd11c_m + cd21_p_cd11c_m + cd21_p_cd11c_p) +
  set_rescor(TRUE)

## Check what are the default priors
default_prior(m_formula, m_data_std)

## Stablish priors
priors <- c(
  set_prior("normal(0, 3)", class = "b", resp = var_response),
  set_prior("normal(0, 3)", class = "Intercept", resp = var_response),
  set_prior("normal(1, 3)", class = "sigma", resp = var_response, lb = 0),
  set_prior("lkj(1)", class = "rescor")
)

# Fit models --------------------------------------------------------------

m_mod_1 <- brm(
  formula = m_formula,
  data = m_data_std,
  family = gaussian(),
  prior = priors,
  chains = 5, cores = 5,
  iter = 4000, warmup = 2000,
  seed = 1234,
  file = "models/m_mod_1.RDS"
)

# Visualize models --------------------------------------------------------

## Inspect model parameters
bayestestR::describe_posterior(m_mod_1)

# ## Plot posterior distributions
# plot(m_mod_1)
# 
# ## Posterior predictive checks
# pp_check(m_mod_1, resp = "alpha", ndraws = 100)
# pp_check(m_mod_1, resp = "beta", ndraws = 100)
# pp_check(m_mod_1, resp = "c", ndraws = 100)
# pp_check(m_mod_1, resp = "lambda", ndraws = 100)
# pp_check(m_mod_1, resp = "phi", ndraws = 100)
# pp_check(m_mod_1, resp = "tau", ndraws = 100)
# pp_check(m_mod_1, resp = "delta", ndraws = 100)


# Create parameters ~ immune plot -----------------------------------------

## Model to data object
p_mod_1 <- as_draws_df(m_mod_1) |> 
  as.data.table()

## Select only principal effects
p_vars <- grep("Intercept|rescor|^lp|sigma", names(p_mod_1), value = TRUE, invert = TRUE)
p_mod_1 <- p_mod_1[, .SD, .SDcols = p_vars]

## Adapt data object to long format for ease of plotting
p_mod_1_long <- melt.data.table(p_mod_1, id.vars = c(".chain", ".iteration", ".draw"))

## Add identifying variables
p_mod_1_long[, `:=`(
  response = fcase(
    grepl("b_alpha_", variable), "alpha",
    grepl("b_beta_", variable), "beta",
    grepl("b_c_", variable), "c",
    grepl("b_lambda_", variable), "lambda",
    grepl("b_phi_", variable), "phi",
    grepl("b_tau_", variable), "tau",
    grepl("b_delta_", variable), "delta",
    default = NA
  ),
  variable = fcase(
    grepl("abc_linfocitos_total$", variable), "abc_linfocitos_total",
    grepl("abc_linfocitosb_total$", variable), "abc_linfocitosb_total",
    grepl("cd21_m_cd11c_p$", variable), "cd21_m_cd11c_p",
    grepl("cd21_m_cd11c_m$", variable), "cd21_m_cd11c_m",
    grepl("cd21_p_cd11c_m$", variable), "cd21_p_cd11c_m",
    grepl("cd21_p_cd11c_p$", variable), "cd21_p_cd11c_p",
    default = NA
  )
)]

## Change response to factor (so the labels appear ordered)
p_mod_1_long[
  j = response := factor(
    x = response, 
    levels = rev(c("alpha", "beta", "c", "lambda", "phi","tau", "delta"))
  )
]

## Make pretty labels for Immune markers
p_mod_1_long[
  j = variable2 := fcase(
    grepl("cd21_m_cd11c_p", variable), "CD21^-''~CD11c^+''",
    grepl("cd21_m_cd11c_m", variable), "CD21^-''~CD11c^-''",
    grepl("cd21_p_cd11c_m", variable), "CD21^+''~CD11c^-''",
    grepl("cd21_p_cd11c_p", variable), "CD21^+''~CD11c^+''",
    default = variable
  )
]

p_mod_1_long[, variable2 := factor(
  x = variable2, 
  levels = c(
    "CD21^+''~CD11c^-''",
    "CD21^+''~CD11c^+''",
    "CD21^-''~CD11c^-''",
    "CD21^-''~CD11c^+''"
  ))]

## Plot
fig_params_vs_immune <- ggplot(p_mod_1_long[!variable %in% c("abc_linfocitos_total", "abc_linfocitosb_total")], aes(value, response)) +
  facet_wrap(~ variable2, nrow = 2, labeller = label_parsed) +
  ggdist::stat_halfeye(aes(fill = after_stat(abs(x) > 0.1)), show.legend = FALSE) +
  theme_classic(base_size = 14, base_rect_size = 2/3, base_line_size = 2/3) +
  scale_y_discrete(labels = scales::label_parse()) +
  scale_fill_manual(values = c("grey30","darkorange")) +
  labs(x = "Standardized effect", y = "Model parameters",
       title = "Effect of Immune Markers",
       subtitle = "On Model Parameters Affecting RRi Dynamics",
       caption = "Unadjusted Effects")

ggsave(filename = "figures/fig_params_vs_immune_simple.pdf", 
       plot = fig_params_vs_immune, 
       width = 7, height = 7)

# Signatures based on Z-quantiles -----------------------------------------

## Function to generate predicted parameters based on input constrained data
generate_epred_params <- function(abc_linfocitos_total = 0,
                                  abc_linfocitosb_total = 0,
                                  cd21_m_cd11c_m = 0,
                                  cd21_m_cd11c_p = 0,
                                  cd21_p_cd11c_m = 0,
                                  cd21_p_cd11c_p = 0) {
  
  ## Using `epred_draws` function to generate predicted draws
  pred_data <- tidybayes::epred_draws(object = m_mod_1, 
                         ## Using expand grid to input sim data
                         newdata = expand.grid(
                           abc_linfocitos_total = abc_linfocitos_total,
                           abc_linfocitosb_total = abc_linfocitosb_total,
                           cd21_m_cd11c_m = cd21_m_cd11c_m,
                           cd21_m_cd11c_p = cd21_m_cd11c_p,
                           cd21_p_cd11c_m = cd21_p_cd11c_m,
                           cd21_p_cd11c_p = cd21_p_cd11c_p
                         )) |> 
    as.data.table()
  
  ## Filter columns with more than 1 unique value
  cols <- sapply(pred_data, function(i) length(x = unique(i)) > 1)
  pred_data <- pred_data[, .SD, .SDcols = cols]
  
  return(pred_data)
}

## Unstandardize expected parameters
unstandardize_epreds <- function(std_epred) {
  ## Obtain SD for each parameter
  sd_table <- m_data[, lapply(.SD, sd), .SDcols = alpha:delta] |> 
    melt.data.table(measure.vars = 1:7,
                    value.name = ".sd",
                    variable.name = ".category")
  
  ## Obtain Mean for each parameter
  mean_table <- m_data[, lapply(.SD, mean), .SDcols = alpha:delta] |> 
    melt.data.table(measure.vars = 1:7,
                    value.name = ".mean",
                    variable.name = ".category")
  
  ## Add those statistics to the predicted table per parameter
  std_epred <- std_epred[sd_table, on = ".category"]
  std_epred <- std_epred[mean_table, on = ".category"]
  
  ## Unstandardize each effect
  std_epred[, value := .mean + (.epred * .sd)][]
}

## Simulate RRi curve
sim_rri_curve <- function(unstd_epred) {
  
  ## Identify independent variable
  indp_var <- names(unstd_epred)[[1]]
  dcast_formula <- paste0(indp_var, " + .row + .draw ~ .category")
  
  ## Rearrange to wide format
  plot_data <- dcast.data.table(
    data = unstd_epred, 
    formula = as.formula(dcast_formula), 
    value.var = "value"
  )
  
  ## Generate a Time vector
  t_vec <- seq(0, 15, 0.1)
  
  ## Compute RRi Dynamics
  plot_data[
    j = list(
      time = t_vec,
      RRi = 
        ## Baseline
        alpha - 
        ## Exercise-induced Drop
        beta / (1 + exp(-lambda * (t_vec - tau))) +
        ## Recovery Kinetics
        (c * beta) / (1 + exp(-phi * (t_vec - tau - lambda)))
    ), 
    by = c(indp_var, ".draw")]
}

## Summarise the RRi curve per timestamp
summary_rri_curve <- function(simd_rri) {
  indp_var <- names(simd_rri)[[1]]
  out <- simd_rri[
    j = list(RRi = median(RRi), 
             ci_up = median(RRi) + mad(RRi), 
             ci_low = median(RRi) - mad(RRi)),
    by = c(indp_var, "time")]
  d_names <- names(out)
  d_names[[1L]] <- "pred"
  names(out) <- d_names
  return(out)
}

## Obtain predicted data (on standardized units)
pred_data <- list(
  cd21_m_cd11c_p = generate_epred_params(cd21_m_cd11c_p = c(-2,0,2)) |> 
    unstandardize_epreds() |>  sim_rri_curve() |>  summary_rri_curve(),
  cd21_m_cd11c_m = generate_epred_params(cd21_m_cd11c_m = c(-2,0,2)) |> 
    unstandardize_epreds() |>  sim_rri_curve() |>  summary_rri_curve(),
  cd21_p_cd11c_m = generate_epred_params(cd21_p_cd11c_m = c(-2,0,2)) |> 
    unstandardize_epreds() |>  sim_rri_curve() |>  summary_rri_curve(),
  cd21_p_cd11c_p = generate_epred_params(cd21_p_cd11c_p = c(-2,0,2)) |> 
    unstandardize_epreds() |>  sim_rri_curve() |>  summary_rri_curve()
) 

plot_data <- rbindlist(pred_data, idcol = "var")

plot_data[
  j = var2 := fcase(
    grepl("cd21_m_cd11c_p", var), "CD21^-''~CD11c^+''",
    grepl("cd21_m_cd11c_m", var), "CD21^-''~CD11c^-''",
    grepl("cd21_p_cd11c_m", var), "CD21^+''~CD11c^-''",
    grepl("cd21_p_cd11c_p", var), "CD21^+''~CD11c^+''",
    default = var
  )
]

plot_data[, var2 := factor(
  x = var2, 
  levels = c(
    "CD21^+''~CD11c^-''",
    "CD21^+''~CD11c^+''",
    "CD21^-''~CD11c^-''",
    "CD21^-''~CD11c^+''"
  ))]

fig_rri_signatures <- ggplot(plot_data, aes(time, RRi, group = pred)) +
  facet_wrap(~var2, labeller = label_parsed) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_up, fill = ordered(pred)), 
              alpha = 0.3) +
  geom_line(aes(col = ordered(pred)), linewidth = 1) +
  ggsci::scale_color_bmj(aesthetics = c("color", "fill")) +
  labs(x = "Time (min)", y = "RRi (ms)", 
       fill = "Cell count (Z-score)", 
       col = "Cell count (Z-score)",
       title = "Dependent RRi Signatures",
       subtitle = "On CD21 and CD11c Immune Cell Markers",
       caption = "Unadjusted Effects") +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(filename = "figures/fig_rri_signatures_simple.pdf",
       plot = fig_rri_signatures,
       width = 7, height = 9)
