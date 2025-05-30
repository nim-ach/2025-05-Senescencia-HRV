
# Prepare workspace -------------------------------------------------------

## Load libraries
library(data.table)
library(brms)
library(ggplot2)

## Load data
data("senescence_d")

senescence_d[, plot(RRi_std ~ time, pch = ".", col = rgb(0,0,0,.1))]

# Prepare data objects ----------------------------------------------------

std_covariates <- senescence_d[, .SD[1], id
             ][, j = datawizard::standardise(.SD, exclude = "id"), 
               .SDcols = c("fat_total","muscle_total","sex","age","id",
                           "abc_linfocitos_total", "abc_linfocitosb_total",
                                  "cd21_m_cd11c_p", "cd21_m_cd11c_m",
                                  "cd21_p_cd11c_m", "cd21_p_cd11c_p")]

m_data_std <- senescence_d[, list(RRi_std, RRi_mean, RRi_sd, time, id)
                           ][std_covariates, on = "id"]


# Intercept only model ----------------------------------------------------

i_formula <- bf(RRi_std ~ 
                  alpha - beta / (1 + exp(-lambda * (time - tau))) +
                  (c * beta) / (1 + exp(-phi * (time - tau - delta))),
                alpha + beta + c + lambda + phi + tau + delta ~ 1 + (1|id),
                nl = TRUE)

default_prior(i_formula, m_data_std)

nl_par <- c("alpha", "beta", "c", "lambda", "phi", "tau", "delta")

i_prior <- c(
  set_prior("normal(1, 2)", class = "sigma", lb = 0),
  set_prior("normal(.85, 1.0)", coef = "Intercept", nlpar = "alpha"),
  set_prior("normal(3.0, 0.5)", coef = "Intercept", nlpar = "beta"),
  set_prior("normal(.85, 0.25)", coef = "Intercept", nlpar = "c"),
  set_prior("normal(2.0, 0.5)", coef = "Intercept", nlpar = "lambda"),
  set_prior("normal(2.0, 0.5)", coef = "Intercept", nlpar = "phi"),
  set_prior("normal(6, 1)", coef = "Intercept", nlpar = "tau"),
  set_prior("normal(3, 1)", coef = "Intercept", nlpar = "delta"),
  set_prior("normal(0, 0.5)", class = "sd", nlpar = "alpha", lb = 0),
  set_prior("normal(0, 0.25)", class = "sd", nlpar = "beta", lb = 0),
  set_prior("normal(0, 0.125)", class = "sd", nlpar = "c", lb = 0),
  set_prior("normal(0, 0.25)", class = "sd", nlpar = "lambda", lb = 0),
  set_prior("normal(0, 0.25)", class = "sd", nlpar = "phi", lb = 0),
  set_prior("normal(0, 0.5)", class = "sd", nlpar = "tau", lb = 0),
  set_prior("normal(0, 0.5)", class = "sd", nlpar = "delta", lb = 0)
)

m_mod_int_3 <- brm(
  formula = i_formula,
  data = m_data_std,
  family = gaussian(),
  prior = i_prior,
  chains = 5, cores = 5,
  iter = 4000, warmup = 2000,
  seed = 1234,
  file = "models/m_mod_int_3.RDS"
)

conditional_effects(m_mod_int_3)

# Joint Hierarchical Model ------------------------------------------------

h_formula <- bf(RRi_std ~ 
     alpha - beta / (1 + exp(-lambda * (time - tau))) +
     (c * beta) / (1 + exp(-phi * (time - tau - delta))),
   alpha + beta + c + lambda + phi + tau + delta ~ 
     abc_linfocitos_total + abc_linfocitosb_total + cd21_m_cd11c_p + 
     cd21_m_cd11c_m + cd21_p_cd11c_m + cd21_p_cd11c_p + 
     sex + age + fat_total + muscle_total + (1|id),
   nl = TRUE)

default_prior(h_formula, m_data_std)

nl_par <- c("alpha", "beta", "c", "lambda", "phi", "tau", "delta")

h_prior <- c(
  set_prior("normal(1, 2)", class = "sigma", lb = 0),
  set_prior("normal(0, 0.1)", class = "b", nlpar = "alpha"),
  set_prior("normal(0, 0.1)", class = "b", nlpar = "beta"),
  set_prior("normal(0, 0.1)", class = "b", nlpar = "c"),
  set_prior("normal(0, 0.1)", class = "b", nlpar = "lambda"),
  set_prior("normal(0, 0.1)", class = "b", nlpar = "phi"),
  set_prior("normal(0, 0.1)", class = "b", nlpar = "tau"),
  set_prior("normal(0, 0.1)", class = "b", nlpar = "delta"),
  set_prior("normal(.85, 0.1)", coef = "Intercept", nlpar = "alpha"),
  set_prior("normal(3.0, 0.1)", coef = "Intercept", nlpar = "beta"),
  set_prior("normal(.85, 0.1)", coef = "Intercept", nlpar = "c"),
  set_prior("normal(2.5, 0.1)", coef = "Intercept", nlpar = "lambda"),
  set_prior("normal(2.0, 0.1)", coef = "Intercept", nlpar = "phi"),
  set_prior("normal(6.6, 0.1)", coef = "Intercept", nlpar = "tau"),
  set_prior("normal(2.5, 0.1)", coef = "Intercept", nlpar = "delta")
  , set_prior("normal(0, 0.1)", class = "sd", nlpar = "alpha", lb = 0)
  , set_prior("normal(0, 0.1)", class = "sd", nlpar = "beta", lb = 0)
  , set_prior("normal(0, 0.1)", class = "sd", nlpar = "c", lb = 0)
  , set_prior("normal(0, 0.1)", class = "sd", nlpar = "lambda", lb = 0)
  , set_prior("normal(0, 0.1)", class = "sd", nlpar = "phi", lb = 0)
  , set_prior("normal(0, 0.1)", class = "sd", nlpar = "tau", lb = 0)
  , set_prior("normal(0, 0.1)", class = "sd", nlpar = "delta", lb = 0)
)

m_mod_3 <- brm(
  formula = h_formula,
  data = m_data_std,
  family = gaussian(),
  prior = h_prior,
  chains = 5, cores = 5,
  iter = 2000, warmup = 1000,
  seed = 1234,
  #file = "models/m_mod_3.RDS",
  sample_prior = "only"
)

conditional_effects(m_mod_3)
