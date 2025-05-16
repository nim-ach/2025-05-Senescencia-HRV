
# Prepare workspace -------------------------------------------------------

## Load libraries
library(CardioCurveR)

## Load previous data processing file
source("data-raw/1-basic-processing.r")

# Estimate RRi-vs-time parameters -----------------------------------------

## Estimate prior parameters for brms priors
start_params <- list(alpha = 0.85, beta = -3, c = 0.85, lambda = -2.5, phi = -2, tau = 6.6, delta = 2.5)
upper_lim <- list(alpha = 2, beta = -0.1, c = 2.0, lambda = -0.5, phi = -0.5, tau = 10, delta = 6)
lower_lim <- list(alpha = 0, beta = -5, c = 0.3, lambda = -6, phi = -6, tau = 3, delta = 1)

params_fit <- senescence_d[, {
  as.list(estimate_RRi_curve(time, RRi_std, start_params, lower_lim, upper_lim)$parameters)
}, id]

theta_each <- params_fit[, lapply(.SD[,-1L], median)]

t <- seq(0,15, 0.1)
rri_each <- dual_logistic(t, theta_each)

with(senescence_d, {
  plot(RRi_std ~ time, pch = 16, cex = 0.15, col = rgb(0,0,0,1))
  lines(t, rri_each, col = "white", lwd = 4)
})

descale_params <- senescence_d[, list(mean = unique(RRi_mean), sd = unique(RRi_sd)), id]

params_fit <- descale_params[params_fit, on = "id"]

params_fit <- params_fit[, list(alpha = alpha * sd + mean, 
                  beta = abs(beta * sd),
                  c, 
                  lambda = abs(lambda), 
                  phi = abs(phi), 
                  tau, delta), id]

senescence_d <- senescence_d[params_fit, on = "id"]

# Export the data ---------------------------------------------------------

## Save RData object
save(senescence_d, file = "data/senescence_d.RData")

## Save CSV file
fwrite(senescence_d, file = "data/senescence_d.csv")
