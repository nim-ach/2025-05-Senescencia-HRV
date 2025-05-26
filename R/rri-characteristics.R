
# Prepare workspace -------------------------------------------------------

## Import libraries
library(data.table)
library(brms)

## Import custom functions
source("R/_functions.R")

## Import data
data("senescence_d")

## Import adjusted fitted model
m_mod_2 <- readRDS("models/m_mod_2.RDS")

## And non-standardized data
m_data <- readRDS("models/m_data.RDS")

## Then, estimate mean and sd tables to transform standardized
## model parameters back to unstandardized parameters for reporting
## purposes
## 1) Obtain SD for each parameter
sd_table <- m_data[, lapply(.SD, sd), .SDcols = alpha:delta] |> 
  melt.data.table(measure.vars = 1:7,
                  value.name = ".sd",
                  variable.name = ".category")

## 2) Obtain Mean for each parameter
mean_table <- m_data[, lapply(.SD, mean), .SDcols = alpha:delta] |> 
  melt.data.table(measure.vars = 1:7,
                  value.name = ".mean",
                  variable.name = ".category")
## 3) Merge both tables into one
normal_table <- mean_table[sd_table, on = ".category"]
## 4) Remove residual data objects
rm(mean_table, sd_table)

# Extract predicted parameters --------------------------------------------

## Using `epred_draws` function to generate predicted draws
pred_data <- tidybayes::epred_draws(object = m_mod_2, 
                                    ## Using expand grid to input sim data
                                    newdata = expand.grid(
                                      abc_linfocitos_total = 0,
                                      abc_linfocitosb_total = 0,
                                      cd21_m_cd11c_m = 0,
                                      cd21_m_cd11c_p = 0,
                                      cd21_p_cd11c_m = 0,
                                      cd21_p_cd11c_p = 0,
                                      sex = NA, age = 0
                                    ),
                                    allow_new_levels = TRUE) |> 
  as.data.table()

## Filter columns with more than 1 unique value
cols <- sapply(pred_data, function(i) length(x = unique(i)) > 1)
pred_data <- pred_data[, .SD, .SDcols = cols]

## Remove intermediate `cols` object
rm(cols)

## Compute unstandardized model parameters
pred_data <- 
  pred_data[normal_table, on = ".category"
            ][, value := .mean + (.epred * .sd)][]

pred_data[, list(Estimate = round(x = median(value), digits = 2),
                 SE = round(x = mad(value), digits = 2),
                 CI = paste0("[", paste0(
                   round(x = tidybayes::hdi(value), digits = 2), 
                   collapse = ", "), "]")), 
          list(Parameters = .category)]
#>    Parameters Estimate    SE               CI
#>        <fctr>    <num> <num>           <char>
#> 1:      alpha   876.89 15.29 [845.45, 906.55]
#> 2:       beta   386.64 17.99 [350.97, 420.27]
#> 3:          c     0.90  0.02     [0.86, 0.93]
#> 4:     lambda     2.96  0.16     [2.63, 3.26]
#> 5:        phi     2.29  0.20     [1.91, 2.69]
#> 6:        tau     6.75  0.08     [6.58, 6.91]
#> 7:      delta     2.47  0.08      [2.3, 2.62]
