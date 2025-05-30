
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
                                      sex = NA, age = 0,
                                      muscle_total = 0, fat_total = 0
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
#> 1:      alpha   884.77 21.60 [840.27, 925.77]
#> 2:       beta   403.99 24.57  [355.6, 451.58]
#> 3:          c     0.89  0.02     [0.84, 0.94]
#> 4:     lambda     3.01  0.22     [2.57, 3.44]
#> 5:        phi     2.18  0.27     [1.63, 2.71]
#> 6:        tau     6.89  0.11     [6.66, 7.11]
#> 7:      delta     2.50  0.11     [2.28, 2.72]
