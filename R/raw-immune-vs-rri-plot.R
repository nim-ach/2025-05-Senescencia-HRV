
# Prepare workspace -------------------------------------------------------

## Import libraries
library(data.table)
library(ggplot2)
library(brms)

## Import data 
data("senescence_d")

## Import model
m_mod_1 <- readRDS("models/m_mod_1.RDS")
m_mod_2 <- readRDS("models/m_mod_2.RDS")

# Prepare data objects ----------------------------------------------------

# Response variables
var_response <- c("alpha","beta","c","lambda","phi","tau","delta")

# Predictors
var_predictors <- c("abc_linfocitos_total", "abc_linfocitosb_total",
                    "cd21_m_cd11c_p","cd21_p_cd11c_p",
                    "cd21_p_cd11c_m","cd21_m_cd11c_m")

# Covariables
var_cov <- c("sex", "age", "weight", "height", "imc", "muscle_total", "fat_total")

var_all <- c(var_response, var_predictors, var_cov)

## Get model data, only one observation per subject
m_data <- senescence_d[, .SD[1], .SDcols = var_all, id]

## Standardize response variables
m_data_std <- datawizard::standardize(m_data)


# Plot signatures ---------------------------------------------------------

plot_signatures <- function(subpopulation = NULL) {
  
  j_var <- substitute({
    list(var = cut_number(subpopulation, n = 4), id)
  })
    
  groups <- senescence_d[, .SD[1], id][, eval(j_var)]
  
  lvls <- groups[, levels(var)[c(1,4)]]
  
  plot_data <- senescence_d[groups, on = "id"][var %in% lvls]
  
  plot_data[, var := factor(var, levels = unique(var), labels = c("Q1", "Q4"))]
  
  title_raw <- deparse(substitute(subpopulation))
  
  title_lab <- fcase(
    grepl("cd21_p_cd11c_p", title_raw), "CD21+ CD11c+",
    grepl("cd21_m_cd11c_p", title_raw), "CD21- CD11c+",
    grepl("cd21_p_cd11c_m", title_raw), "CD21+ CD11c-",
    grepl("cd21_m_cd11c_m", title_raw), "CD21- CD11c-"
  )
  
  ggplot(plot_data, aes(time, RRi)) +
    # facet_grid(cols = vars(var)) +
    geom_line(aes(group = id, col = var), alpha = 0.1, na.rm = TRUE) +
    geom_smooth(method = "gam",
                aes(group = var, col = var), na.rm = TRUE) +
    theme_classic() +
    labs(col = "Cell Count Quartile", y = "RRi (ms)", x = "Time (min)",
         title = title_lab) +
    ggsci::scale_color_lancet() +
    theme(legend.position = "top")
}

fig_raw <- ggpubr::ggarrange(
  plot_signatures(subpopulation = cd21_p_cd11c_m),
  plot_signatures(subpopulation = cd21_p_cd11c_p),
  plot_signatures(subpopulation = cd21_m_cd11c_m),
  plot_signatures(subpopulation = cd21_m_cd11c_p),
  ncol = 2, nrow = 2,
  common.legend = TRUE
)

plot(fig_raw)
