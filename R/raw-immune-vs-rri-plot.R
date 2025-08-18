
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

plot_signatures <- function(subpopulation = NULL, d = senescence_d, return_data = FALSE) {
  
  j_var <- substitute({
    list(var = cut_number(subpopulation, n = 4), id)
  })
    
  groups <- d[, .SD[1], id][, eval(j_var)]
  
  lvls <- groups[, levels(var)[c(1,4)]]
  
  plot_data <- d[groups, on = "id"][var %in% lvls]
  
  plot_data[, var := factor(var, levels = unique(var), labels = c("Q1", "Q4"))][]
  
  if (isTRUE(return_data)) {
    return(plot_data)
  }
  
  title_raw <- deparse(substitute(subpopulation))
  
  title_lab <- fcase(
    grepl("cd21_p_cd11c_p", title_raw), "CD21+ CD11c+",
    grepl("cd21_m_cd11c_p", title_raw), "CD21- CD11c+",
    grepl("cd21_p_cd11c_m", title_raw), "CD21+ CD11c-",
    grepl("cd21_m_cd11c_m", title_raw), "CD21- CD11c-"
  )
  
  p <- ggplot(plot_data, aes(time, RRi)) +
    # facet_grid(cols = vars(var)) +
    geom_line(aes(group = id, col = var), alpha = 0.1, na.rm = TRUE) +
    geom_smooth(method = "gam",
                aes(group = var, col = var), na.rm = TRUE) +
    theme_classic() +
    labs(col = "Cell Count Quartile", y = "RRi (ms)", x = "Time (min)",
         title = title_lab) +
    ggsci::scale_color_lancet() +
    theme(legend.position = "top")
  
  return(p)
  
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

ggsave(filename = "figures/raw-immune-vs-rri-plot.pdf", 
       plot = fig_raw, width = 7, height = 7)


# -------------------------------------------------------------------------

id_table <- senescence_d[, .SD[1], id
             ][, list(
               id, sex, age,
               p_p = cut_number(cd21_p_cd11c_p, 4, labels = FALSE),
               m_m = cut_number(cd21_m_cd11c_m, 4, labels = FALSE)
               )][abs(p_p - m_m) > 2]

senescence_d[, .SD[1], id
][, list(
  id, sex, age,
  p_p = cut_number(cd21_p_cd11c_p, 4, labels = FALSE),
  m_m = cut_number(cd21_m_cd11c_m, 4, labels = FALSE)
)][abs(p_p - m_m) > 2]

try_id <- id_table$id

senescence_d[, list(RRi, time, sex, age, id)][id_table, on = "id"] |> 
  ggplot(aes(time, RRi)) +
  facet_wrap(~ sex) +
  geom_line(aes(group = id, col = ordered(p_p)), alpha = 0.5) +
  geom_smooth(aes(col = ordered(p_p)))

id_phase <- senescence_d[id %in% try_id, .SD[1], id][, list(id)]

fread(input = "data-raw/phase_rut.csv")[id_phase, on = "id"] |> 
  fwrite(file = "")

# -------------------------------------------------------------------------

p_data <- plot_signatures(subpopulation = cd21_p_cd11c_m, 
                d = senescence_d,
                return_data = TRUE)


try_id <- c(4, 27, 47, 48, 75, 91)

ggpubr::ggarrange(
  plot_signatures(subpopulation = cd21_p_cd11c_p, 
                  d = senescence_d[id %in% try_id]),
  plot_signatures(subpopulation = cd21_m_cd11c_m, 
                  d = senescence_d[id %in% try_id]),
  nrow = 2
)


senescence_d[id %in% c(4, 10, 27, 45, 47, 48, 75, 91)][, .SD[1], id]