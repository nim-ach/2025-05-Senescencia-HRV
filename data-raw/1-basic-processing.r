
# Prepare workspace -------------------------------------------------------

## Load libraries
library(data.table)

## Import the immune data
data <- fread(file = "data-raw/raw-data.csv")

# Append the RRi data -----------------------------------------------------

## Standardize the identifying variables to bind the datasets
data[, `:=`(
  first_name = tools::toTitleCase(text = tolower(x = trimws(first_name))),
  last_name = tools::toTitleCase(text = tolower(x = trimws(last_name)))
)]

data[, name := paste(first_name, last_name)]

## Import RRi data ----

## Identify RRi files
rri_files <- list.files("data-raw/rri-data")
rri_paths <- paste0("data-raw/rri-data/", rri_files)

## And subject names
rri_names <- gsub("\\.txt", "", rri_files)
names(rri_paths) <- rri_names

## Then, import the data into a single dataset
rri_data <- lapply(rri_paths, CardioCurveR::import_RRi_txt, warn = FALSE) |> 
  rbindlist(idcol = "name")

## Finally, bind the RRi data to the immune data
data <- merge.data.table(data, rri_data, by = "name", all = TRUE)


# Data formatting ---------------------------------------------------------

## Remove HRV variables (Given we are working with RRi)
rm_vars <- grep("pre|peri|post", names(data), value = TRUE)
data[, (rm_vars) := NULL]

## Remove verification variables
rm_vars <- grep("complete", names(data), value = TRUE)
data[, (rm_vars) := NULL]

## Create ID variable based on name
data[, id := rleid(name)]

## Remove identifying variables
rm_vars <- grep("email|rut|phone|record|code|redcap|name", names(data), value = TRUE)
data[, (rm_vars) := NULL]

## Remove missing cases from the immune dataset
## 
## This means, keeping only the subjects from the
## immmune dataset that have RR-interval data
senescence_d <- data[!is.na(sex) & !is.na(RRi) & time <= 15]

## Sex as a factor
senescence_d[, sex := factor(sex)]

## Muscle as numeric variable
senescence_d[, muscle_left_arm := as.numeric(muscle_left_arm)]

## Create appendicular muscle mass (AMM) variable
senescence_d[, appendicular_muscle := 
               muscle_left_arm + muscle_right_arm + 
               muscle_left_leg + muscle_right_leg]

## Borg variable only have two unique values, so remove it
senescence_d[, borg_1 := NULL]

## Choose only individuals whom have at least 12 min of recording in the 
## two-minute step test
ind_id <- senescence_d[, diff(range(time)) >= 12, id][V1 == TRUE,id]
senescence_d <- senescence_d[id %in% ind_id]

## Add a standardized RRi for each subject to ease further modeling work
senescence_d[, RRi_std := scale(RRi), id]
senescence_d[, RRi_mean := mean(RRi, na.rm = TRUE), id]
senescence_d[, RRi_sd := sd(RRi, na.rm = TRUE), id]

## Remove percentage cells
senescence_d[, abc_linfocitos_porcentaje := NULL]
senescence_d[, abc_linfocitosb_porcentaje := NULL]
senescence_d[, cd21_p_cd11c_p_porcentaje := NULL]
senescence_d[, cd21_m_cd11c_m_porcentaje := NULL]
senescence_d[, cd21_p_cd11c_m_porcentaje := NULL]
senescence_d[, cd21_m_cd11c_p_porcentaje := NULL]

## Ratio between p_p and m_m 
senescence_d[, cd21_cd11c_p_m_ratio := cd21_p_cd11c_p / cd21_m_cd11c_m]

## Change subject with 57 to 67 (typo)
senescence_d[age < 60, age := age + 10]