# emulate_stochastic_estimate.R
# emulate stochastic estimates (short-circuit approximation method)

# ------------------------------------------------------------------------------
# load libraries
library (tictoc)
library (data.table)
library (stringr)
library (ggplot2)
library (scales)
library (doParallel)
library (foreach)
library (countrycode)
library (ggpubr)
library (lhs) 
library (truncnorm)

# remove all objects from workspace
rm (list = ls ())
# ------------------------------------------------------------------------------

# move to base directory (run code from source directory)
source_wd <- getwd ()
setwd ("../")

# ------------------------------------------------------------------------------
# generate proportional variation of cases' estimates for the stochastic runs  
# based on a stochastic estimate for a single country (India) and central 
# estimate for a single country (India) for a single scenario 
# "campaign-default" # MCV1 & MCV2 and SIAs
# input: central estimate file + stochastic estimate file for a single country
# input:  psa_variables.csv
#         run_id | take1_input |	take2_input |	take3_input |	mortality_input
# output: psa_variables_casesprop.csv 
#         run_id | take1_input |	take2_input |	take3_input |	mortality_input | cases_prop
# ------------------------------------------------------------------------------
generate_cases_proportion <- function () {
  
  # file names
  psa_variables_filename           <- "input/psa_variables.csv"
  psa_variables_casesprop_filename <- "input/psa_variables_casesprop.csv"
  
  central_estimate_filename    <- "central_burden_estimate/201910gavi/central_burden_estimate_measles-LSHTM-Jit-campaign-default_Portnoy.csv"
  stochastic_estimate_filename <- "stochastic_burden_estimate/Portnoy/stochastic_burden_estimate_measles-LSHTM-Jit-campaign-default_Portnoy.csv"
  country_code                 <- "IND"  # single country
  
  # start and end years
  start_year <- 2000
  # end_year   <- 2100
  end_year   <- 2030
  
  # read psa variables file
  psa_dt <- fread (psa_variables_filename)
  
  # read central estimate
  central_dt <- fread (central_estimate_filename)
  central_dt <- central_dt [year >= start_year & year <= end_year]
  
  # read stochastic estimate for a single country 
  stochastic_dt <- fread (stochastic_estimate_filename)
  stochastic_dt <- stochastic_dt [year >= start_year & year <= end_year]
  
  # extract central estimate for a single country
  central_dt <- central_dt [country == country_code,]
  
  # total cases across all years in central estimate for a single country
  central_total_cases <- sum (central_dt [, cases])
  
  # total cases across all years in stochastic estimate for a single country for each run
  stochastic_dt <- stochastic_dt [, .(total_cases = sum (cases)), by = run_id]
  
  # add column for central total cases
  stochastic_dt [, central_total_cases := central_total_cases]
  
  # add column for proportional variation of cases' estimates for
  # stochastic runs in comparison to central run estimate
  stochastic_dt [, cases_prop := total_cases / central_total_cases]
  
  # merge proportional variation of cases' estimates to psa variables table
  psa_dt <- psa_dt [stochastic_dt, 
                    on = .(run_id = run_id)]
  
  # drop un-required columns
   psa_dt [, c("total_cases", "central_total_cases") := NULL]
  
  # save file -- psa variables with cases proportion across runs
  fwrite (x    = psa_dt, 
          file = psa_variables_casesprop_filename)
  
  return (psa_dt)
  
}
# ------------------------------------------------------------------------------


hist     (psa_dt$cases_prop)
summary  (psa_dt$cases_prop)
quantile (psa_dt$cases_prop, c(0.025, 0.5, 0.975))

plot (psa_dt$take1_input, psa_dt$cases_prop)
plot (psa_dt$take2_input, psa_dt$cases_prop)
plot (psa_dt$take3_input, psa_dt$cases_prop)

# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------

# return to source directory
setwd (source_wd)

