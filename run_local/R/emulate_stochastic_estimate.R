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

# start time
print (Sys.time ())
# ------------------------------------------------------------------------------


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
# output: psa_vac_coverage_low_high_file.csv
#         iso3 | vac_coverage (low or high)
# ------------------------------------------------------------------------------
generate_cases_proportion <- function (psa_variables_filename,
                                       psa_variables_casesprop_filename,
                                       central_estimate_filename,
                                       vac_coverage_level,
                                       stochastic_estimate_single_country_filename,
                                       country_code,
                                       psa_vac_coverage_low_high_file,
                                       start_year, 
                                       end_year) {
  
  # read psa variables file
  psa_dt <- fread (psa_variables_filename)
  
  # loop through vaccination coverage levels (low and high)
  for (vac_coverage in vac_coverage_level) {
    
    # read stochastic estimate for a single country 
    stochastic_dt <- fread (stochastic_estimate_single_country_filename [vac_coverage])
    stochastic_dt <- stochastic_dt [year >= start_year & year <= end_year]
    
    # read central estimate
    central_dt <- fread (central_estimate_filename)
    central_dt <- central_dt [year >= start_year & year <= end_year]
    
    # extract central estimate for a single country
    central_dt <- central_dt [country == country_code [vac_coverage]]
    
    # total cases across all years in central estimate for a single country
    central_total_cases <- sum (central_dt [, cases])
    
    # total cases across all years in stochastic estimate for a single country for each run
    stochastic_dt <- stochastic_dt [, .(total_cases = sum (cases)), by = run_id]
    
    # add column for central total cases
    stochastic_dt [, central_total_cases := central_total_cases]
    
    # add column for proportional variation of cases' estimates for
    # stochastic runs in comparison to central run estimate
    stochastic_dt [, cases_prop := total_cases / central_total_cases]
    
    # drop un-required columns
    stochastic_dt [, c("total_cases", "central_total_cases") := NULL]
    
    # suffix cases_prop column with vaccination coverage level (low or high)
    setnames (x = stochastic_dt,
              old = "cases_prop",
              new = paste0 ("cases_prop", "_", vac_coverage) )
    
    # merge proportional variation of cases' estimates to psa variables table
    psa_dt <- psa_dt [stochastic_dt, 
                      on = .(run_id = run_id)]
  }
  
  # save file -- psa variables with cases proportion across runs
  fwrite (x    = psa_dt, 
          file = psa_variables_casesprop_filename)
  
  return (psa_dt)
  
} # end of function -- generate_cases_proportion
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# basic plots
# ------------------------------------------------------------------------------
basic_plots <- function (psa_dt) {
  
  hist (psa_dt$take1_input)
  hist (psa_dt$take2_input)
  hist (psa_dt$take3_input)
  
  hist(psa_dt$cases_prop_low)
  print (summary  (psa_dt$cases_prop_low))
  print (quantile (psa_dt$cases_prop_low, c(0.025, 0.5, 0.975)))
  
  plot (psa_dt$take1_input, psa_dt$cases_prop_low)
  plot (psa_dt$take2_input, psa_dt$cases_prop_low)
  plot (psa_dt$take3_input, psa_dt$cases_prop_low)
  
  hist(psa_dt$cases_prop_high)
  print (summary  (psa_dt$cases_prop_high))
  print (quantile (psa_dt$cases_prop_high, c(0.025, 0.5, 0.975)))
  
  plot (psa_dt$take1_input, psa_dt$cases_prop_high)
  plot (psa_dt$take2_input, psa_dt$cases_prop_high)
  plot (psa_dt$take3_input, psa_dt$cases_prop_high)
  
  return ()
  
} # end of function -- basic_plots
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# assign central estimate and stochastic estimate file names
# ------------------------------------------------------------------------------
assign_estimates_filenames <- function () {
  
  # scenarios
  scenarios <- c("campaign-only-bestcase",  # 1  SIAs only
                 "mcv2-bestcase",           # 2  MCV1 & MCV2
                 "campaign-bestcase",       # 3  MCV1 & MCV2 and SIAs
                 "mcv1-bestcase",           # 4  MCV1 only
                 "campaign-only-default",   # 5  SIAs only
                 "mcv2-default",            # 6  MCV1 & MCV2
                 "campaign-default",        # 7  MCV1 & MCV2 and SIAs
                 "mcv1-default",            # 8  MCV1 only
                 "no-vaccination",          # 9  no vaccination (set vaccination and using_sia to 0)
                 "stop"                     # 10 MCV1 & MCV2 and SIAs
                 )
  
  # flag -- vaccination or no vaccination scenario
  vaccine_flag <- c (T, T, T, T, T, T, T, T, F, T)
  
  central_estimate_folder <- "central_burden_estimate/201910gavi/"
  
  central_estimate_files <- paste0 (central_estimate_folder, 
                                    "central_burden_estimate_measles-LSHTM-Jit-", 
                                    scenarios, 
                                    "_Portnoy.csv")
  
  stochastic_estimate_folder <- "stochastic_burden_estimate/emulator/"
  
  stochastic_estimate_files <- paste0 (stochastic_estimate_folder,
                                       "stochastic_burden_estimate_measles-LSHTM-Jit-", 
                                       scenarios, 
                                       "_Portnoy.csv")
  
  return (list (central      = central_estimate_files, 
                stochastic   = stochastic_estimate_files, 
                vaccine_flag = vaccine_flag, 
                scenarios    = scenarios))
  
} # end of function -- assign_estimates_filenames
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# assign countries with low or high vaccination coverage based of MCV1 coverage in 2020
# low vaccination coverage < 85%
# high vaccination coverage >= 85%
# ------------------------------------------------------------------------------
set_countries_vac_coverage_level <- function (vac_coverage_file,
                                              psa_vac_coverage_low_high_file) {
  
  # read vaccine coverage file
  vac_coverage_level_dt <- fread (vac_coverage_file)
  
  # extract MCV1 vaccination coverage data for 2020
  vac_coverage_level_dt <- vac_coverage_level_dt [year == 2020 & vaccine == "MCV1", 
                                                  c("country_code", "coverage")]
  
  # assign countries with low or high vaccination coverage based of MCV1 coverage in 2020
  # low vaccination coverage < 85%
  # high vaccination coverage >= 85%
  vac_coverage_level_dt [coverage <  0.85, vac_coverage_level := "low"]
  vac_coverage_level_dt [coverage >= 0.85, vac_coverage_level := "high"]
  
  # drop un-required columns
  vac_coverage_level_dt [, coverage := NULL]
  
  # save file -- countries with low or high vaccination coverage based of MCV1 coverage in 2020
  fwrite (x    = vac_coverage_level_dt, 
          file = psa_vac_coverage_low_high_file)
  
  return (vac_coverage_level_dt)
  
} # end of function -- set_countries_vac_coverage_level
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# emulate stochastic estimates
# ------------------------------------------------------------------------------
emulate_stochastic_estimate <- function (central_file,
                                         stochastic_file,
                                         vaccine_flag,
                                         lexp_remain_file,
                                         psa_dt,
                                         vac_coverage_level_dt,
                                         countryCodes = -1,
                                         psa_runs) {
  
  # read in central burden results
  central_burden <- fread (central_file)
  
  # initialise empty table to save stochastic results in vimc format
  header <- central_burden [0, ]                # empty table with requisite 
                                                #   columns
  header [, run_id := numeric ()]               # add column for run id
  setcolorder (header, c("disease", "run_id"))  # reorder columns
  
  # save header to file for stochastic estimates of disease burden
  fwrite (x      = header, 
          file   = stochastic_file, 
          append = FALSE)
  
  # get country codes
  if (countryCodes [1] == -1) {
    countryCodes <- unique (central_burden [, country])
  }
  
  # ----------------------------------------------------------------------------
  # remaining life expectancy
  
  # read data for remaining life expectancy
  lexp_remain <- fread (lexp_remain_file)
  
  # add data column for remaining life expectancy
  central_burden <- lexp_remain [central_burden,
                           .(disease, year, age, i.country, country_name, cohort_size, cases, dalys, deaths, value),
                           on = .(country_code = country,
                                  age_from    <= age,
                                  age_to      >= age,
                                  year         = year) ]
  
  # rename columns
  setnames (central_burden, 
            old = c ("i.country", "value"),  
            new = c ("country",   "remain_lexp") )
  # ----------------------------------------------------------------------------
  
  # add data column for case fatality rates (cfr)
  central_burden [, cfr := 0]
  central_burden [cases > 0, cfr := deaths / cases]
  # cases in older years are 0 -- division by zero is avoided
  
  # ----------------------------------------------------------------------------
  
  # generate stochastic burden estimates for each country
  for (country_code in countryCodes) {
    
    # get central burden estimates for current country 
    central_burden_country <- central_burden [country == country_code]
    
    # iso3 country code 
    iso3_country_code <- country_code
    
    # check if country with low or high vaccination coverage (based of MCV1 coverage in 2020)
    vac_coverage_level <- vac_coverage_level_dt [country_code == iso3_country_code, 
                                                 vac_coverage_level]
    
    # initialise header for stochastic burden estimates
    stochastic_burden_country <- header
    
    # generate post-hoc stochastic estimates for each run
    for (run_number in 1:psa_runs) {
      
      # initialise burden estimate to central burden estimates for current country
      burden <- copy (central_burden_country)
      
      # add column for run id
      burden [, run_id := run_number]
      
      # get proportional variation of cases' estimates for this (stochastic) run(s)
      if (vaccine_flag) {  # vaccination scenario
        
        if (vac_coverage_level == "low") {
          cases_prop <- psa_dt [run_id == run_number, 
                                cases_prop_low]
          
        } else if (vac_coverage_level == "high") {
          cases_prop <- psa_dt [run_id == run_number, 
                                cases_prop_high]
        }
        
      } else {  # no vaccination scenario
        
        # no proportional variation of cases since they vary only by vaccine efficacy
        cases_prop <- 1  
      }
      
      # get proportional change in cfr (case fatality rate)
      cfr_prop <- psa_dt [run_id == run_number, 
                          mortality_input]
      
      # ------------------------------------------------------------------------
      # sensitivity analysis
      # estimate proportional variation in cases, deaths, and dalys
      
      # CASES
      # proportionally change the cases' estimates for current run for current country
      burden [, cases := cases * cases_prop]
      
      # DEATHS
      # estimate deaths -- apply cfr proportion for current run for current country
      # proportionally change the deaths for current run for current country
      burden [, deaths := cases * cfr * cfr_prop]
      
      # DALYs
      # calculate dalys = (ylds) + (ylls)
      burden [, dalys := ((cases - deaths) * 0.002) + (deaths * remain_lexp)]
      
      # drop columns -- cfr and remaining life expectancy
      burden [, c("cfr", "remain_lexp") := NULL]
      # ------------------------------------------------------------------------
      
      # add stochastic burden estimate adjusted by sensitivity analysis to the
      # stochastic burden table of all runs
      stochastic_burden_country <- rbind (stochastic_burden_country, 
                                          burden,
                                          use.names = TRUE)
      
    } # end of -- for (run_number in 1:psa_runs)
    
    # save stochastic burden estimates of current country to larger file with 
    # stochastic burden estimates of all countries
    fwrite (x      = stochastic_burden_country, 
            file   = stochastic_file, 
            append = TRUE)
    
  } # end of -- for (country_code in countryCodes) 
  
  return ()
  
} # end of function -- emulate_stochastic_estimate
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# main program -- start
# ------------------------------------------------------------------------------

# move to base directory (run code from source directory)
source_wd <- getwd ()
setwd ("../")

# ------------------------------------------------------------------------------
# file names
psa_variables_filename           <- "input/psa_variables.csv"
psa_variables_casesprop_filename <- "input/psa_variables_casesprop.csv"

# vaccination coverage labels
vac_coverage_level <- c(low = "low", high = "high")

# single country based of low and high vaccination coverage
country_code <- c ("NGA", "IND")
names (country_code) <- vac_coverage_level

# burden estimates
central_estimate_filename    <- "central_burden_estimate/201910gavi/central_burden_estimate_measles-LSHTM-Jit-campaign-default_Portnoy.csv"
stochastic_estimate_single_country_filename <- 
  c("stochastic_burden_estimate/Portnoy/stochastic_burden_estimate_measles-LSHTM-Jit-campaign-default_Portnoy_NGA.csv", 
    "stochastic_burden_estimate/Portnoy/stochastic_burden_estimate_measles-LSHTM-Jit-campaign-default_Portnoy_IND.csv")
names (stochastic_estimate_single_country_filename) <- vac_coverage_level

# vaccination coverae
vac_coverage_file <- "vaccine_coverage/coverage_201910gavi-5_measles-campaign-default.csv"

# remaining life expectancy file
lexp_remain_file <- "input/demographicdata2019/201910gavi-4_dds-201910_2_life_ex_both_full.csv"

# countries with low or high vaccination coverage based of MCV1 coverage in 2020
# low vaccination coverage < 85%
# high vaccination coverage >= 85%
psa_vac_coverage_low_high_file <- "input/psa_vac_coverage_low_high.csv"
# ------------------------------------------------------------------------------

# generate proportional variation of cases' estimates for the stochastic runs
psa_dt <- generate_cases_proportion (
  psa_variables_filename                      = psa_variables_filename,
  psa_variables_casesprop_filename            = psa_variables_casesprop_filename,
  central_estimate_filename                   = central_estimate_filename,
  vac_coverage_level                          = vac_coverage_level,
  stochastic_estimate_single_country_filename = stochastic_estimate_single_country_filename,
  country_code                                = country_code,
  start_year                                  = 2000,
  end_year                                    = 2030)


# basic plots
# basic_plots (psa_dt = psa_dt)

# assign countries with low or high vaccination coverage based of MCV1 coverage in 2020
# low vaccination coverage < 85%
# high vaccination coverage >= 85% # CONTINUE
vac_coverage_level_dt <- set_countries_vac_coverage_level (
  vac_coverage_file              = vac_coverage_file,
  psa_vac_coverage_low_high_file = psa_vac_coverage_low_high_file)
# ------------------------------------------------------------------------------

# assign central estimate and stochastic estimate file names
estimate_files <- assign_estimates_filenames ()

# loop through scenarios
for (i in 1:length(estimate_files$central)) {
  
  tic ()
  
  # emulate stochastic estimates
  emulate_stochastic_estimate (
    central_file          = estimate_files$central [i],
    stochastic_file       = estimate_files$stochastic [i],
    vaccine_flag          = estimate_files$vaccine_flag [i], 
    lexp_remain_file      = lexp_remain_file, 
    psa_dt                = psa_dt,
    vac_coverage_level_dt = vac_coverage_level_dt,
    countryCodes          = -1,  # all countries: -1; specific: c("AFG", "IND")
    psa_runs              = 200)
  
  toc ()
}

# return to source directory
setwd (source_wd)

# end time
print (Sys.time ())

# ------------------------------------------------------------------------------
# main program -- end
# ------------------------------------------------------------------------------
