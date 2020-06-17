# run_dynamice.R
# run it from the source directory

# DynaMICE (Dynamic Measles Immunisation Calculation Engine)
# Measles vaccine impact model
#   To estimate the health impact of measles vaccination for a given set of 
#   countries.

# ------------------------------------------------------------------------------
# load libraries
library (tictoc)
library (data.table)
library (stringr)
library (ggplot2)
library (scales)
library (doParallel)
library (foreach)

# remove all objects from workspace
rm (list = ls ())
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# include functions
source ("functions.R")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# main program -- start
# ------------------------------------------------------------------------------

# start time
print (Sys.time ())
tic ()

# set seed for random number generator
set.seed (1) 

# move to base directory (run code from source directory)
source_wd <- getwd ()
setwd ("../")

# set values for global variables
var <- list (
  # vaccine coverage
  vaccine_coverage_folder           = "vaccine_coverage/",
  coverage_prefix                   = "coverage", 
  touchstone                        = "_202005covid-4_",
  antigen                           = "measles-",
  vaccine_coverage_subfolder        = "scenarios/",
  
  # disease burden
  burden_template                   = "burden_template/central-burden-template.202005covid-4.Measles_LSHTM-Jit_standard.csv",
  central_burden_estimate_folder    = "central_burden_estimate/",
  stochastic_burden_estimate_folder = "stochastic_burden_estimate/", 
  
  # diagnostic plots folder
  plot_folder                       = "plots/",
  
  # modelling group name
  group_name                        = "LSHTM-Jit-",
  
  # countries - specify iso3 codes to analyse only these countries
  #             or set it to "all" to analyse all included countries 
  countries                         = c("all"),  # debug -- c("ETH"), 
  
  cluster_cores                     = 2,  # number of cores
  psa                               = 0   # psa runs; 0 for single run
  )

# scenarios
scenarios <- c("counterfactual-bau-scenario1", 
               "disruption-scenario2-sia2021", 
               "disruption-scenario3-sia2022", 
               "disruption-scenario4-25rout",
               "disruption-scenario5-25rout-sia2021",
               "disruption-scenario6-25rout-sia2022",
               "disruption-scenario7-50rout",
               "disruption-scenario8-50rout-sia2021",
               "disruption-scenario9-50rout-sia2022",
               "disruption-scenario10-25rout-25sia"
               )

# debug
# scenarios <- scenarios [1]

# create remaining life expectancy file for each year across all age intervals
create_life_expectancy_remaining_full ()

# loop through scenarios
for (index in 1:length(scenarios)) {

  # set scenario name 
  scenario_name   <- scenarios [index]
  
  # set scenario number in order to save it in the right folder
  #   scenarioXX, XX = 01, 02, ... ,10, 11, 12, ...
  #   this is set to 10 characters in the fortran model code
  scenario_number <- sprintf ("scenario%02d", index)
  
  # ----------------------------------------------------------------------------
  # generate 2 vaccine coverage files per scenario for routine and SIA from 
  # VIMC vaccine coverage file: 
  #   paste0 (vaccine_coverage_folder, prefix, touchstone, scenario, ".csv")
  create_vaccine_coverage_routine_sia (
    vaccine_coverage_folder    = var$vaccine_coverage_folder, 
    coverage_prefix            = var$coverage_prefix,
    touchstone                 = var$touchstone,
    antigen                    = var$antigen,
    scenario_name              = scenario_name,
    vaccine_coverage_subfolder = var$vaccine_coverage_subfolder
    )
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  # estimate cases
  #
  # run scenario -- get burden estimates -- primarily cases
  # return burden estimate file name where estimates are saved
  burden_estimate_file <- runScenario (
    vaccine_coverage_folder    = var$vaccine_coverage_folder,
    coverage_prefix            = var$coverage_prefix,
    touchstone                 = var$touchstone,
    antigen                    = var$antigen,
    scenario_name              = scenario_name,
    scenario_number            = scenario_number,
    vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
    burden_template            = var$burden_template,
    burden_estimate_folder     = var$central_burden_estimate_folder,
    group_name                 = var$group_name,
    countries                  = var$countries,
    cluster_cores              = var$cluster_cores,
    psa                        = var$psa,
    vaccination                = 2,
    using_sia                  = 1 
    )
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  # estimate deaths and DALYs
  #
  # apply CFR (case fatality rates) to estimate deaths -- Wolfson
  #   save results in corresponding cfr_option subfolder
  #   append cfr_option to results file
  estimateDeathsDalys (cfr_option             = "Wolfson",
                       burden_estimate_file   = burden_estimate_file,
                       burden_estimate_folder = var$central_burden_estimate_folder)

  # apply CFR (case fatality rates) to estimate deaths -- Portnoy
  #   save results in corresponding cfr_option subfolder
  #   append cfr_option to results file
  estimateDeathsDalys (cfr_option             = "Portnoy",
                       burden_estimate_file   = burden_estimate_file,
                       burden_estimate_folder = var$central_burden_estimate_folder)
  # ----------------------------------------------------------------------------

} # end of loop -- for (scenario in scenarios)

# ------------------------------------------------------------------------------
# diagnostic plots of vaccine coverage and burden estimates (cases, deaths, dalys)
diagnostic_plots (
  vaccine_coverage_folder    = var$vaccine_coverage_folder,
  coverage_prefix            = var$coverage_prefix,
  touchstone                 = var$touchstone,
  antigen                    = var$antigen,
  scenarios                  = scenarios,
  burden_estimate_folder     = var$central_burden_estimate_folder,
  plot_folder                = var$plot_folder,
  group_name                 = var$group_name,
  countries                  = var$countries,
  cfr_options                = c("Wolfson", "Portnoy"),
  psa                        = var$psa,
  start_year                 = 2015,
  end_year                   = 2030
)
# ------------------------------------------------------------------------------


# return to source directory
setwd (source_wd)

# end time
print (Sys.time ())
toc ()

# ------------------------------------------------------------------------------
# main program -- stop
# ------------------------------------------------------------------------------