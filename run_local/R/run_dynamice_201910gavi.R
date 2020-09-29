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
library (countrycode)
library (ggpubr)

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
  touchstone                        = "_201910gavi-5_",
  antigen                           = "measles-",
  vaccine_coverage_subfolder        = "scenarios/",
  
  # disease burden
  burden_template                   = "burden_template/central-burden-template.201910gavi-5.Measles_LSHTM-Jit_standard.csv",
  central_burden_estimate_folder    = "central_burden_estimate/",
  stochastic_burden_estimate_folder = "stochastic_burden_estimate/", 
  
  # diagnostic plots folder
  plot_folder                       = "plots/",
  
  # modelling group name
  group_name                        = "LSHTM-Jit-",
  
  # countries - specify iso3 codes to analyse only these countries
  #             or set it to "all" to analyse all included countries
  countries                         = c ("all"),
  # countries                         = c("AFG", "BDI"),  # debug -- c("BGD", "ETH") / "all"

  cluster_cores                     = 2,  # number of cores
  psa                               = 0   # psa runs; 0 for single run
  )

# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# specify scenarios to run
first_scenario <- 1
last_scenario  <- length (scenarios)
# debug
# first_scenario <- 1
# last_scenario  <- 1
# ------------------------------------------------------------------------------

# set SIAs and vaccination parameters for each scenario to minimize errors for running
set.sia         <- c (1, 0, 1, 0, 1, 0, 1, 0, 0, 1)
set.vaccination <- c (0, 2, 2, 1, 0, 2, 2, 1, 0, 2)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Create vaccine coverage file (0% coverage) for no vaccination scenario using 
# another vaccination scenario. This is done because VIMC vaccine coverage file
# for no vaccination scenario is empty
create_no_vaccination_coverage_file (
  no_vaccination_coverage_file = "vaccine_coverage/coverage_201910gavi-5_measles-no-vaccination.csv", 
  vaccination_coverage_file    = "vaccine_coverage/coverage_201910gavi-5_measles-mcv1-default.csv"
  )
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Create vaccine coverage file for campaign only vaccination scenario using 
# (routine + campaign) vaccination scenario. 
# Set routine coverage to zero. This is done because routine coverage values are
# needed even if they are only zero to run campaign only vaccination scenario.

# scenario: campaign-only-bestcase
create_campaign_vaccination_coverage_file (
  campaign_only_vaccination_coverage_file    = "vaccine_coverage/coverage_201910gavi-5_measles-campaign-only-bestcase.csv",
  routine_campaign_vaccination_coverage_file = "vaccine_coverage/coverage_201910gavi-5_measles-campaign-bestcase.csv"
  )
 
# scenario: campaign-only-default
create_campaign_vaccination_coverage_file (
  campaign_only_vaccination_coverage_file    = "vaccine_coverage/coverage_201910gavi-5_measles-campaign-only-default.csv",
  routine_campaign_vaccination_coverage_file = "vaccine_coverage/coverage_201910gavi-5_measles-campaign-default.csv"
)
# ------------------------------------------------------------------------------


# create remaining life expectancy file for each year across all age intervals
create_life_expectancy_remaining_full ()

# loop through scenarios
for (index in first_scenario:last_scenario) {

  # set scenario name 
  scenario_name   <- scenarios [index]
  print (scenario_name)
  tic ()
  
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
    vaccination                = set.vaccination [index],
    using_sia                  = set.sia         [index],
    measles_model              = "vaccine2019_sia_singlematrix.exe",
    debug_model                = FALSE
    )
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  # estimate deaths and DALYs
  #
  # apply CFR (case fatality rates) to estimate deaths -- Wolfson
  #   save results in corresponding cfr_option subfolder
  #   append cfr_option to results file
  # uncomment below to run with cfr_option = "Wolfson"
  # estimateDeathsDalys (cfr_option             = "Wolfson",
  #                      burden_estimate_file   = burden_estimate_file,
  #                      burden_estimate_folder = var$central_burden_estimate_folder)

  # apply CFR (case fatality rates) to estimate deaths -- Portnoy
  #   save results in corresponding cfr_option subfolder
  #   append cfr_option to results file
  estimateDeathsDalys (cfr_option             = "Portnoy",
                       burden_estimate_file   = burden_estimate_file,
                       burden_estimate_folder = var$central_burden_estimate_folder,
                       vimc_scenario          = scenario_name,
                       portnoy_scenario       = "s6"  # portnoy scenario 6
                       )
  # ----------------------------------------------------------------------------

  toc ()
} # end of loop -- for (scenario in scenarios)


# base scenario for comparison
base_scenario <- "no-vaccination"

# TEMP comment
# # ------------------------------------------------------------------------------
# # diagnostic plots of vaccine coverage and burden estimates (cases, deaths, dalys)
diagnostic_plots (
  vaccine_coverage_folder    = var$vaccine_coverage_folder,
  coverage_prefix            = var$coverage_prefix,
  touchstone                 = var$touchstone,
  antigen                    = var$antigen,
  scenarios                  = scenarios [first_scenario:last_scenario],
  base_scenario              = base_scenario,
  burden_estimate_folder     = var$central_burden_estimate_folder,
  plot_folder                = var$plot_folder,
  group_name                 = var$group_name,
  countries                  = var$countries,
  # cfr_options                = c("Wolfson", "Portnoy"),
  cfr_options                = c ("Portnoy"),
  psa                        = var$psa,
  start_year                 = 2000,
  end_year                   = 2100,
  compare_plots              = FALSE
  )
# # ------------------------------------------------------------------------------


# return to source directory
setwd (source_wd)

# end time
print (Sys.time ())

# ------------------------------------------------------------------------------
# main program -- end
# ------------------------------------------------------------------------------