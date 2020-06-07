# run_measles_vimc.R
# main script for measles VIMC runs

# load libraries
library (tictoc)
library (data.table)
library (doParallel)
library (foreach)

# clear workspace
rm (list = ls ())

tic ()  # start timer
print (Sys.time ())

# global options in which R computes and displays its results.
options ("scipen" = 999, 
         "digits" = 14)

# ------------------------------------------------------------------------------
# source files containing functions

# create remaining life expectancy file for each year across all age intervals
source ("lexp.R")

# generate 2 vaccine coverage files per scenario for routine and SIA from VIMC vaccine coverage files
source ("generate_sia_vimc2019_0401.R")

# source functions -- setup
source ("gavi_setup_v201910version4_singlematrix.R")

# Run country funnction -- can be called in the foreach loop or in a for loop. 
source ("gavi_runcountry_v201910version4_singlematrix.R")

# Convert output files (burden estimates) provided by Allison with age-varying cfrs to vimc format
source ("allison_out_final.R")

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# create remaining life expectancy file for each year across all age intervals
create_life_expectancy_remaining_full ()

# generate 2 vaccine coverage files per scenario for routine and SIA from VIMC vaccine coverage files
#   input folder:  /generate_coverage/VIMC-coverage-files/
#   output folder: /input/scenarios2019/
#      note: vaccine coverage files
#            "...mcv2..." from vimc download is renamed to "...mcv12..."
generate_vaccine_coverage_routine_sia ()


# source gavi_v201910version4_det_local_scenario1_singlematrix.R
source ("gavi_v201910version4_det_local_scenario1_singlematrix _allison_in.R")

# Convert output files (burden estimates) provided by Allison with time-varying cfrs to vimc format
# Deaths and dalys will need to be computed using cases, cfrs and remaining life expectancy
create_central_burden_estimates ()

# ------------------------------------------------------------------------------

print (Sys.time ())
toc ()  # stop timer 
