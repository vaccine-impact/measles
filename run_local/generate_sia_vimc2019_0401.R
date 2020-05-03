# Title:       Generate coverage_routine and coverage_sia
# Author:      Petra Klepac, London School of Hygiene & Tropical Medicine
# Date:        22/11/2019
# Description: Read in scenario files and format them for DynaMICE model.
# Data:        needs campaign (SIA) data:  https://montagu.vaccineimpact.org 
#              >> Scenario >> Download combined coverage set data
# 
# Change - there are many problems with age_range_verbatim column, so following WHO 2019 estimates
#          this one uses only age_first and age_last columns to determine the age range of SIA campaign


# generate 2 vaccine coverage files per scenario for routine and SIA from 
# VIMC vaccine coverage files
#   note: vaccine coverage files
#     "...mcv2..." from vimc download is renamed to "...mcv12..."
generate_vaccine_coverage_routine_sia <- function () {
  
  # set working directory 
  wd <- getwd ()
  
  # total number of scenarios
  scenarios <- 10
  
  # loop through different scenarios
  for (index in 1:scenarios) {
    
    # scenarios
    scenario.name <- c("campaign-only-bestcase",           # 1  SIAs only
                       "mcv12-bestcase",                   # 2  MCV1, MCV2        # "...mcv2..." from vimc download is renamed to "...mcv12..."
                       "campaign-bestcase",                # 3  MCV1, MCV2, SIAs
                       "mcv1-bestcase",                    # 4  MCV1 only
                       "campaign-only-default",            # 5  SIAs only
                       "mcv12-default",                    # 6  MCV1, MCV2        # "...mcv2..." from vimc download is renamed to "...mcv12..."
                       "campaign-default",                 # 7  MCV1, MCV2, SIAs
                       "mcv1-default",                     # 8  MCV1 only
                       "no-vaccination",                   # 9  no vaccination (set vaccination and using_sia to 0)
                       "stop" )                            # 10 MCV1, MCV2, SIAs 
    
    # specific scenario
    scenario <- scenario.name [index]
    
    # index variable for scenario to save it to a correct folder (scenarioXX, XX = 01,02,..... ,10)
    counter <- c(paste0 ("0", c(1:9)), "10")
    save.scenario <- paste0("scenario", counter[index])
    
    # vaccination coverage filename
    fname      <- paste0 (wd, "/generate_coverage/VIMC-coverage-files/")
    input_sia  <- paste0 (fname, "coverage_201910gavi-4_measles-", scenario, ".csv" ) # full path to coverage file
    output_dir <- paste0 (wd, "/input/scenarios2019/") 
    
    # read data
    sia <- fread (input_sia, stringsAsFactors = F, na.string = "<NA>")
    
    # extract routine vaccination out of those
    routine <- sia[activity_type != "campaign",]
    
    # this line is redundant -- vaccine is already set to MCV1 or MCV2 in the corresponding rows
    routine [, vaccine := substr (set_name, 10, 13)]  
    
    # ----------------------------------------------------------------------------
    # only select campaigns
    sia2 <- sia [activity_type == "campaign" & ( (!is.na(target) & target != 0) | (!is.na(coverage) ) ), ]
    
    # why sias with zero coverage are included -- shouldn't this be (check this)
    # sia2 <- sia [activity_type == "campaign" & ( (!is.na(target) & target != 0) | (!is.na(coverage) & coverage != 0) ), ]
    # ----------------------------------------------------------------------------
    
    # remove unused columns
    discard.cols <- c ("scenario", "set_name", "gavi_support", "activity_type")
    sia2 <- sia2 [, .SD, .SDcols = -discard.cols]
    
    # remove all whitespace from age_range_verbatim, put all in lowercase
    sia2$age_range_verbatim <- tolower (gsub("[[:space:]]", "", sia2$age_range_verbatim))
    
    # use more inprecise values if age range verbatim is empty
    temp.cols2 <- c("t_combined", "t_l_a0", "t_l_a1", "t_n_a0", "t_n_a1")
    
    # age_first and age_last is what is used
    sia2 [, (temp.cols2) := list (0, "y", "y", age_first, age_last)]
    
    # calculate values for a0 and a1
    # for age <= 3 years, age is in weeks
    # for age >= 3 years, (156 weeks for upto age 3) + (remaining years above 3)
    #      age 1 year ~ 52; 2 year ~ 104; 3 year ~ 156; 4 year ~ 157; 5 year ~ 158
    find.a <- function (x) {
      t0 <- ifelse (x <= 3, round (x * 52), round ((x - 3) + (3 * 52) ) )
      return (t0)
    }
    
    # set start and end ages for SIA
    sia2 [, `:=` (a0 = find.a(t_n_a0), a1 = find.a(t_n_a1))]
    
    # age 0 ~ 1 week
    sia2 [a0 == 0, a0 := 1]  
    sia2 [a1 == 0, a1 := 1]
    
    # not sure what next 2 lines are for? (check)
    sia2 [t_l_a0 != "y",]
    sia2 [t_l_a1 != "y",]
    
    # remove temporary columns
    sia2 <- sia2 [, .SD, .SDcols = -temp.cols2]
    
    # FROM GUIDELINES FOR 2019 RUNS: 
    # *coverage* shows the level of vaccination coverage, usually ranging from 0 (0%)
    # to 1 (100%). In some cases, this value may be greater than 1. For example, if a
    # campaign originally targets 1 million people but ends up vaccinating 1.1
    # million people, the coverage would be shown as 1.1 (equating to 110%). Coverage
    # and target population are now always specified at a national level. For
    # example, where a campaign targets all ages in Region A (population 1,000,000)
    # and achieves 90% coverage, and where the population of the whole country is
    # 5,000,000, the coverage would appear on Montagu as 0.18 (18%) and the target
    # population as 5,000,000. This way of specifying coverage and target population
    # applies in both past and future years. *target* is the number of individuals in
    # the target population. This is always shown for campaigns, and is now specified
    # at a national level. (See ‘coverage’ section above.) For routine, target is
    # shown as NA, which means you should assume the target population matches the
    # population shown in the demographic data downloads for the corresponding ages
    # (age_first and age_last). 
    
    # infer reached from coverage
    sia2 [, reached := as.numeric(target) * as.numeric(coverage)]
    
    
    # add additional columns to make file similar to original
    sia2 [, c("Geography", "Gavi73", "Gavi-supported", "Activity","Extent")] <- NA
    
    # reorder columns to make file similar to original
    sia2 <- sia2 [ , c("Geography", "country_code", "vaccine", "Gavi73", "Gavi-supported",
                       "Activity", "vaccine", "year", "a0", "a1", 
                       "Extent", "target", "reached", "coverage") ]
    
    sia2 [, vaccine := NULL]
    
    # write vaccine coverage data for SIAs and routine vaccination
    fwrite (x         = sia2, 
            file      = paste0 (output_dir, "coverage_sia_", scenario,".csv"),
            row.names = F)
    
    fwrite (x    = routine, 
            file = paste0 (output_dir,"coverage_routine_", scenario,".csv") )
    
  }
  
} # end of function -- generate_vaccine_coverage_routine_sia

