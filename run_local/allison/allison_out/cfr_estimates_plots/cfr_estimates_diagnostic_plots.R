# cfr_estimates_diagnostic_plots.R

# create: Single CFR file with cfr estimates for various VIMC scenarios and Portnoy scenarios 4 and 6
# create diagnostic cfr plots

# ------------------------------------------------------------------------------
# The output files of burden estimates in this folder are for Portnoy scenario 4.
# Scenario 4  (Portnou provided these values) : New CFRs estimated 2000 to 2030 
# including PSU- and LSHTM-specific incidence for best case and no vax scenarios
# 
# The streamlined CFR file is created for Portnoy scenario 6.
# Scenario 6 (To be calculated): New CFRs estimated 2000 to 2018 and CFRs held 
# constant 2018 to 2030 including PSU- and LSHTM-specific incidence for best case 
# and no vax scenarios.
# ------------------------------------------------------------------------------

# libraries
library (data.table)
library (ggplot2)
library (tictoc)

# clear workspace
rm (list = ls ())

# ------------------------------------------------------------------------------
# create: Single CFR file with cfr estimates for various VIMC scenarios and Portnoy scenarios 4 and 6
# create diagnostic cfr plots
cfr_estimates_diagnostic_plots <- function (cfr_estimates_file, 
                                            cfr_plot_file) {
  
  # diagnostic plot file -- cfrs at regional and country level for 
  # various VIMC scenarios and Portnoy scenarios 4 and 6
  pdf (cfr_plot_file)
  
  # data table with cfrs for various VIMC scenarios and Portnoy scenarios 4 and 6
  cfr_scenarios <- NULL
  
  # vimc scenarios -- touchstone 201910gavi
  vimc_scenarios <- c("campaign-only-bestcase",      # 1  SIAs only
                      "mcv2-bestcase",               # 2  MCV1&2 ** note the slight change in scenarios below
                      "campaign-bestcase",           # 3  MCV1&2 and SIAs
                      "mcv1-bestcase",               # 4  MCV1 only
                      "campaign-only-default",       # 5  SIAs only
                      "mcv2-default",                # 6  MCV1&2 ** note the slight change in scenarios below
                      "campaign-default",            # 7  MCV1&2 and SIAs
                      "mcv1-default",                # 8  MCV1 only
                      "no-vaccination",              # 9  no vaccination (set vaccination and using_sia to 0)
                      "stop"                         # 10 MCV1&2 and SIAs
  )
  
  
  # scenarios
  # change these when new scenarios are released:
  scenarios <- c("campaign-only-bestcase",           # 1  SIAs only
                 "mcv12-bestcase",                   # 2  MCV1&2
                 "campaign-bestcase",                # 3  MCV1&2 and SIAs
                 "mcv1-bestcase",                    # 4  MCV1 only
                 "campaign-only-default",            # 5  SIAs only
                 "mcv12-default",                    # 6  MCV1&2
                 "campaign-default",                 # 7  MCV1&2 and SIAs
                 "mcv1-default",                     # 8  MCV1 only
                 "no-vaccination",                   # 9  no vaccination (set vaccination and using_sia to 0)
                 "stop"                              # 10 MCV1&2 and SIAs
  )
  
  # scenario filenames
  scenarios <- paste0 ("cfr_Measles-LSHTM-Jit-", 
                       scenarios, ".csv")
  
  counter <- 0
  # loop through scenarios
  for (scenario in scenarios) {
    
    # loop counter
    counter <- counter + 1
    
    # read burden estimates along with cfr values
    burden <- fread (scenario)
    
    # ----------------------------------------------------------------------------  
    # copy CFR values of (Jordan to Palestine) and (Albania to Kosovo)
    
    countries           <- c("PSE", "XK")
    alternate_countries <- c("JOR", "ALB")
    
    for (i in 1:length (countries)) {
      
      burden_iso <- burden [country == countries [i]]
      burden_iso <- burden_iso [, ':=' (region = NULL, 
                                        U5MR = NULL, 
                                        cfr = NULL)
      ]
      
      # copy cfr values of alternative country
      burden_iso <- burden_iso [burden [country == alternate_countries [i]], 
                                on = .(year = year,
                                       age = age)
      ]
      
      # keep requisite columns
      burden_iso <- burden_iso [, c("disease", "year", "age", 
                                    "country", "country_name", 
                                    "cohort_size", "cases", "dalys", "deaths", 
                                    "region", "U5MR", "cfr")
      ] 
      
      # remove country with CFR - NA values
      burden <- burden [country != countries [i]]
      
      # add country with CFR values
      burden <- rbindlist (list (burden, burden_iso), 
                           use.names = TRUE)
    }
    # ----------------------------------------------------------------------------  
    
    # # ----------------------------------------------------------------------------  
    # # diagnostic plots
    # pdf (paste0 (wd, "/cfr_country_diagnostics_portnoy_scenario_4.pdf"))
    # for (country_code in unique (burden [, country])) {
    #   
    #   # diagnostic plot -- country-specific cfr plots
    #   p <- ggplot () +
    #     geom_line (data = burden [country == country_code & age < 5],
    #                aes (x = year, y = cfr), colour = "red") +
    #     geom_line (data = burden [country == country_code & age >= 5],
    #                aes (x = year, y = cfr), colour = "blue") +
    #     labs (title = country_code)
    #   
    #   print (p)
    # }
    # dev.off ()
    # 
    # pdf (paste0 (wd, "/cfr_region_diagnostics_scenario_4.pdf"))
    # for (region_name in unique (burden [, region])) {
    #   
    #   # diagnostic plot -- region-specific cfr plots
    #   p <- ggplot () +
    #     geom_line (data = burden [region == region_name & age < 5 & U5MR == ">=50"], 
    #                aes (x = year, y = cfr), colour = "red") +
    #     geom_line (data = burden [region == region_name & age < 5 & U5MR == "<50"], 
    #                aes (x = year, y = cfr), colour = "red") +
    #     geom_line (data = burden [region == region_name & age >= 5 & U5MR == ">=50"], 
    #                aes (x = year, y = cfr), colour = "blue") + 
    #     geom_line (data = burden [region == region_name & age >= 5 & U5MR == "<50"], 
    #                aes (x = year, y = cfr), colour = "blue") +
    #     labs (title = region_name)
    #   
    #   print (p)
    # }
    # dev.off ()
    # # ----------------------------------------------------------------------------  
    
    
    # ----------------------------------------------------------------------------
    # Note: CFR estimates for ages < 5 are the same for country (based of region and U5MR)
    # Similarly, CFR estimates for ages >= 5 are the same for country (based of region and U5MR)
    
    # CFR estimates for age < 5 and age >= 5 years
    cfr_less5 <- burden [age <  5]
    cfr_over5 <- burden [age >= 5]
    
    # get CFR summary estimates for country and year for ages <5 and >= 5 years
    cfr_less5 <- cfr_less5 [, .(under5 = mean(cfr)), 
                            by = .(country_name, country, region, U5MR, year)] 
    
    cfr_over5 <- cfr_over5 [, .(over5 = mean(cfr)), 
                            by = .(country_name, country, region, U5MR, year)] 
    
    # combine CFR estimates columns for age < 5 and age >= 5 years
    cfr <- cfr_less5 [cfr_over5, 
                      on = .(country_name = country_name, 
                             country      = country, 
                             region       = region,
                             U5MR         = U5MR,
                             year         = year)
    ]
    
    # ----------------------------------------------------------------------------
    
    
    # # ----------------------------------------------------------------------------  
    # # diagnostic plots
    # pdf (paste0 (wd, "/cfr_summary_country_diagnostics_portnoy_scenario_4.pdf"))
    # for (country_code in unique (cfr [, country])) {
    #   
    #   # diagnostic plot -- country-specific cfr plots
    #   p <- ggplot (data = cfr [country == country_code]) +
    #     geom_line (aes (x = year, y = under5), colour = "red") +
    #     geom_line (aes (x = year, y = over5), colour = "blue") +
    #     labs (title = country_code)
    #   
    #   print (p)
    # }
    # dev.off ()
    # # ----------------------------------------------------------------------------  
    
    
    # ----------------------------------------------------------------------------  
    # Portnoy CFRs (scenarios 4 and 6)
    
    # Scenario 4: New CFRs estimated 2000 to 2030 
    # including PSU- and LSHTM-specific incidence for best case and no vax scenarios
    
    # initialise cfr columns for portnoy scenario 6 (by copying values of scenario 4)
    cfr [, ':=' (under5_s6 = under5, 
                 over5_s6  = over5)]
    
    # Scenario 6: New CFRs estimated 2000 to 2018 and CFRs held 
    # constant 2018 to 2030 including PSU- and LSHTM-specific incidence for best case 
    # and no vax scenarios.
    for (iso_code in unique (cfr [, country])) {
      
      cfr [country == iso_code & year > 2018, 
           ':=' (under5_s6 = cfr [country == iso_code & year == 2018, under5_s6], 
                 over5_s6  = cfr [country == iso_code & year == 2018, over5_s6] )
      ]
    }
    
    
    # ---------------------------------------------------------------------------- 
    # add vimc scenario name to cfr column name (and portnoy scenario 4)
    setnames (x = cfr,
              old = c ("under5",
                       "over5"
              ), 
              new = c (paste0 ("cfr_under5_", vimc_scenarios [counter], "_s4"),
                       paste0 ("cfr_over5_" , vimc_scenarios [counter], "_s4") 
              )
    )
    
    # add vimc scenario name to cfr column name (and portnoy scenario 6)
    setnames (x = cfr,
              old = c ("under5_s6",
                       "over5_s6"
              ), 
              new = c (paste0 ("cfr_under5_", vimc_scenarios [counter], "_s6"),
                       paste0 ("cfr_over5_" , vimc_scenarios [counter], "_s6") 
              )
    )
    # ---------------------------------------------------------------------------- 
    
    
    # add cfrs of current scenarios to full cfr table of 
    # various VIMC scenarios and Portnoy scenarios 4 and 6
    if (is.null (cfr_scenarios)) {
      cfr_scenarios <- cfr
      
    } else {
      
      cfr_scenarios <- cfr_scenarios [cfr, 
                                      on = .(country_name = country_name, 
                                             country      = country, 
                                             region       = region,
                                             U5MR         = U5MR,
                                             year         = year)
      ]
    }
    
    
    # ---------------------------------------------------------------------------- 
    # diagnostic plots -- regional level
    for (region_name in unique (cfr [, region])) {
      
      # diagnostic plot -- region-specific cfr plots
      p <- ggplot () +
        
        geom_line (data = cfr [region == region_name & U5MR == ">=50"], 
                   aes (x = year, 
                        y = eval (as.name (paste0 ( "cfr_over5_", vimc_scenarios [counter], "_s4"))),     
                        colour = paste0 ( "U5MR>=50_over5_", vimc_scenarios [counter], "_s4"))) +
        geom_line (data = cfr [region == region_name & U5MR == "<50"], 
                   aes (x = year, 
                        y = eval (as.name (paste0 ( "cfr_over5_", vimc_scenarios [counter], "_s4"))),     
                        colour = paste0 ( "U5MR<50_over5_", vimc_scenarios [counter], "_s4"))) +
        geom_line (data = cfr [region == region_name & U5MR == ">=50"], 
                   aes (x = year, 
                        y = eval (as.name (paste0 ( "cfr_under5_", vimc_scenarios [counter], "_s4"))),     
                        colour = paste0 ( "U5MR>=50_under5_", vimc_scenarios [counter], "_s4"))) +
        geom_line (data = cfr [region == region_name & U5MR == "<50"],
                   aes (x = year, 
                        y = eval (as.name (paste0 ( "cfr_under5_", vimc_scenarios [counter], "_s4"))),     
                        colour = paste0 ( "U5MR<50_under5_", vimc_scenarios [counter], "_s4"))) +
        
        geom_line (data = cfr [region == region_name & U5MR == ">=50"], 
                   aes (x = year, 
                        y = eval (as.name (paste0 ( "cfr_over5_", vimc_scenarios [counter], "_s6"))),     
                        colour = paste0 ( "U5MR>=50_over5_", vimc_scenarios [counter], "_s6"))) +
        geom_line (data = cfr [region == region_name & U5MR == "<50"], 
                   aes (x = year, 
                        y = eval (as.name (paste0 ( "cfr_over5_", vimc_scenarios [counter], "_s6"))),     
                        colour = paste0 ( "U5MR<50_over5_", vimc_scenarios [counter], "_s6"))) +
        geom_line (data = cfr [region == region_name & U5MR == ">=50"], 
                   aes (x = year, 
                        y = eval (as.name (paste0 ( "cfr_under5_", vimc_scenarios [counter], "_s6"))),     
                        colour = paste0 ( "U5MR>=50_under5_", vimc_scenarios [counter], "_s6"))) +
        geom_line (data = cfr [region == region_name & U5MR == "<50"],
                   aes (x = year, 
                        y = eval (as.name (paste0 ( "cfr_under5_", vimc_scenarios [counter], "_s6"))),     
                        colour = paste0 ( "U5MR<50_under5_", vimc_scenarios [counter], "_s6"))) +
        
        labs (title    = region_name, 
              subtitle = vimc_scenarios [counter],
              y        = "CFR", 
              colour   = "measles CFR")
      
      print (p)
    }
    
    
    # diagnostic plots -- country level
    for (country_code in unique (cfr [, country])) {
      
      # diagnostic plot -- country-specific cfr plots
      p <- ggplot (data = cfr [country == country_code]) +
        geom_line (aes (x = year, 
                        y = eval (as.name (paste0 ( "cfr_over5_", vimc_scenarios [counter], "_s4"))),     
                        colour = paste0 ( "over5_", vimc_scenarios [counter], "_s4"))) +
        geom_line (aes (x = year, 
                        y = eval (as.name (paste0 ( "cfr_under5_", vimc_scenarios [counter], "_s4"))),     
                        colour = paste0 ( "under5_", vimc_scenarios [counter], "_s4"))) +
        geom_line (aes (x = year, 
                        y = eval (as.name (paste0 ( "cfr_over5_", vimc_scenarios [counter], "_s6"))),     
                        colour = paste0 ( "over5_", vimc_scenarios [counter], "_s6"))) +
        geom_line (aes (x = year, 
                        y = eval (as.name (paste0 ( "cfr_under5_", vimc_scenarios [counter], "_s6"))),     
                        colour = paste0 ( "under5_", vimc_scenarios [counter], "_s6"))) +
        labs (title    = country_code, 
              subtitle = vimc_scenarios [counter], 
              y        = "CFR", 
              colour   = "measles CFR")
      
      print (p)
    }
    # ---------------------------------------------------------------------------- 
    
  }
  
  # rename columns -- (country to country_code) and (region to GBD_region)
  setnames (x = cfr_scenarios, 
            old = c ("country",      "region"), 
            new = c ("country_code", "GBD_region") 
            )
  
  # save single CFR file with cfr estimates for various VIMC scenarios 
  # and Portnoy scenarios 4 and 6
  fwrite (x    = cfr_scenarios, 
          file = cfr_estimates_file)
  
  # close plot file
  dev.off ()
  
  return ()
  
} # end of function -- cfr_estimates_diagnostic_plots
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# main program -- start
# ------------------------------------------------------------------------------
tic ()
print (Sys.time ())

# working directory
wd <- getwd ()
setwd ("../")

# create: Single CFR file with cfr estimates for various VIMC scenarios and Portnoy scenarios 4 and 6
# create diagnostic cfr plots
cfr_estimates_diagnostic_plots (cfr_estimates_file = paste0 (wd, "/cfr_scenarios.csv"),
                                cfr_plot_file      = paste0 (wd, "/cfr_region_country_s4_s6.pdf")
)

# working directory
setwd (wd)

toc ()
print (Sys.time ())
# ------------------------------------------------------------------------------
# main program -- end
# ------------------------------------------------------------------------------
