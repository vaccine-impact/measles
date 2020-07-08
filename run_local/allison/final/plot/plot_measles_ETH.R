# plot script - Ethiopia measles scenario

library (data.table)
library (ggplot2)

rm (list = ls ())

# plot vaccine coverage and disease burden estimates of comparative scenarios
plot_coverage_burden <- function () {
  
  # option
  for (opt in c("vimc", "matt") ) {

    # ----------------------------------------------------------------------------
    # Vaccine coverage
    all_vac_cov <- NULL
    
    # vaccine coverage directory
    vac_cov_dir <- paste0 ("../", opt)
    
    # scenarios
    for (s in 1:4) {
      
      # vaccine coverage filename
      filename_vac_cov <- paste0 (vac_cov_dir, "/coverage_201910gavi-4_measles-campaign-default_ETH_s", s, "_", opt, ".csv")
      
      # read vaccine coverage
      vac_cov <- fread (filename_vac_cov)
      
      # add scenario number # note: "scenario" is already a column in vaccine coverage file
      vac_cov [, scenario_s := s]
      
      # create data table of vaccine coverage estimates for all scenarios
      if (is.null(all_vac_cov)) {
        all_vac_cov <- vac_cov
        
      } else {
        all_vac_cov <- rbindlist (list (all_vac_cov, vac_cov), use.names = TRUE)
      } 
    }
    # ----------------------------------------------------------------------------
    
    
    # ----------------------------------------------------------------------------
    # Burden estimates -- central
    all_burden <- NULL
    
    # burden directory
    burden_dir <- paste0 ("../", opt)
    
    # scenarios
    for (s in 1:4) {
      
      # burden filename
      filename_bur <- paste0 (burden_dir, "/central_burden_estimate_Measles-LSHTM-Jit-campaign-default_ETH_s", s, "_", opt, ".csv") 
      
      # read burden estimates
      burden <- fread (filename_bur)
      
      # add birth cohort year
      burden [, birthcohort_year := year - age]
      
      # add scenario number
      burden [, scenario := s]
      
      # create data table of burden estimates for all scenarios
      if (is.null(all_burden)) {
        all_burden <- burden
        
      } else {
        all_burden <- rbindlist (list (all_burden, burden), use.names = TRUE)
      } 
    } 
    # ----------------------------------------------------------------------------
    
    
    # ----------------------------------------------------------------------------
    # plots for comparative analysis across scenarios
    
    # plot results file
    pdf (paste0 ("comparative_scenarios_ETH_", opt, ".pdf"))
    
    # plot vaccine coverage (routine)
    p <- ggplot (all_vac_cov [year >= 2015 & year <= 2030 & vaccine == "MCV1" & activity_type == "routine"], 
                 aes (year, y = coverage, group = scenario_s, color = factor (scenario_s))) + 
      geom_line () + 
      ylab ("MCV1 routine coverage") +
      labs (color='scenario') +  
      theme_bw ()

    print (p)
    
    # plot vaccine coverage (campaign)
    p <- ggplot (all_vac_cov [year >= 2015 & year <= 2030 & vaccine == "Measles" & activity_type == "campaign"], 
                 aes (year, y = coverage, group = scenario_s, color = factor (scenario_s))) + 
      geom_line () + 
      ylab ("Measles campaign coverage") +
      labs (color='scenario') +  
      theme_bw ()
    
    print (p)
    
    # plot burden -- cases, deaths, dalys
    plotwhat = c("cases", "deaths", "dalys")
    
    for (i in 1:length(plotwhat)){
      
      toplot = plotwhat[i]
      
      p <- ggplot(all_burden [year >= 2015 & year <= 2030], 
                  aes(x = year, y = get(toplot), group = scenario, color = factor (scenario))) +
        stat_summary (fun = sum, geom = "line") +
        ylab (toplot) +
        labs (color='scenario') + 
        theme_bw ()
      
      print (p)
      
    }
    
    dev.off ()
    
  } # end of -- for (opt in c("vimc", "matt") )
  
} # end of function -- plot_coverage_burden


# plot vaccine coverage and disease burden estimates of comparative scenarios
plot_coverage_burden ()



# # plots of disease burden by birth cohort year (linear scale)
# plotwhat = c("cohort_size", "cases", "deaths", "dalys")
# 
# for (i in 1:length(plotwhat)){
#   
#   toplot = plotwhat[i]
#   
#   p <- ggplot (burden, aes(x = birthcohort_year)) +
#     geom_col (aes(y = get(toplot), col=age)) +
#     scale_colour_gradientn(colours=rev(rainbow(5))) +
#     theme_bw(base_size = 8) +
#     labs (
#       x="birth cohort year",
#       y=toplot)
#   
#   print (p)
#   
#   # for (j in 1:6) {
#   #   print (p + facet_wrap_paginate(country_name ~ ., nrow=4, ncol=5, page=j, scales="free"))
#   # }
# }
# 
# dev.off ()
# 
# }



