tic ()

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
  countries                         = c("ETH"),  # debug -- c("all"), 
  
  cluster_cores                     = 2,  # number of cores
  psa                               = 0   # psa runs; 0 for single run
)


# ------------------------------------------------------------------------------
# diagnostic_plots
#
# diagnostic plots of vaccine coverage and burden estimates (cases, deaths, dalys)
# ------------------------------------------------------------------------------
diagnostic_plots <- function (vaccine_coverage_folder,
                              coverage_prefix,
                              touchstone,
                              antigen,
                              scenarios,
                              burden_estimate_folder,
                              plot_folder,
                              group_name,
                              countries,
                              cfr_options, 
                              psa, 
                              start_year = -1,
                              end_year   = -1) {
  
  # burden estimate type -- central or stochastic
  if (psa > 0) {
    burden_estimate_type <- "stochastic_burden_estimate_"
  } else {
    burden_estimate_type <- "central_burden_estimate_"
  }
  
  # cfr (case fatality rate) options (Wolfson and/or Portnoy)
  for (cfr_option in cfr_options) {
    
    # diagnostic plots filename
    pdf (paste0 (plot_folder,
                 "diagnostic_plot_",
                 antigen,
                 cfr_option, 
                 ".pdf"))
    
    # burden estimates of all scenarios
    all_burden <- NULL
    
    # scenarios
    for (scenario_name in scenarios) {
      
      # vaccine coverage file
      vaccine_coverage_file <- paste0 (vaccine_coverage_folder, 
                                       coverage_prefix, 
                                       touchstone, 
                                       antigen,
                                       scenario_name, 
                                       ".csv")
      
      # burden estimate filename
      burden_estimate_file <- paste0 (burden_estimate_type,
                                      antigen, 
                                      group_name, 
                                      scenario_name, 
                                      ".csv")
      
      
      # append/suffix cfr_option to the end of filename
      updated_burden_estimate_file <- 
        str_replace (string      = burden_estimate_file, 
                     pattern     = ".csv", 
                     replacement = paste0 ("_", cfr_option, ".csv"))
      
      # burden file with folder
      burden_file = paste0 (burden_estimate_folder,
                            cfr_option, "/", 
                            updated_burden_estimate_file)
      
      # read data -- vaccine coverage and burden estimates
      vaccine_coverage <- fread (vaccine_coverage_file)
      burden_estimate  <- fread (burden_file)
      
      # set start and end year if not set
      if (start_year == -1) {
        start_year <- min (burden_estimate [, year])
      }
      if (end_year == -1) {
        end_year <- max (burden_estimate [, year])
      }
      
      # extract vaccine coverage and burden estimates between start and end years
      vaccine_coverage <- vaccine_coverage [year >= start_year & year <= end_year]
      burden_estimate  <- burden_estimate  [year >= start_year & year <= end_year]
      
      # set vaccine (type) to SIA for campaign immunisation
      vaccine_coverage [activity_type == "campaign", vaccine := "SIA"]
      
      # add scenario name to burden estimate data table
      burden_estimate [, scenario := scenario_name]
      
      # combine burden estimates of all scenarios
      if (is.null (all_burden)) {
        all_burden <- burden_estimate
      } else {
        all_burden <- rbindlist (list (all_burden, 
                                       burden_estimate), 
                                 use.names = TRUE)
      }
      
      # if countries are specified to all, then set countries to all countries in coverage file
      if (countries == "all") {
        countries	<- as.character (unique (burden_estimate [, country] ) )  
      }
      
      # iso3 country codes
      country_iso3_codes <- countries
      
      # plot for each country
      for (country_iso3_code in country_iso3_codes) {
        
        # plot vaccine coverage
        coverage_plot <- ggplot (data = vaccine_coverage [country_code == country_iso3_code], 
                                 aes (x = year,
                                      y = coverage * 100, 
                                      color = factor (vaccine))) +  
          scale_x_continuous (breaks = pretty_breaks ()) + 
          geom_line () + 
          labs (title = countrycode (sourcevar   = country_iso3_code, 
                                     origin      = "iso3c", 
                                     destination = "country.name"),
                x = "Year", 
                y = "Vaccine coverage (%)",
                colour = "vaccine") + 
          theme_bw ()
        
        
        # plot burden -- cases, deaths, dalys
        plotwhat       <- c("cases", "deaths", "dalys")
        plotwhat_label <- c("Cases", "Deaths", "DALYs")
        
        burden_plot_list <- lapply (1:length(plotwhat), function (i) {
          
          toplot = plotwhat[i]
          
          p <- ggplot(burden_estimate [country == country_iso3_code], 
                      aes(x = year, y = get(toplot))) +
            scale_x_continuous (breaks = pretty_breaks ()) + 
            stat_summary (fun = sum, geom = "line") +
            ylab (toplot) + 
            labs (title = countrycode (sourcevar   = country_iso3_code, 
                                       origin      = "iso3c", 
                                       destination = "country.name"),
                  x = "Year", 
                  y = plotwhat_label [i]) + 
            theme_bw ()
        })
        
        # list of plots
        plot_list <- list (coverage_plot, 
                           burden_plot_list [[1]], 
                           burden_plot_list [[2]], 
                           burden_plot_list [[3]])
        
        # arrange plots in a single page
        plots <- ggarrange (plotlist = plot_list, 
                            ncol = 2, 
                            nrow = 2)
        
        # print plots
        print (annotate_figure (plots,
                                top = text_grob (scenario_name,
                                                 color = "black",
                                                 size = 9)))
        
      } # end of loop -- for (country_iso3_code in country_iso3_codes)
      
    } # end of loop -- for (scenario_name in scenarios)
    
    # add comparative plot across all scenarios for each country
    for (country_iso3_code in country_iso3_codes) {
      
      # plot burden -- cases, deaths, dalys
      plotwhat       <- c("cases", "deaths", "dalys")
      plotwhat_label <- c("Cases", "Deaths", "DALYs")
      
      for (i in 1:length (plotwhat)) {
        
        toplot = plotwhat [i]
        
        p <- ggplot(all_burden [country == country_iso3_code], 
                    aes(x     = year, 
                        y     = get (toplot), 
                        group = scenario, 
                        color = factor (scenario))) +
          scale_x_continuous (breaks = pretty_breaks ()) + 
          stat_summary (fun = sum, geom = "line") +
          ylab (toplot) + 
          labs (title = countrycode (sourcevar   = country_iso3_code, 
                                     origin      = "iso3c", 
                                     destination = "country.name"),
                x = "Year", 
                y = plotwhat_label [i], 
                colour = "Scenario") + 
          theme_bw ()
        
        print (p)
      }
      
      
    } # end of loop -- for (country_iso3_code in country_iso3_codes)
    
    dev.off ()
    
  } # end of loop -- for (cfr_option in cfr_options)

} # end of function -- diagnostic_plots
# ------------------------------------------------------------------------------


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

# scenarios <- c("counterfactual-bau-scenario1")

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
  # cfr_options                = c("Wolfson"),
  psa                        = var$psa,
  start_year                 = 2015,
  end_year                   = 2030
)
# ------------------------------------------------------------------------------

toc ()


