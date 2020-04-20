library (data.table)
library (scales)


#################### Burden estimates -- central ####################

rm (list = ls ())

filename_bur_baseline   <- "baseline.csv" 
filename_bur_comparator <- "comparator.csv"

burden_baseline   <- fread (filename_bur_baseline)
burden_comparator <- fread (filename_bur_comparator)

# add birth cohort year
burden_baseline   [, birthcohort_year := year - age]
burden_comparator [, birthcohort_year := year - age]

# add scenario
burden_baseline     [, scenario := "baseline"]
burden_comparator   [, scenario := "comparator"]

# row bind to single table -- burden estimates of baseline and comparator scenarios
burden <- rbindlist (list (burden_baseline, 
                           burden_comparator), 
                     use.names = T)

# pdf ("central_estimate_baseline.pdf")
pdf ("central_estimate.pdf")

# plots of disease burden by birth cohort year (linear scale)
plotwhat = c("cohort_size", "cases", "deaths", "dalys")

for (i in 1:length(plotwhat)){
  
  toplot = plotwhat[i]
  
  p <- ggplot (burden [year >2019 & year <2031], aes(x = year)) +
    # geom_col (aes(y = get(toplot), col=age, fill = age)) +
    geom_col (aes(y = get(toplot), fill = age)) +
    scale_colour_gradientn(colours=rev(rainbow(5))) +
    theme_bw(base_size = 8) +
    labs (
      # x="birth cohort year",
      y=toplot) +
    facet_wrap (scenario ~ .)
  
  print (p)
  
  # for (j in 1:6) {
  #   print (p + facet_wrap_paginate(country_name ~ ., nrow=4, ncol=5, page=j, scales="free"))
  # }
}



# column bind to single table -- burden estimates of baseline and comparator scenarios
burden <-  burden_baseline [burden_comparator, 
                            on = .(disease = disease, 
                                   year    = year, 
                                   age     = age,
                                   country = country,
                                   country_name = country_name, 
                                   birthcohort_year = birthcohort_year)]

burden [, cohort_size_impact := i.cohort_size - cohort_size]  # check is 0

burden [, cases_impact  := i.cases  - cases]
burden [, deaths_impact := i.deaths - deaths]
burden [, dalys_impact  := i.dalys  - dalys]

# plots of disease burden by birth cohort year (linear scale)
plotwhat = c("cases_impact", "deaths_impact", "dalys_impact")

for (i in 1:length(plotwhat)){
  
  toplot = plotwhat[i]
  
  p <- ggplot (burden [year >2019 & year <2031, sum (get(toplot)), by = year], aes(x = year, y = V1)) +
    # geom_col (aes(y = get(toplot), col=age, fill = age)) +
    geom_bar (stat = "identity") +
    scale_colour_gradientn(colours=rev(rainbow(5))) +
    theme_bw(base_size = 8) +
    labs (
      # x="birth cohort year",
      y=toplot)
  
  print (p)
  
  # for (j in 1:6) {
  #   print (p + facet_wrap_paginate(country_name ~ ., nrow=4, ncol=5, page=j, scales="free"))
  # }
}




dev.off ()