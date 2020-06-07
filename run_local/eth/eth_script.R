# set up 4 scenario files


library (data.table)

rm (list = ls ())

# read vaccine coverage file -- campaign-default (MCV1, MCV2, SIA)
wd <- getwd ()
setwd ("C:/Users/kajam/OneDrive/Documents/GitHub/measles/run_local/generate_coverage/VIMC-coverage-files")

# ------------------------------------------------------------------------------
# 1-time set-up to get ETH coverage
# coverage <- fread ("coverage_201910gavi-4_measles-campaign-default.csv")
# eth <- coverage [country_code == "ETH"]
# save base file and update with Matt's SIA numbers
# fwrite (eth, file = "eth_coverage.csv")
# ------------------------------------------------------------------------------

# coverage scenarios
# 1) Baseline (2019 RI coverage + SIA in 2020)
# 2) 50% reduction for 12 months in RI coverage
# 3) 50% reduction for 12 months in RI coverage + SIA postponed for 1 year
# 4) 50% reduction for 12 months in RI coverage + SIA cancelled


fname <- "coverage_201910gavi-4_measles-campaign-default"

# scenario 1) Baseline (2019 RI coverage + SIA in 2020)
eth <- fread ("eth_coverage_matt.csv")
fwrite (eth, file = paste0 (fname, "_ETH_s1.csv"))

# scenario 2) 50% reduction for 12 months in RI coverage
eth [year == 2020 & vaccine == "MCV1", coverage := coverage * 0.5]

fwrite (eth, file = paste0 (fname, "_ETH_s2.csv"))


# 3) 50% reduction for 12 months in RI coverage + SIA postponed for 1 year
eth [year == 2021 & activity_type == "campaign", `:=` (target    = eth [year == 2020 & activity_type == "campaign", target], 
                                                       coverage  = eth [year == 2020 & activity_type == "campaign", coverage], 
                                                       age_first = eth [year == 2020 & activity_type == "campaign", age_first],
                                                       age_last =  eth [year == 2020 & activity_type == "campaign", age_last]
                                                       )]

eth [year == 2020 & activity_type == "campaign", `:=` (target = 0, coverage = 0)]

fwrite (eth, file = paste0 (fname, "_ETH_s3.csv"))

# 4) 50% reduction for 12 months in RI coverage + SIA cancelled
eth [year == 2021 & activity_type == "campaign", `:=` (target = 0, coverage = 0)]

fwrite (eth, file = paste0 (fname, "_ETH_s4.csv"))


setwd (wd)