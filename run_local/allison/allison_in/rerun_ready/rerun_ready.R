# get input burden estimate files of Allison into regular format

library (data.table)

rm (list = ls ())

# scenarios
scenarios <- c("campaign-only-bestcase",  # 1  SIAs only
               "mcv2-bestcase",           # 2  MCV1&2
               "campaign-bestcase",       # 3  MCV1&2 and SIAs
               "mcv1-bestcase",           # 4  MCV1 only
               "campaign-only-default",   # 5  SIAs only
               "mcv2-default",            # 6  MCV1&2
               "campaign-default",        # 7  MCV1&2 and SIAs
               "mcv1-default",            # 8  MCV1 only
               "no-vaccination",          # 9  no vaccination (set vaccination and using_sia to 0)
               "stop"                    # 10 MCV1&2 and SIAs
               )

burden_files <- paste0 ("burden_estimate_Measles-LSHTM-Jit-", scenarios, ".csv")

# loop through files
for (burden_file in burden_files) {
  
  # open file
  dt <- fread (burden_file)
  
  # change columns
  setnames (dt, "disease", "remain_lexp")
  dt [, disease := "Measles"]
  dt [, MCV1 := NULL]
  setcolorder (dt, "disease")

  # save file
  central_burden_file <- paste0 ("central_", burden_file)
  fwrite (dt, central_burden_file)
  
}
