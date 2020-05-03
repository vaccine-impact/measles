# lexp.R

# create remaining life expectancy file for each year across all age intervals
create_life_expectancy_remaining_full <- function () {
  
  # read data file for remaining life expectancy
  lexp_remain       <- fread ("input/demographicdata2019/201910gavi-4_dds-201910_2_life_ex_both.csv")
  
  lexp_remain_full <- copy (lexp_remain)
  
  # Since year is 5-year intervals, add data for all in-between years
  for (i in 1:4) {
    
    dt <- copy (lexp_remain) 
    dt [, year := year + i]
    
    # add table to full table of remaining life expectancy 
    lexp_remain_full <- rbindlist (list (lexp_remain_full, dt), 
                                   use.names = T)
  }
  
  # add data for 2100 for remaining life expectancy
  # copy values for year 2095 to year 2100 
  lexp_remain_2100 <- lexp_remain [year == 2095]
  lexp_remain_2100 [, year := 2100]
  
  lexp_remain_full <- rbindlist (list (lexp_remain_full, lexp_remain_2100), 
                                 use.names = T)
  
  # save file
  fwrite (lexp_remain_full, 
          file = "input/demographicdata2019/201910gavi-4_dds-201910_2_life_ex_both_full.csv")
  
} # end of function -- create_life_expectancy_remaining_full
