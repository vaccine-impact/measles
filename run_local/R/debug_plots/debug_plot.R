# debug plots (measles)

library (data.table)
library (ggplot2)

rm (list = ls ())

# read debug data file
dt <- fread ("ETH_debug_model.txt")

dt <- dt [V5 != "MCV1_eff"]

# convert age
dt [V4 <= 156, age := (V4 - 52) / 52]
dt [V4  > 156, age := (V4 - 154)]

# add year
dt [, year:= (V1 - 1) / 1000]

# compartment by percentage
dt [, N := sum (V6), by = year]
dt [, V6_pct := V6 / N * 100]

pdf ("debug_plot.pdf")

# absolute values
p <- ggplot (data = dt, 
        aes (x = year, y = V6, color = age)) + 
  geom_point () + 
  facet_wrap ( ~ V5) + 
  theme_bw()

print (p)

# > unique (dt$V5)
# [1] "M"        "S"        "I"        "R"        "VS"       "VI"       "VR"   
# "MCV1_eff" "V2S"      "V2I"      "V2R"     
# [12] "V3S"      "V3I"      "V3R" 

# compartmental states
states <- unique (dt$V5)

# loop through compartmental states
for (state in states) {
  
  state_dt <- dt [V5 == state]

  p <- ggplot (data = dt [V5 == state], 
          aes (x = year, y = V6, color = age)) + 
    geom_point () + 
    labs (title = state)
  
  print (p)
}

# percentage
p <- ggplot (data = dt, 
             aes (x = year, y = V6_pct, color = age)) + 
  geom_point () + 
  facet_wrap ( ~ V5) + 
  theme_bw()

print (p)

# > unique (dt$V5)
# [1] "M"        "S"        "I"        "R"        "VS"       "VI"       "VR"   
# "MCV1_eff" "V2S"      "V2I"      "V2R"     
# [12] "V3S"      "V3I"      "V3R" 

# compartmental states
states <- unique (dt$V5)

# loop through compartmental states
for (state in states) {
  
  state_dt <- dt [V5 == state]
  
  p <- ggplot (data = dt [V5 == state], 
               aes (x = year, y = V6_pct, color = age)) + 
    geom_point () + 
    labs (title = state)
  
  print (p)
  
  
}

dev.off ()




