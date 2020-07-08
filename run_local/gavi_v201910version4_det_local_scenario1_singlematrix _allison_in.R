# gavi_v201910version4_det_local_scenario1_singlematrix _allison_in.R

# load libraries
library (tictoc)
library (data.table)
library (doParallel)

# clear workspace 
rm (list = ls ())

tic ()  # start timer
print (Sys.time ())

# Title:		Run mathematical model for measles to update estimates for Gavi
# Author:		Ed Jones, London School of Hygiene & Tropical Medicine. Edited by Kevin van Zandvoort, London School of Hygiene & Tropical Medicine, Petra Klepac, LSHTM
# Date:			20/11/2019
# Description:	Main file to run VIMC 2019 

# Changelog:	This file has been heavily changed compared to old file by Ed Jones:
# 				1.	There is only one file in which different settings can be changed to run different scenario's, rather than 8 different folders for different scenarios;
# 				2.	Variables that are created are no longer manually changed. Before, it was necessary to adapt multiple values in order to allow to run the script for
#					different number of years and/or countries. This now happens dynamically based on the coverage data.
# 				3. Age stratified output for Montagu system is now automatically generated.
#				4. No longer does a new vaccine.exe needs to be compiled when the number of years that is modelled changes. This is now passed directly to the model in shell.

# changes 2019: 
#       a) add scenario folder pathname for clarity and for easier running on the cluster (fortran code changed too)
#       b) change contact matrix to fix augumented mixing in <3 yo and to reflect country specific mixing (rescale polymod by local population structure to keep contacts reciprocal)
#          - uses newly rescaled contact matrix to expand polymod to yearly ages up to 100
#       c) time-varying CFRs (postprocessing)
#       d) MCV1 and MCV2 dependency - give MCV2 to those who received MCV1 for MCV2 coverage <= MCV1, if MCV2 coverage > MCV1 distribute the remainder randomly
#       e) MCV1 and SIA dependency - use linear regression of data from Portnoy et al 2018 to inform the proportion of zero-dose children reached by SIA campaign
#       f) due to issues with age_range_verbatim field in the scenario (coverage) input files, these runs use age_first and age_last

# ------------------------------------------------------------------------------
# USER-DEFINED OPTIONS	
# I  Working directory for Measles model in General (full-path including trailing slash)
# ------------------------------------------------------------------------------

run.local <- TRUE

version <- "201910version4"

# set working directory
if (run.local) {
  wd <- paste0 (getwd (), "/")
  # wd <- paste0("~/Dropbox/VIMC/measles/201910version4/run_local/")
  # wd.noslash <- paste0("~/Dropbox/VIMC/measles/201910version4/run_local")
} else {
  wd <- paste0("/Users/petra/OneDrive - London School of Hygiene and Tropical Medicine/VIMC/VIMC-2019runs/version", 
               version, "/")
  wd.noslash <- paste0("/Users/petra/OneDrive - London School of Hygiene and Tropical Medicine/VIMC/VIMC-2019runs/version", 
                       version)
}

## if setting directories for the first time, this creates figures and outcome folders
#  project.name <- "VIMC-2019runs"
# 
# # create folders to save figures and outputfiles
# fig.path <- dir.create(file.path(wd.noslash, project.name, "figures"))
# output.path <- dir.create(file.path(wd.noslash, project.name, "outcome"))
# 


# scenario index to run
index <- 2
# for (index in 1:10) {  # debug #
for (index in 9:10) {  # debug #
  
  
  #  change these when new scenarios are released:
  scenario.name <- c("campaign-only-bestcase",           # 1  SIAs only
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
  
  scenario <- scenario.name [index]
  
  # index variable for scenario to save it to a correct folder (scenarioXX, XX = 01,02,..... ,10)
  counter       <- c(paste0 ("0", c(1:9)), "10")
  save.scenario <- paste0 ("scenario", counter[index ])  # debug #
  # save.scenario <- paste0 ("scenario", as.character(index))
  
  # set SIAs and vaccination parameters for each scenario to minimize errors for running
  set.sia         <- c(1, 0, 1, 0, 1, 0, 1, 0, 0, 1)
  set.vaccination <- c(0, 2, 2, 1, 0, 2, 2, 1, 0, 2)
  
  # ------------------------------------------------------------------------------
  # II		Options
  # ------------------------------------------------------------------------------
  vaccination	<- set.vaccination [index]  # Whether children are vaccinated. 0: No vaccination; 1: Only MCV1; 2: MCV1 and MCV2
  using_sia	  <- set.sia [index]     	    # Whether supplementary immunization campaigns are used. 0: no SIA; 1: with SIA
  ages 		    <- c (0:100)	              # Numeric vector with age-strata that are reported (Model ALWAYS models 101 age-groups; 0-2yo in weekly age-strata, and 3-100yo in annual age-strata)
  psa 		    <- 0		                    # Number of PSAs to run (0 to use deterministic).
  
  
  # ------------------------------------------------------------------------------
  # III	 Input-data
  # ------------------------------------------------------------------------------
  if (index == 9) {
    
    data_coverage_routine 	<- paste0("scenarios2019/coverage_routine_", scenario.name[index-1], ".csv")
    # file with SIA data (csv; in ./input/)
    data_coverage_sia 		<- paste0("scenarios2019/coverage_sia_", scenario.name[index-1], ".csv") # "SIA_Gavi_1.csv" #"coverage_sia_imputed.csv"
    
  } else {
    
    data_coverage_routine 	<- paste0("scenarios2019/coverage_routine_", scenario, ".csv")
    # file with SIA data (csv; in ./input/)
    data_coverage_sia 		<- paste0("scenarios2019/coverage_sia_", scenario, ".csv")  
    
  }
  
  if (index %in% c(1,5)) {
    data_coverage_routine 	<- paste0("scenarios2019/coverage_routine_", scenario.name[2], ".csv")
  }
  
  # file with timeliness data (csv; in ./input/)
  data_timeliness <- "timeliness_new.csv"
  
  # file with R0 (csv; in ./input/)
  data_r0 		<- "r0.csv" #"r0est.csv"
  
  # file with life-expectancies (csv; in ./input/)
  data_life_exp 	<- "demographicdata2019/201910gavi-4_dds-201910_lx0_both.csv"
  #data_life_exp 	<- "new/lx_old.csv"
  
  # template Gavi (used to check if all data is present)
  data_gavi_template <- "central-burden-template_2019104.csv" 
  
  # file with CFRs (csv; in ./input/)
  data_cfr 		<- "cfrs_new_noage.csv"
  
  # file with population sizes (csv; in ./input/)
  data_pop 		<- "demographicdata2019/201910gavi-4_dds-201910_2_int_pop_both.csv"
  #data_pop 		<- "new/pop_old.csv"
  
  # file with contact_data
  data_contact	<- "contact/uk_polymod_physical_101.csv"  # use newly rescaled polymod to yearly age bands to 100
  
  # data with PSA variables
  # should be the same for each scenario. Will be CREATED if doesnt exist
  data_psa 		<- "psa_variables.csv"
  
  # expected remaining years of life
  # data_life_exp_remain <- "demographicdata2019/201910gavi-4_dds-201910_2_life_ex_both.csv"
  data_life_exp_remain <- "demographicdata2019/201910gavi-4_dds-201910_2_life_ex_both_full.csv"
  
  
  # ------------------------------------------------------------------------------
  # IV		Advanced options (only change these if you know what you're doing!)
  # ------------------------------------------------------------------------------
  
  sia.method	<- 1				             # 1=variable take, 2=variable degree
  dinf		    <- 14				             # duration of infection (days)
  amplitude 	<- 0.05			             # amplitude for seasonality
  take 		    <- c (0.85, 0.95, 0.98)  # vaccine efficacy for take1 (take dose 1, before age 1), take2 (dose 1, after age 1) & take3 (dose 2). Note that dose2 only has an effect if vaccine==2.
  degree 		  <- c (0.85, 0.95, 0.98)  # vaccine efficacy for degree1 (degree dose 1, before age 1), deree2 (dose 1, after age 1) & degree3 (dose 2). Note that dose2 only has an effect if vaccine==2.			
  tstep			  <- 1000				           # Number of time steps in a year
  
  # filename of compiled fortran-model (in ./model/compiled/)
  # measles_model <- "vaccine2019_sia_singlematrix" # change this to reflect the right version of the fortran code
  measles_model <- "vaccine2019_sia_singlematrix.exe" # change this to reflect the right version of the fortran code
  
  # number of clusters to use
  # if larger than 1, country-specific model runs are distributed over specified number of clusters
  # note that model uses a lot of memory, so might not want to max out all clusters
  use_cluster  <- 3   # debug #
  remove_files <- T
  
  # may want to process results after generating all data. Note OUTPUT files are not removed if remove_files == TRUE and process_results == FALSE
  process_results <- F
  run_model       <- TRUE
  # folder will be created if not given - should usually be commented out, except when run model is FALSE
  
  
  # ------------------------------------------------------------------------------
  # V		Debug
  # ------------------------------------------------------------------------------
  
  debug_country		  <- "*"			#ISO3 codes of country to debug, * to debug all countries
  debug_spinup		  <- FALSE		#TRUE/FALSE: If true, generate data for spin-up period of model
  debug_model			  <- FALSE		#TRUE/FALSE: If true: generate data for period after spin-up
  debug_compartments<- 0			  #TRUE/FALSE: If true: output size of each compartment. If false: output number of cases. If 2: debug vaccinated
  debug_age         <- 0        #0-2. If 0: output all in annual age-strata. If 1: output age 0-2 in weekly age-strata, 3-100 in annual age-strata. If 2: sum all age-strata.
  debug_timepoints	<- 0			  #0-2. If 0: output per year. If 1: output per timepoint and report first 25% of timepoints. If 2: output per timepoint and report all timepoints.
  debug_relative		<- FALSE		#If true: output proportion of new cases. If false, output absolute number of new cases.
  
  # START OF MODEL  
  
  # setup
  # wd <- paste0(wd, project.name, "/")
  setwd (wd)			
  source (paste0("./gavi_setup_v", version, "_singlematrix.R"))
  
  # log
  writelog("gavi_log",paste0("Main; gavi.r started"))
  if(using_openmpi){
    writelog("gavi_log",paste0("Main; OpenMPI enabled"))
  } else {
    writelog("gavi_log",paste0("Main; OpenMPI disabled"))
  }
  
  # read_data
  coverage_routine	<- fread (paste0 ("input/", data_coverage_routine))
  coverage_sia		  <- fread (paste0 ("input/", data_coverage_sia))
  timeliness			  <- fread (paste0 ("input/", data_timeliness))
  cfr					      <- fread (paste0 ("input/", data_cfr))
  rnought			    	<- fread (paste0 ("input/", data_r0))
  population		  	<- fread (paste0 ("input/", data_pop))
  lexp 		      		<- fread (paste0 ("input/", data_life_exp))
  template		    	<- fread (paste0 ("input/", data_gavi_template))
  contact			    	<- fread (paste0 ("input/", data_contact))
  lexp_remain       <- fread (paste0 ("input/", data_life_exp_remain)) 
  
  # coverage_sia has multiple entries per year for some countries, take only the first entry
  coverage_sia <- coverage_sia[, .SD[1], by = c("country_code", "year")]
  
  if (psa > 0) {
    if (file.exists(paste0("input/",data_psa))){
      # read csv if file already exists
      psa_var <- fread(paste0("input/",data_psa))
      
      # check if psa_var corresponds with psa
      if(nrow(psa_var) != psa){
        stop(paste0("Number of runs in PSA file is not the same as those specified! Variables used in the same run in each scenario should be similar. Delete or rename './input/",data_psa," if a new file needs to be created, or change the number of PSAs to ",nrow(psa_var)," if the same file should be used."))
      }
    } else {
      
      # create csv if file does not exist
      psa_var <- data.table(
        run_id=c(1:psa),
        take1_input=rep(NA,psa),
        take2_input=rep(NA,psa),
        take3_input=rep(NA,psa),
        mortality_input=rep(NA,psa)
      )
      
      for(t in 1:3){
        # use 5% difference in take
        tk <- runif(n=psa,min=(take[t]*0.95),max=(take[t]*1.05))
        # vaccine efficacy is bounded between 0 and 1
        tk [which(tk < 0)] <- 0
        tk [which(tk > 1)] <- 1
        psa_var[,paste0("take",t,"_input")] <- tk
      }
      
      # use 25% difference in mortality, will be multiplied by CFR in each country
      psa_var [,"mortality_input"] <- runif(n=psa,min=0.75,max=1.25)
      fwrite (psa_var,paste0("input/",data_psa))
    }
  }
  
  # process data
  countries	<- as.character(unique(coverage_routine[,country_code]))  
  # countries	<- as.character(unique(coverage_routine[,country_code]))[1]  # debug #
  # countries	<- c ("ETH", "BGD", "SSD", "NGA")
  # countries	<- c ("ETH")
  # countries	<- as.character(unique(coverage_routine[,country_code]))
  
  # change to run only for a set of specified missing countries
  # countries <- missing_countries
  
  # years 		<- sort(unique(as.numeric(coverage_routine[,year]))) this no longer works as some imput files don't have any years
  years <- as.numeric(c(1980:2100))
  
  for(c in countries){
    if(psa > 1){
      for(p in 1:psa){
        write(
          paste0(
            c,
            " ",
            paste0(c(rep(0,(nchar(psa)-nchar(p))),p),collapse=""),
            " 0"
          ),
          file="gavi_progress",
          append=TRUE
        )
      }
    } else {
      write(
        paste0(
          c,
          " ",
          "1 0"
        ),
        file="gavi_progress",
        append=TRUE
      )
    }
  }
  
  contact <- contact[,-"contact.age"]
  contact <- as.matrix(contact)
  
  
  r0_basic <- Re (eigen (contact*dinf, only.values=T)$values[1])
  gamma    <- 1/(dinf*tstep/365)
  
  # Run model
  writelog ("gavi_log", paste0 ("Main; Foldername: ", foldername))
  
  if(using_openmpi | use_cluster > 1){
    if(!using_openmpi){
      if( use_cluster != detectCores() ){
        warning(paste0(detectCores(), " cores detected but ", use_cluster, " specified."))
      }
      cl <- makeCluster(use_cluster)
      registerDoParallel(cl)
    } else {
      print(paste0("Clustersize: ",clusterSize(cl)))
    }
  }
  
  # foreach will run countries and PSA runs in parralel if a parallel backend 
  # is registered, and sequentially otherwise
  combine <- foreach (
    ii = 1:length(countries),
    .packages = c("data.table"),
    .errorhandling="pass"
  ) %:% foreach (
    r = 1:runs,
    .errorhandling = "pass"
  ) %dopar% {
    
    out_run <- runCountry (ii, 
                           countries, 
                           years,
                           vaccination, 
                           using_sia, 
                           dinf, 
                           gamma, 
                           r0_basic, 
                           amplitude, 
                           take, 
                           degree, 
                           sia_method,
                           coverage_routine, 
                           coverage_sia, 
                           timeliness, 
                           contact, 
                           rnought, 
                           population, 
                           lexp, 
                           cfr,
                           tstep,save.scenario, 
                           foldername,  
                           measles_model, 
                           debug_model, 
                           debug_spinup, 
                           debug_age, 
                           debug_compartments, 
                           debug_country, 
                           debug_relative, 
                           debug_timepoints,
                           r,
                           runs, 
                           psa, 
                           psa_var,
                           process_results, 
                           run_model, 
                           remove_files
                           )
    return(out_run)
  }
  
  if(using_openmpi | use_cluster > 1){
    if(!using_openmpi){
      stopCluster(cl)
    } else {
      closeCluster(cl)
    }
  }
  
  # check for errors
  errorcount <- 0
  
  for (i in 1:length(combine)){
    if("error" %in% class(combine[[i]])){
      errormessage <- paste0("Error in task ",i,": ",combine[[i]])
      warning(errormessage)
      writelog(paste0(gavi.dir,"gavi_log"),errormessage)
      errorcount <- errorcount + 1
      #remove from data
      combine[[i]] <- NULL
    }
  }
  #if(errorcount > 0){
  #	stop(paste0("There were ",errorcount," errors."))
  #}
  
  for(ii in 1:length(countries)){
    for(r in 1:runs){
      if(ii==1 & r==1){
        out_gavi_agestratified <- combine[[1]][[1]]
      } else {
        out_gavi_agestratified <- rbindlist(
          list(
            out_gavi_agestratified,
            combine[[ii]][[r]]
          )
        )
      }
    }
  }
  if(process_results){
    #clean environment
    if(remove_files){
      if(psa >0){
        for(p in 1:psa){
          if(p<10){
            p <- paste0("00",p)
          } else if(p<100){
            p <- paste0("0",p)
          }
          do.call(file.remove, list(list.files(paste0("./outcome/",save.scenario, "/",foldername,"/output/run",p,"/"), full.names = TRUE)))
        }
      } else {
        do.call(file.remove, list(list.files(paste0("./outcome/",save.scenario, "/",foldername,"/output/"), full.names = TRUE)))
      }
      if(Sys.info()[["sysname"]] == "Windows"){
        do.call(file.remove, list("./model/compiled/fort.6"))
      }
    }
    
    # write output
    if(psa >0){
      det_stoch <- "stochastic"
    } else {
      det_stoch <- "deterministic"
    }
    
    # make output similar to template
    # only report years to be reported
    out_gavi_agestratified <- out_gavi_agestratified[year %in% unique(template[, year])]
    out_gavi_agestratified[, country_name := unique(template[country == country_code, country_name]), by = country_code]
    out_gavi_agestratified[, "disease"] <- "Measles"
    colnames(out_gavi_agestratified)[3] <- "country"
    out_gavi_agestratified <- out_gavi_agestratified[, colnames(template), with=F]
    setorder(out_gavi_agestratified, country, year, age)
    
    # write data
    filename <- paste0(
      "outcome/",save.scenario, "/",
      foldername,
      "/",
      format(Sys.time(),format="%Y%m%d"),
      "_gavi_measles_",
      c("vax_none","vax_mcv1","vax_mcv1_mcv2")[vaccination+1],
      c("_sia_no","_sia_yes")[using_sia+1],
      "_",
      det_stoch,
      "_agestratified.csv"
    )
    z <- 1
    while(file.exists(filename)){
      warning(paste0(filename)," already exists, adding _",z)
      filename <- paste0(
        "outcome/",save.scenario, "/",
        foldername,
        "/",
        format(Sys.time(),format="%Y%m%d"),
        "_gavi_measles_",
        c("vax_none","vax_mcv1","vax_mcv1_mcv2")[vaccination+1],
        c("_sia_no","_sia_yes")[using_sia+1],
        "_",
        det_stoch,
        "_agestratified_",
        z,
        ".csv"
      )
      if(!file.exists(filename)){
        break
      }
      z <- z + 1
    }
    fwrite(
      out_gavi_agestratified,
      filename,
      row.names=FALSE
    )
  }
  
  # write file stochastic variables
  if(psa > 0){
    stoch_file <- paste0(format(Sys.time(),format="%Y%m"),"_out_measles_stochastic_variables.csv")
    if(file.exists(paste0("outcome/",stoch_file))){
      #read csv if file already exists
      stoch_file <- read.csv(paste0("input/",data_psa))
      #check if psa_var corresponds with psa
      if(nrow(stoch_file) != psa){
        stop(paste0("Number of runs in PSA output file is not the same as those specified! File is NOT overwritten, please delete or rename the old file if a new file needs to be generated"))
        writelog("gavi_log",paste0("Main; gavi.r aborted PSA error"))
        if(using_openmpi){
          mpi.quit()
        }
      }
    } else {
      #create csv if file does not exist
      print("Creating stochastic parameters file")
      stoch_file <- psa_var[,c("run_id","take1_input","take2_input","take3_input")]
      #get country specific CFR
      mortality <- as.numeric(psa_var[,"mortality_input"])
      for(c in 1:length(countries)){
        cfr <- Crit[which(Crit$country==as.character(countries[c])),"CFR"]
        stoch_file[,paste0(countries[c],"_CFR")] <- cfr*mortality
      }
      fwrite(
        stoch_file,
        paste0(
          "outcome/",
          format(Sys.time(),format="%Y%m"),
          "_out_measles_stochastic_variables.csv"
        ),
        row.names=FALSE
      )
    }
  }
  toc()
  
  
  # process results
  
  report_years <- sort(unique(template$year))
  ages         <- sort(unique(template$age))
  
  # get years that the model was run for (can be different from reporting years)
  years 		   <- sort(unique(as.numeric(coverage_routine[,year])))
  
  
  # set cfrs before 2000 to 2000
  cfr.year <- rbindlist(lapply(1980:1999, function(i) copy(cfr[Year ==2000,])[,Year := i]))
  cfr <- rbind(cfr.year, cfr)
  # set Kosovo mortaity to Serbia as it was once Serbia
  cfr.xk <- rbindlist(lapply("XK", function(i) copy(cfr[Code == "SRB",])[,Code := i]))
  cfr.xk[, Country := "Kosovo"]
  cfr <- rbind(cfr, cfr.xk)
  
  # expand cfr by age so that it uses the right value for <5 and 5-9 and 0 for >=10
  cfr.year.all <- rbindlist(lapply(0:100, function(i) copy(cfr)[, age:=i]))
  over10 <- T
  if (over10){
    cfr.year.all[, cfr.value :=  if (age < 5) under5 else if (age > 4 & age <10) over5 else 0, 
                 by = c("Code", "Year", "age") ] 
  } else {
    cfr.year.all[, cfr.value :=  if (age < 5) under5 else if (age > 4 ) over5, 
                 by = c("Code", "Year", "age") ]
  }
  cfr.year.all$Year = as.integer(cfr.year.all$Year)
  
  
  # ANALYSE RUNS 
  
  # get the names of all the files in all the subfolders the stochastic runs specified in scenarioname; separate cases and popsize files
  # will use those to read them all in at once and put in a data.table
  myfiles <- list.files(path = paste0(wd,"outcome/",save.scenario, "/",foldername,"/output/"), 
                        recursive = T, pattern = "cases", full.names = T)
  myfiles.popsize <- list.files(path = paste0(wd,"outcome/",save.scenario, "/",foldername,"/output"), 
                                recursive = T, pattern = "popsize", full.names = T)
  
  # read in all those files specified in myfiles and put them in a single data table
  all_cases  <- rbindlist(lapply(myfiles, function(fn, ...) {
    res <- withCallingHandlers(
      fread(fn, stringsAsFactors=F, check.names = F, fill = T, 
            col.names = as.character(c(0:100))),
      warning = function(w) { warning(w, fn); }
    )
    res[, country := gsub("^.+/(\\w+)_age.+$","\\1", fn) ]             # get the country code from the filename
    res[, year := years]                                               # add year of simulation (from coverage data file)
  })) 
  
  # switch back Kosovo to XK otherwise montagu returns error on upload
  all_cases[country == "XKX", country := "XK"]
  
  all_popsize  <- rbindlist(lapply(myfiles.popsize, function(fn, ...) {
    res <- fread(fn, stringsAsFactors=F, check.names = F, col.names = as.character(c(0:100)))
    res[, country := gsub("^.+/(\\w+)_age.+$","\\1", fn) ]             # get the country code from the filename
    res[, year := years]                                               # add year of model was run for (from coverage data file)
  })) 
  all_popsize[country == "XKX", country := "XK"]
  
  all_cases.m <- melt(all_cases, id.vars = c("year","country"), 
                      measure.vars = c(as.character(0:100)), 
                      variable.name = "age", value.name = "cases", variable.factor = F)
  
  all_popsize.m <- melt(all_popsize, id.vars = c("year", "country"),
                        measure.vars = c(as.character(0:100)), 
                        variable.name = "age", value.name = "cohort_size", variable.factor = F)
  # melt produces character variable when variable.factor is set to FALSE - change it to integer
  all_cases.m[, age := lapply(.SD, as.integer), .SDcols = "age"]
  all_popsize.m[, age := lapply(.SD, as.integer), .SDcols = "age"]
  # merge cases and cohort sizes
  all_runs = merge(all_cases.m, all_popsize.m, by = c("year", "age", "country"), all.x = T)
  
  # add country_name, life expectancy and disease (Measles) to match template file
  country_names <- unique(subset(template, select = c("country", "country_name")))
  c_names <- country_names$country_name; names(c_names) = country_names$country
  life.exp2 <- rbindlist(lapply(2100, function(i) copy(lexp[year == 2099])[, year:=i]))
  life.exp <- rbind(lexp, life.exp2)
  life.exp.all <- rbindlist(lapply(0:100, function(i) copy(life.exp)[, age:=i]))
  colnames(life.exp)[8] <- "LE"
  #all_runs[, c("country_name", "disease", "LE") := list(c_names[country], "Measles", life.exp[country])]
  
  all_runs[, c("country_name", "disease") := list(c_names[country],"Measles")]
  
  
  # merge mortality rates from CFR file
  all_runs <- merge(all_runs, subset(cfr.year.all, 
                                     select = c( "Code", "Year", "age", "cfr.value")),
                    by.x = c("country", "year", "age" ), by.y = c("Code", "Year", "age"), 
                    sort= F, all.x = T)
  
  all_runs <- merge (all_runs, 
                     subset (life.exp, select= c("country_code", "year", "LE")), 
                     by.x = c("country", "year"), 
                     by.y = c("country_code", "year"))
  
  # ----------------------------------------------------------------------------
  # add data column for remaining life expectancy
  all_runs <- lexp_remain [all_runs,
                           .(i.country, year, age, cases, cohort_size, country_name, disease, cfr.value, LE, value),
                           on = .(country_code = country,
                                  age_from    <= age,
                                  age_to      >= age,
                                  year         = year)
                           ]
  
  # rename column "i.country" to "country"
  setnames (x = all_runs, old = "i.country", new = "country")
  
  # save a copy of remaining life expectancy values in disease column
  all_runs [, disease := value]
  all_runs [, value   := NULL]
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  # MCV1 coverage
  coverage_routine_MCV1 <- coverage_routine [vaccine == "MCV1"]
  
  # add MCV1 column
  all_runs <- coverage_routine_MCV1 [all_runs, 
                                     .(i.country, i.year, age, cases, cohort_size, country_name, disease, cfr.value, LE, coverage),
                                     on = .(country_code = country,
                                            year         = year)
                                     ]
  
  # rename column "coverage" to "MCV1"
  setnames (x = all_runs, 
            old = c("i.country", "i.year", "coverage"), 
            new = c("country",   "year",   "MCV1"))
  # ----------------------------------------------------------------------------
  
  # calculate deaths and dalys
  all_runs[, deaths := cases*cfr.value]
  #all_runs[, dalys := (cases - deaths)*0.002 + deaths*(LE - 4)] # assumes all deaths are age 4, 
  all_runs [, dalys := (cases - deaths) * 0.002 + deaths * (LE - age)] # to account for actual age of death
  
  # ------------------------------------------------------------------------------
  # YLL calculation above is fine since most measles burden < 5 years
  # remaining life expectancy estimation can be improved by using remaining life
  # years by age, year and country
  # ------------------------------------------------------------------------------
  
  # OUTPUT RUNS
  
  # don't need all of these columns for VIMC, save only ones that are needed
  save.cols <- c(colnames(template))
  
  # ----------------------------------------------------------------------------
  save.cols <- c(save.cols, "MCV1")
  # ----------------------------------------------------------------------------
  
  output_runs <- subset(all_runs, year %in% report_years, select = save.cols)
  
  # write all output_runs to a file - need to specify the name (at the begining of the code), can be automated later
  fwrite(output_runs[order(country, year, age)], 
         paste0(wd,"outcome/",save.scenario, "/",foldername, "/", "burden_estimate_Measles-LSHTM-Jit-", scenario, ".csv"))
  
}

# clean environment
writelog ("gavi_log", paste0 ("Main; gavi.r finished"))

if (using_openmpi){
  mpi.quit()
}

print (Sys.time ())
toc ()  # stop timer 

