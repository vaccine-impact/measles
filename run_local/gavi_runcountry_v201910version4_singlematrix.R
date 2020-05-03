# Title:		Function ruiso3() for Measles Model Gavi
# Author:		Ed Jones, London School of Hygiene & Tropical Medicine. 
#           Edited by Kevin van Zandvoort and Petra Klepac, London School of Hygiene & Tropical Medicine
# Date:			11/25/2019
# Description:	Run country can be called in the foreach loop or in a for loop. It is wrapped in a function mainly to keep the environment clean.

# Changes 2019:
# 1) expanded matrix is rescaled to keep the total number of
# contacts in weekly age groups correct (checked - R0 of rescaled matrix is now
# the same as prior to expanding it to weekly ages)
# 2) projects the contact matrix to represent country's demography

runCountry <- function (
	#variables specific for loop
	ii,
	countries,
	years,
	
	#infection dynamic variables
	vaccination,
	using_sia,
	dinf,
	gamma,
	r0_basic,
	amplitude,
	take,
	degree,
	sia_method,
	
	#input data
	coverage_routine,
	coverage_sia,
	timeliness,
	contact,
	rnought,
	population,
	lexp,
	cfr,
	
	#dynaMICE model options
	tstep,
	save.scenario,
	foldername,
	measles_model,
	debug_model,
	debug_spinup,
	debug_age,
	debug_compartment,
	debug_country,
	debug_relative,
	debug_timepoints,
	
	#PSA options
	r,
	runs,
	psa,
	psa_var,
	
	#additional options
	process_results,
	run_model,
	remove_files = TRUE
) {
  
	iso3 <- countries[ii]
	#temporarily assign 3-letter-ISO code to Kosovo until Kosovo is assigned official ISO3-code
	if(iso3 == "XK"){
		fortran_country_code <- "XKX"
	} else {
		fortran_country_code <- iso3
	}
	
	if(psa > 0){
	  p <- r
	  if(p<10){
	    p <- paste0("00",p)
	  } else if(p<100){
	    p <- paste0("0",p)
	  }
	  psadir <- paste0("run",p,"/")
	} else {
	  psadir <- ""
	}
	
	if(run_model){
	  if(using_sia == 1){
	  	if(nrow(coverage_sia[country_code == iso3]) > 0){
	  		sia_a0 <- coverage_sia[country_code == iso3, a0]
	  		sia_year <- coverage_sia[country_code == iso3, year]
	  		sia_a1 <- coverage_sia[country_code == iso3, a1]
	  		sia_coverage <- coverage_sia[country_code == iso3, coverage]
	  		
	  		#for(s in 1:length((sia_coverage))){
	  		#	if(coverage_sia[country_code == iso3, Reached][s] > 0){
	  		#		#if reached is given, use that to estimate proportion vaccinated
	  		#		sia_coverage[s] <- coverage_sia[country_code == iso3, Reached][s] / sum(population[country_code == iso3 & year == sia_year[s] & age_from >= sia_a0[s] & age_from <= sia_a1[s], value])
	  		#	} else {
	  		#		#if reached is not given, assume 95% of target population has been vaccinated
	  		#		sia_coverage[s] <- 0.95*coverage_sia[country_code == iso3, Reached][s]/ sum(population[country_code == iso3 & year == sia_year[s] & age_from >= sia_a0[s] & age_from <= sia_a1[s], value])
	  		#	}
	  		#}
	  		sia_coverage[which(sia_coverage>1)] <- 0.95
	  	} else {
	  		sia_a0 			<- 0
	  		sia_a1 			<- 0
	  		sia_coverage	<- 0
	  		sia_year 		<- years[1]
	  	}
	  } else {
	  	sia_a0 			<- 0
	  	sia_a1 			<- 0
	  	sia_coverage	<- 0
	  	sia_year 		<- years[1]
	  }
	  
	  #t_cov is x year of model (model is ran for t_cov timesteps)
	  t_cov <- length(years) * tstep
	  t_end <- t_cov * 2
  
	  #t_run: for each year that is modelled,  get timepoint at which the year starts
	  t_run <- t_end-(length(years):1)*tstep+1
	  #get timepoints for year sia, assuming that year 1 of interest is the year 2000
	  t_sia=(t_run[1]+1000*(sia_year - years[1])) + c(1:length(sia_year))
	  if(using_sia == 0){
	  	t_sia <- t_sia*-1
	  }
	  	
	  print("Generating data for model...")
	  writelog("gavi_log",paste0(iso3, "; Run ",r,"/",runs,"; Start generating data"))
	  updateProgress(iso3,ii,runs,r,1)
	  
	  #Vaccine parameters
	  if(psa>0){
	  	take <- as.numeric(psa_var[r,c("take1_input","take2_input","take3_input")])
	  }
	  take1 <- c(take[1],0)[sia.method]    #vaccine take (dose 1, before age 1)
	  take2 <- c(take[2],0)[sia.method]    #vaccine take (dose 1, after age 1)
	  take3 <- c(take[3],0)[sia.method]    #vaccine take (dose 2)
	  degree1 <- c(0,degree[1])[sia.method]  #vaccine degree (dose 1, before age 1)
	  degree2 <- c(0,degree[2])[sia.method]  #vaccine degree (dose 1, after age 1)
	  degree3 <- c(0,degree[3])[sia.method]  #vaccine degree (dose 2)
	  
	  #country_specific contact matrix
	  r0_target <- rnought[country_code == iso3, r0]
	  q <- r0_target / r0_basic #proportionality factor (infectivity, underreporting)
	  gamma <- 1 / ( dinf * tstep / 365 ) #rate of losing infection
	  contact_day <- contact*q
	  contact_tstep <- contact_day * (365/tstep)
	  r0_tstep <- Re(eigen(contact_tstep, only.values=T)$values[1])
	  #country specific timeliness curve
	  country_timeliness <- timeliness[country_code == iso3 & !is.na(age), timeliness]
	  timeliness_ages <- timeliness[country_code == iso3 & !is.na(age), age]
	  
	  #Beta only a single file
	  s = 52 # number of finer stages within an age band (weekly ages, so 52)
	  jt = 3 # how many ages to expand to s (or to weeks)
	  beta_full <- matrix(0, ncol=254, nrow=254) # expanding to 100 years; contacts> 80 set to previously 0
	  
	  beta <- contact_tstep
	  r0_tstep <- Re(eigen(contact_tstep, only.values=T)$values[1]) # will need it to make sure R0 is the same after rescaling
	  
	  # contact_tstep_country <- (contact_tstep + contact_tstep*(1/pop.vector)%*%t(pop.vector))/2
	  # correction <- Re(eigen(contact_tstep, only.values=T)$values[1])/
	  #   Re(eigen(contact_tstep_country, only.values=T)$values[1])
	  # contact_tstep_country <-  contact_tstep_country*correction
	  # beta <- sweep(contact_tstep_country, 2, population[country_code == iso3 & year == y & age_to <= 80, value], "/")
	  #beta <- contact_tstep_country
	  #old script uses same popsize for ages 1-4 in 'pip'
	  #can be commented out when using interpolated pop
	  #beta <- sweep(contact_tstep, 2, c(rep(population[country_code == iso3 & year == y & age_to == 0, value],4), population[country_code == iso3 & year == y & age_to >3 & age_to <= 80, value]), "/")
	  
	  #create a new contact matrix for each age stratum in model
	  
	  #expand contact of first 3 age-strata with itself
	  
	  beta_full[(1:(s*jt)), (1:(s*jt))] <- expandMatrix(
	    A = beta[1:jt, 1:jt]/s,  # needs to be divided by 52 so that the mean total number of contacts stays the same,
	    expand_rows =  s, expand_cols =  s,
	    rescale_rows = FALSE, rescale_cols = FALSE
	  )
	  #expand contact for first 3 age-strata with all other contacts
	  beta_full[1:(s*jt),((s*jt)+1):(ncol(beta_full))] <- expandMatrix(
	    A = beta[1:jt,(jt+1):ncol(beta)]/s,
	    expand_rows = s, expand_cols = 1, 
	    rescale_rows = F, rescale_cols = F)
	  
	  beta_full[((s*jt)+1):(nrow(beta_full)), 1:(s*jt)] <- expandMatrix(
	    A = beta[(jt+1):nrow(beta),1:jt],
	    expand_rows = 1, expand_cols = s,
	    rescale_rows = F, rescale_cols = F)
	  
	  beta_full[((s*jt)+1):(nrow(beta_full)), ((s*jt)+1):(ncol(beta_full))] <-  beta[(jt+1):nrow(beta),(jt+1):ncol(beta)]
	  
	  # make the matrix symmetric
	  beta_full <- (beta_full +t(beta_full))/2
	  #	beta_full_country <- (beta_full + t(beta_full)*(pop.vector_full%*%t(1/pop.vector_full)))/2
	  
	  #make sure the R0 is what it should be
	  beta_full <- (r0_tstep/Re(eigen(beta_full, only.values=T)$values[1]))*beta_full
	  
	  #Reduce file-writes
	  #beta_full <- sweep(beta_full, 2, pop.vector_full, "/") # these are the infection rates from Wallinga et al
	  # here, rows should be divided by the appropriate population vector, not columns, but Fortran uses rows, rather than columns to calculate force of infection so it works out.
	  
	  beta_file <- paste0("outcome/",save.scenario, "/",foldername,"/input/",
	                      psadir,fortran_country_code,"_001beta.txt")
	  
	  #generate DynaMICE input files for each year
	  for (y in years){
	  	
	  	#divide number of effective contacts per timestep by the population size, account for 70+ being grouped in 70-80
	  	#this gives the proportion of individuals aged j contacted by someone aged i, per timestep
	  	
	  	#Change - used to use age-specific popsize, now pool all 70-80 together
	  	#beta <- sweep(contact_tstep, 2, c(population[country_code == iso3 & year == y & age_to <= 69, value], rep(sum(population[country_code == iso3 & year == y & age_to >= 79, value]), 11)), "/")
	    #beta <- sweep(contact_tstep, 2, c(population[country_code == iso3 & year == y & age_to <= 69, value], rep(sum(population[country_code == iso3 & year == y & age_to >= 70 & age_to <= 74, value]), 5), rep(sum(population[country_code == iso3 & year == y & age_to >= 75, value]), 6)), "/")
	  	
	  	#old script groups those aged 70-80, but division is by actual popsize
	    pop.vector <- population[country_code == iso3 & year == y, value]
	    
	    # first expand polymod matrix (contact_tstep) and poopulation vector and
	    # then divide by population sizes, otherwise it doesn't work.
	    
	    pop.vector_full <- c(rep(pop.vector[1:jt]/s, each = s), pop.vector[(jt +1):length(pop.vector)])
	    # change zero values to 1 to avoid division by 0 
	    pop.vector_full[pop.vector_full==0] <- 1
	   
	    
	  	if(vaccination>=1){
	  		#Maximum coverage can (obviously) only be 100%
	  		#To estimate proportion that is vaccinated at each week, we first calculate the number of individuals remaining susceptible
	  		#Then we calculate the number of individuals that should be vaccinated each week, in order to remain 1 - coverage susceptibles at the end of the timeliness data
	  		#In essence, this becomes the inverse of the cumulative timeliness curve
	  		cycov <- coverage_routine[country_code == iso3 & year == y & vaccine == "MCV1", coverage]/timeliness[country_code == iso3 & is.na(age), prop_final_cov]
	  		if(is.na(cycov)){
	  		  cycov <- 0
	  		}
	  	  country_year_timeliness_mcv1 <- 1 - min(
	  	    cycov,
	  			1
	  		) * country_timeliness
	  		
	  		country_year_timeliness_mcv1 <- (
	  		  country_year_timeliness_mcv1[1:(length(country_year_timeliness_mcv1)-1)] - country_year_timeliness_mcv1[2:(length(country_year_timeliness_mcv1))]
	  		)/(country_year_timeliness_mcv1[1:(length(country_year_timeliness_mcv1)-1)])
	  		
	  		#Timeliness is reported by week in the first year, and by month in the second year. Assume there is no vaccination in between
	  		country_year_timeliness_mcv1[is.na(country_year_timeliness_mcv1)] <- 0
	  		country_year_timeliness_mcv1[is.nan(country_year_timeliness_mcv1)] <- 0
	  		country_year_timeliness_mcv1_allages <- rep(0, 254)
	  		country_year_timeliness_mcv1_allages[round(timeliness_ages)] <- country_year_timeliness_mcv1
	  	} else {
	  		country_year_timeliness_mcv1_allages <- rep(0, 254)
	  	}
	  	if(vaccination==2){
	  		country_year_mcv2 <- coverage_routine[country_code == iso3 & year == y & vaccine == "MCV2", coverage]
	  	} else {
	  		country_year_mcv2 <- 0
	  	}
	  	if(is.na(country_year_mcv2)){
	  	  country_year_mcv2 <- 0
	  	}
	  	#First three ages are modelled in weekly strata
	  	pop_year <- c(
	  		rep(population[country_code == iso3 & year == y & age_to == 0, value]/52, 52),
	  		rep(population[country_code == iso3 & year == y & age_to == 1, value]/52, 52),
	  		rep(population[country_code == iso3 & year == y & age_to == 2, value]/52, 52),
	  		population[country_code == iso3 & year == y & age_to >= 3, value]
	  	)
	  	# change any zeros to 0 to avoid division by zero
	  	pop_year[pop_year==0] <- 1
	  	
	  	#write data for each year
	  	i <- which(years == y)
	  	
	  	if (i<10){
	  		#beta_file <- paste0("outcome/",save.scenario, "/",foldername,"/input/",
	  		#                    psadir,fortran_country_code,"_00",i,"beta.txt")
	  		dynamice_input_file <- paste0("outcome/",save.scenario, "/",foldername,"/input/",
	  		                              psadir,fortran_country_code,"_00",i,"measle_data.txt")
	  	} else if(i<100) {
	  		#beta_file <- paste0("outcome/",save.scenario, "/",foldername,"/input/",
	  		#                    psadir,fortran_country_code,"_0",i,"beta.txt")
	  		dynamice_input_file<- paste0("outcome/",save.scenario, "/",foldername,"/input/",
	  		                             psadir,fortran_country_code,"_0",i,"measle_data.txt")
	  	} else {
	  		#beta_file <- paste0("outcome/",save.scenario, "/",foldername,"/input/",
	  		#                    psadir,fortran_country_code,"_",i,"beta.txt")
	  		dynamice_input_file<- paste0("outcome/",save.scenario, "/",foldername,"/input/",
	  		                             psadir,fortran_country_code,"_",i,"measle_data.txt")
	  	}
	  	fwrite(
	  		as.list(
	  			as.data.table(
	  				format(
	  					beta_full,
	  					digits=14,
	  					scientific=FALSE
	  				)
	  			)
	  		),
	  		beta_file,
	  		quote=FALSE,
	  		row.names=FALSE,
	  		col.names=FALSE,
	  		sep=" "
	  	)
	  	
	  	dynamice_input_vector <- c(
	  		t_end,
	  		1,
	  		gamma,
	  		tstep,
	  		take1,
	  		take2,
	  		take3,
	  		degree1,
	  		degree2,
	  		degree3,
	  		length(sia_a0),
	  		sia_a0,
	  		r0_target,
	  		amplitude,
	  		1,
	  		length(sia_a1),
	  		sia_a1,
	  		t_end,
	  		t_cov,
	  		length(t_sia),
	  		t_sia,
	  		length(country_year_timeliness_mcv1_allages),
	  		country_year_timeliness_mcv1_allages,
	  		country_year_mcv2,
	  		length(sia_coverage),
	  		sia_coverage,
	  		length(t_run),
	  		t_run,
	  		format(
	  			pop_year,
	  			digits=14
	  		)
	  	)
	  	fwrite(
	  		list("a"=dynamice_input_vector),
	  		file=dynamice_input_file,
	  		col.names=FALSE
	  	)
	  }
	  
	  #run actual model
	  if(psa==0){
	  	p <- 0
	  } else {
	  	p <- r
	  }
	  writelog("gavi_log",paste0(iso3, "; Run ",r,"/",runs,"; Finished generating data"))
	  
	  #process debugging options
	  if(debug_country != "*" & !(iso3 %in% debug_country)){
	  	debug_debug			<- 0
	  	debug_compartments	<- 0
	  	debug_age			<- 0
	  	debug_timepoints	<- 0
	  	debug_relative		<- 0
	  } else {
	  	debug_debug 		<- as.integer(debug_spinup) + 2*as.integer(debug_model)
	  	debug_compartments	<- as.integer(debug_compartments)
	  	debug_relative		<- as.integer(debug_relative)
	  }
	  
	  #run the model
	  writelog("gavi_log",paste0(iso3, "; Run ",r,"/",runs,"; Start model"))
	  updateProgress(iso3,ii,runs,r,2)
	  if(Sys.info()[["sysname"]] == "Windows"){
	  	model <- paste0(
	  		'cd "',wd,'model/compiled/" & ',
	  		measles_model, " ",
	  		length(years), " ",
	  		tstep, " ",
	  		save.scenario, " ",
	  		foldername, " ",
	  		fortran_country_code, " ", 
	  		p, " ",
	  		"WIN", " ",
	  		debug_debug, " ",
	  		debug_compartments, " ",
	  		debug_age, " ",
	  		debug_timepoints, " ",
	  		debug_relative
	  	)
	  	model <- gsub("/","\\\\",model)
	  	shell(model, intern=TRUE)
	  } else {
	  	model <- paste0(
	  		'cd "',wd,'model/compiled/"; ./',
	  		measles_model, " ",
	  		length(years), " ",
	  		tstep, " ",
	  		save.scenario, " ",
	  		foldername, " ",
	  		fortran_country_code, " ", 
	  		p, " ",
	  		"LIN", " ",
	  		debug_debug, " ",
	  		debug_compartments, " ",
	  		debug_age, " ",
	  		debug_timepoints, " ",
	  		debug_relative
	  	)
	  	system(model, intern=TRUE)
	  }
	}
	writelog("gavi_log",paste0(iso3, "; Run ",r,"/",runs,"; Finished model"))
	updateProgress(iso3,ii,runs,r,3)
	
	#remove inputdata for this country
	#if files are removed around the same time by different workers, position may change
	if(remove_files){
		files <- list.files(
			paste0("./outcome/",save.scenario, "/",foldername,"/input/",psadir), full.names = TRUE
		)
		do.call(
			file.remove, list(
				files[
					grepl(
						fortran_country_code,
						files
					)
				]
			)
		)
	}
	if(process_results){
		#process model output
		cases <- fread(
			paste0(
				"outcome/",save.scenario, "/",
				foldername,
				"/output/",
				psadir,
				fortran_country_code,
				"_age_stratified_cases_byyear.txt"
			)
		)
		colnames(cases) <- as.character(c(1:ncol(cases))-1)
		cases[,"year"] <- years
		cases <- melt(cases, id.vars = "year", variable.factor=FALSE, value.factor = FALSE, variable.name="age", value.name="cases")
		class(cases$year) <- "numeric"
		class(cases$age) <- "numeric"
		class(cases$cases) <- "numeric"
		
		dat_cfr <- cfr[country_code == iso3, CFR]
		if(psa > 0){
			dat_cfr <- dat_cfr * as.numeric(psa_var[r,c("mortality_input")])
		}
		#deaths for those under 5: CFR/100 * cases
		#deaths for those over 5 and under 10: CFR/200 * cases
		#deaths for those over 10: 0
		cases[, "deaths"] <- 0
		cases[age < 10, "deaths"] <- cases[age < 10, cases] * (dat_cfr/200)
		cases[age < 5, "deaths"] <- cases[age < 5, cases] * (dat_cfr/100)
		
		for(y in years){
			dat_life_exp <- lexp[country_code == iso3 & year == y, value]
			z <- 1
			while(length(dat_life_exp) == 0){
			  dat_life_exp <- lexp[country_code == iso3 & year == years[which(years==y) - z], value]
			  z <- z+1
			  if(z == length(years)){
			    break
			  }
			}
			cases[year == y, "cohort_size"] <- population[country_code == iso3 & year == y, value]
			
			#Change 1 - Use life expectancy that changes by year (old lexp was fixed at life expectanct of 2017)
			
			#Change 2 - 
			#In the old version, assume that average age of the deaths was 4 years, and use lexp - 4 to estimate years of life lost for each death
			#In the new version, use the actual age of the child that died and life-expectancy of that year to estimate the effect
			
			#cases[year == y, "dalys"] <- (cases[year == y, cases] - cases[year == y, deaths]) * 0.002 + cases[year == y, deaths] * (dat_life_exp - cases[year == y, age])
			cases[year == y, "dalys"] <- (cases[year == y, cases] - cases[year == y, deaths]) * 0.002 + cases[year == y, deaths] * (dat_life_exp - 4)
		}
		
		cases[, "country_code"] <- iso3
		if(psa > 0){
			cases[, "run_id"] <- r
			cases <- cases[, c("run_id", "year", "age", "country_code", "cohort_size", "deaths", "cases", "dalys"), with=F]
		} else {
			cases <- cases[, c("year", "age", "country_code", "cohort_size", "deaths", "cases", "dalys"), with=F]
		}
		
		setorder(cases, year, age)
		
		writelog("gavi_log",paste0(iso3, "; Run ",r,"/",runs,"; Return results"))
		updateProgress(iso3,ii,runs,r,4)
		return(cases)
	}
}
