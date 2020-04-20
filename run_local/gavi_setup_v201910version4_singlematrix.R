# load libraries
library (tictoc)
library (data.table)
library (foreach)

# global options in which R computes and displays its results.
options ("scipen" = 999, 
         "digits" = 14)

# small function to write message a to logfile logname
writelog <- function(logname, x){
	write(
		paste0(
			format(Sys.time(), "%Y/%m/%d %H:%M:%S"),
			" ",
			x
		),
		file=logname,
		append=TRUE
	)
}

# function that prints a message, and writes to the logfile
logwarning <- function(logname, x){
	writelog(logname, x)
	print(x)
}

# easily expand a Matrix (i.e. duplicate columns and rows in an AxB matrix - and rescale where applicable)
# easy and transparant way to redistribute contact matrix - set default rescale to FALSE
expandMatrix <- function(
  A, expand_rows=1, expand_cols=1, rescale_rows=F, rescale_cols=F
){
  if(!is.matrix(A)){
	stop("A is not a matrix")
  }
  matvals <- numeric(0)
  rows <- nrow(A)
  cols <- ncol(A)
  for(c in 1:cols){
    matvals <- c(
      matvals,
      rep(
        A[,c],
        expand_cols,
        each = expand_rows
      )
    )
  }
  B <- matrix(matvals, nrow=rows*expand_rows, ncol = cols*expand_cols)
  if(rescale_rows & rescale_cols){
    B <- B/(expand_rows*expand_cols) 
  } else if(rescale_rows){
    B <- B/expand_rows
  } else if(rescale_cols){
    B <- B/expand_cols
  }
  return(B)
}

# update file showing progress (temporarily disabled)
updateProgress <- function(country,cnt,runs,r,status){
	##wait until gavi_progress is yours (so cannot be overwritten by two parallel scripts at the same time)
	#while(
	#	!(file.rename("gavi_progress",paste0("gavi_progress_locked_",cnt,"_",r)))
	#){}
	##update line
	#progress <- readLines("gavi_progress_locked",-1)
	#progress[((cnt-1)*runs+r)] <- paste0(
	#	country,
	#	" ",
	#	paste0(c(rep(0,(nchar(runs)-nchar(r))),r),collapse=""),
	#	" ",
	#	status
	#)
	#writeLines(progress,"gavi_progress_locked")
	#file.rename(paste0("gavi_progress_locked_",cnt,"_",r),"gavi_progress")
}

# load correct libraries when on cluster (using open MPI)
if("cluster_using_openmpi" %in% commandArgs()){
	using_openmpi <- TRUE
	require("doMPI")
	require("parallel")
	cl <- startMPIcluster()
	registerDoMPI(cl)
	print("using openMPI")
} else {
	using_openmpi <- FALSE
	print("no openMPI")
}
#load required libraries
if(!using_openmpi){
	require("doParallel")
}




# function for each country
source (paste0(wd,"gavi_runcountry_v", version, "_singlematrix.R"))
#clean environment
if(file.exists("gavi_progress")){
	file.remove("gavi_progress")
}

# are we running a PSA - a deterministic or "stochastic" model?		
if(psa==0){
	det_stoch <- "deter" #deterministic
	runs<-1
} else {
	det_stoch <- "stoch" #stochastic
	runs<-psa
}

# create folder with correct name for in- and output data if not yet exists
# typically foldername should not exist - but may come in handy when only processing results
if( !exists("foldername_analysis") ){
	foldername <- paste0(
		format(Sys.time(),format="%Y%m%d"),
		"_v",
		vaccination,
		"_s",
		using_sia,
		"_",
		det_stoch
	)
	dir.create(
		file.path(
			paste0(
				getwd(),
				"/outcome/",save.scenario, "/",
				foldername
			)
		), recursive = T
	)
	dir.create(
		file.path(
			paste0(
				getwd(),
				"/outcome/",save.scenario, "/",
				foldername,
				"/input"
			)
		), recursive = T
	)
	dir.create(
		file.path(
			paste0(
				getwd(),
				"/outcome/",save.scenario, "/",
				foldername,
				"/output"
			)
		), recursive = T
	)
	if(psa > 0){
		for(p in 1:psa){
			if(p<10){
				p <- paste0("00",p)
			} else if(p<100){
					p <- paste0("0",p)
			}
			dir.create(
				file.path(
					paste0(
						getwd(),
						"/outcome/",save.scenario, "/",
						foldername,
						"/input/",
						paste0("run",p)
					)
				), recursive = T
			)
			dir.create(
				file.path(
					paste0(
						getwd(),
						"/outcome/",save.scenario, "/",
						foldername,
						"/output/",
						paste0("run",p)
					)
				), recursive = T
			)
		}
	}
} else {
  foldername <- foldername_analysis
}

timelinessCov <- function(
	tim_shape,
	target_age,
	max_cov
){
	coverage_byage_cumulative <- c(
      rep(0,(target_age-1)),
      rep(max_cov,(254-target_age+1))
    )
	beta_function <- cumsum(
      dbeta(
        (1:52/52),
        4,
        tim_shape
      )
    )/max(
      cumsum(
        dbeta(
          (1:52/52),
          4,
          tim_shape
        )
      )
    )
	coverage_byage_cumulative[39:(39+51)] <- max_cov*beta_function
    coverage_byage=rep(0,254)
    coverage_byage[2:length(coverage_byage_cumulative)] <- (
      coverage_byage_cumulative[2:length(coverage_byage_cumulative)]
      -coverage_byage_cumulative[1:(length(coverage_byage_cumulative)-1)]
    )/(
      1-coverage_byage_cumulative[1:(length(coverage_byage_cumulative)-1)]
    )
    coverage_byage[which(is.na(coverage_byage))] <- 1
	return(coverage_byage)
}