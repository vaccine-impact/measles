
------------------------------------------------------------------------------------------------------
File 

 - lexp.R
     create remaining life expectancy table for each year matching values from 5-year intervals for 
       country, age, and year (both sexes)

  - allison_out_final.R
       Convert output files (burden estimates) provided by Allison with age-varying cfrs to vimc format.
       Deaths and dalys will need to be computed using cases, cfrs and remaining life expectancy

Note: vaccine coverage files
      "...mcv2..." from vimc download is renamed to "...mcv12..."
          "...mcv2-bestcase..." from vimc download is renamed to "...mcv12-bestcase..."
          "...mcv2-default..."  from vimc download is renamed to "...mcv12-default..."
------------------------------------------------------------------------------------------------------


Kevin:
Everything is on Dropbox, ./VIMC/measles/201910version4/run_local
13:17
There is a new model, ./model/fortran/vaccine2019_sia_singlematrix.F95 which will need to be recompiled and referred to in the Rscript
13:18
key is to read in the runcountry_*singlematrix.R file in the respective gavi_setup_.R
13:18
Sorry these need to be * stars
13:18
But got converted in bold-text
13:19
But yes the runcountry_X_singlematrix.R file should be sourced (edited) 
13:19
In the gavi_setup_X file that is being used
13:19
And in the gavi_runcountry_X script, the correct model needs to be specified, as well as the correct setup
13:21
Also attaching the files here
--------------------------------------------------------------------------------------------------------------------

README

This folder has all the files needed to run the stochastic and deterministic runs locally. The thing that needs to be changed manually is “index” parameter that sets the scenario to run (on line 53), in the file 

gavi_v201910version4_det_local_scenario1.R (for deterministic runs)
 
or

gavi_v201910version4_stoch_local_scenario1.R (for stochastic runs). 

Other files are called within that file and there is no need to manipulate them. 

To change the list of countries to run, change the line 204, to specify missing countries.  E.g.

countries <- missing_countries 
# missing countries is a list of country codes that need to be run.


These R codes call fortran code “vaccine2019_sia”. Before you can run this, you need to compile the fortran code. To do that, open the terminal and go to the folder with the fortran code (../model/fortran) and compile using the command: 

gfortran -o ../compiled/vaccine2019_sia vaccine2019_sia.F95 -fcheck=bounds -ffree-line-length-0
gfortran -o ../compiled/vaccine2019_sia_singlematrix.exe vaccine2019_sia_singlematrix.F95 -fcheck=bounds -ffree-line-length-0

All the other input files are in the subfolder "input" (cfr_new_noage.csv, demographic files, life expectancies, coverage scenarios, new R0 values, new timeliness values, new templates, polymod matrix extending to all ages). 

___________________________________________________________________________

Changes to this code:
#       a) add scenario folder pathname for clarity and for easier running on the cluster (fortran code changed too)
#       b) change contact matrix to fix augumented mixing in <3 yo and to reflect country specific mixing (rescale polymod by local population structure to keep contacts reciprocal)
#          - uses newly rescaled contact matrix to expand polymod to yearly ages up to 100
#       c) time-varying CFRs (postprocessing)
#       d) MCV1 and MCV2 dependency - give MCV2 to those who received MCV1 for MCV2 coverage <= MCV1, if MCV2 coverage > MCV1 distribute the remainder randomly (modified fortran code)
#       e) MCV1 and SIA dependency - use linear regression of data from Portnoy et al 2018 to inform the proportion of zero-dose children reached by SIA campaign (modified fortran code)
#       f) due to issues with age_range_verbatim field in the scenario (coverage) input files, these runs use age_first and age_last




