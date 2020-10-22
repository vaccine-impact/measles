To run measles model (execute):
- run_dynamice.R
(run_dynamice_201910gavi.R)

The above file sources:
- functions.R -- contains all requisite functions

Folders:
- R
    R code/programs

- vaccine_coverage
    VIMC vaccine coverage files (for no vaccination, copy another file with all rows and set coverage to 0)
- vaccine_coverage/scenarios
    includes subfolder scenarios which splits two files for routine and campaign immunisation

- burden_template
    disease burden template (VIMC format)

- input
    input files

- model
    fortran model code
       subfolder: fortran & compiled
       to compile fortran code: 
           gfortran -o ../compiled/vaccine2019_sia_singlematrix.exe vaccine2019_sia_singlematrix.F95 -fcheck=bounds -ffree-line-length-0

- outcome
    intermediary output folders and files

- central_burden_estimate
    subfolders Portnoy & Wolfson: central burden estimates for Portnoy and Wolfson CFRs
- central_burden_estimate/Wolfson
- central_burden_estimate/Portnoy
Note: To generate central burden estimates, set psa runs = 0

- plots
    subfolder with touchstone name contains diagnostic plots

- stochastic_burden_estimate
    subfolders Portnoy & Wolfson: stochastic burden estimates for Portnoy and Wolfson CFRs
- stochastic_burden_estimate/Wolfson
- stochastic_burden_estimate/Portnoy
Note: To generate stochastic burden estimates, set psa runs = 200 (or any number of runs)

-----------------------------------------------------------------------------------------------------------------------------------------------

README notes for emulator run of measles

1. Generate central estimates for all (10) scenarios.

2. Generate stochastic estimates (regular stochastic runs) for two countries (India and Nigeria) for campaign-default scenario
   India indicates high vaccination coverage (based on MCV1 coverage in 2015 > 85%) countries.
   Nigeria indicates low vaccination coverage (based on MCV1 coverage in 2015 < 85%) countries.

To emulate stochastic estimates (execute):
emulate_stochastic_estimate.R

3. Generate proportional variation of cases' estimates for the stochastic runs  based on a stochastic estimate for a 
   single country (India/Nigeria) and central estimate for a single country (India/Nigeria) for a single scenario.

4. Based on results of step 3, emulate stochastic estimates for all countries for all scenarios
   (i) apply proportional changes to cases
  (ii) cfr is calculated from central estimate file (deaths / cases)
 (iii) apply proportional changes to deaths
  (iv) calculate dalys

-----------------------------------------------------------------------------------------------------------------------------------------------



