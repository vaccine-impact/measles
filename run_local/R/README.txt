To run measles model (execute):
- run_dynamice.R

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

- plots
    subfolder with touchstone name contains diagnostic plots



