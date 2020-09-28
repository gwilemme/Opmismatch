# Optimal Taxation to Correct Job Mismatching
Data and computer code associated to the article `Optimal Taxation to Correct Job Mismatching', by [Guillaume Wilemme](http://www.gwilemme.com/) (Univ. of Leicester).

The model is calibrated and simulated with Julia, using moments built with Stata from the Current Population Survey.


## Stata: From the raw tables to the moments
The raw data are the Basic Monthly CPS public files from January 2017 to December 2018.

The data are downloaded in the `data/raw` folder, going from `jan17pub.dat` to `dec18pub.dat`.
The dictionnary file used is `cpsbjan2015.dct` and is downloaded in the `code/stata` folder.

In the `code/stata` folder, the file `raw2ready.do` sequentially opens the .dat tables, select the key variables and save the new tables in the folder `data/ready` with the names `cpsjan17.dta` to `cpsdec18.dta`.

In the same folder, the file `ready2stats.do` generates a file `stats.csv` in the `data/moments` folder.
The `stats.csv` file reports the average wage and the share of employment for the 22 occupations.
The file `ready2transitions.do` generates a file `rates.csv` in the `data/moments` folder.
The `rates.csv` file reports the average monthly transition rates between the 22 occupations and unemployment.

These two csv files are in this repository.



## Julia: calibration and simulation
The Julia code runs with Julia 1.0.5 on Linux. It uses the NBInclude package to display the code as a Jupyter notebook. Each code is commented with equations if needed. If GitHub does not manage to display a notebook, I suggest using the url on GitHub with the viewer https://nbviewer.jupyter.org/.


### Routines and module Opmismatch
All The routines are in the `code/julia/src` folder.
First the file `Opmismatch.ipynb` generates a Julia module, calling all the other source codes of the folder.


`structure.ipynb` defines particular types for the allocation object or the tax object.
`data.ipynb` opens the tables `stats.csv` and `rates.csv` and define the moments in Julia.
`equilibrium.ipynb` provides routines to compute the equilibrium.
`calibration.ipynb` provides routines to calibrate the model by distance-to-moment minimisation.
`desc.ipynb` provides functions to compute statistics from the simulaions.

### Initialisation
The file `startup.ipynb` initialises the parameters before calibration and before generating tables.

### Calibration
The calibration of the model is the most time-consuming. Each of the nine versions in the paper and appendix follows the same calibration strategy. The three first steps serve as determining the best initial values for the coefficients, The fourth step calibrates the full model with these initial values. 

First, a segmented model with many parameters fixed is calibrated (`calib_step1.ipynb`).
This step takes less than 3 hours per version on the high-performance computer of the University of Leicester.
Second, the mismatching coefficients are calibrated (`calib_step2.ipynb`) for one hour per version.
Third, the occupation-specific frictions are calibrated (`calib_step3.ipynb`) for one hour per version.

The last step (`calib_step4.ipynb`) calibrates all the parameters together for 19 hours per version. Extending this computing time does not significatnly improve the convergence.

Notice that these codes are run as .jl files. The calibrated are saved in the `code/julia/calibrated` folder as `step4_v1.csv` to `step4_v9.csv`.


### Simulations and final tables
The file `stats_main.ipynb` generates the two tables in the main text of the article, for the three main versions. It takes 3 hours in total on the high-performance computer.
The file `stats_secondary.ipynb` generates the two tables in the appendix for the six other versions.
It takes 6 hours in total on the high-performance computer.


