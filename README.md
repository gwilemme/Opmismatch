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


 The code is written with Jupyter notebooks.
The code runs with Julia 1.1.

The folder `code` contains the Jupyter notebooks. GitHub can directly display these notebooks.
Sometimes GitHub does not manage to display a notebook. 
We suggest using the url on GitHub with the viewer https://nbviewer.jupyter.org/.

