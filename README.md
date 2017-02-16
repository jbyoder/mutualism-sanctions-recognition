Theory and individual-based simulations of mutualistic symbiosis mediated by sanctions and partner recognition
==============================================================================================================

Files in this repository support the manuscript 

> Yoder JB and P Tiffin. Sanctions, partner recognition, and variation in mutualistic symbiosis

Contents now reflect revisions for resubmission to the *American Naturalist*, with extensive changes to the simulation scripts.

Contents
--------

The contents of the directories in this repository are as follows

### data

Contains summarized output from the simulation scripts, which underlies all figures in the manuscript. Data files supporting the current revision have the date-tag `20170117` in their names.

### figures

Contains all figures presented in the manuscript

### Mathematica

Contains Mathematica notebooks (and PDF versions of the same) presenting the three analyitic models evaluated in the manuscript.

### ms

Contains the manuscript text, in both a Markdown/LaTeX flat text file and a typeset PDF document generated from the same.

### Rscripts

Contains the following R scripts

- `sim_functions.R` --- Functions needed in the individual-based simulations described in the manuscript, called using `source()` in the primary simulation script.
- `sim_par_example.R` --- File containing parameters for the simulations, as described in the manuscript, which can be modified to change the metapopulation structure and population genetic parameters of the simulations. This is formatted as R script, and implemented in the primary simulation script by calling with `source()`.
- `simulations_par_soft2.R` --- Primary simulation script, which runs replicate simulations of the different coevolutionary scenarios described in the manuscript. Relies on `sim_functions.R` and `sim_par_example.R`. To run, use the following format at the command line:

> R --vanilla --args {sim} {rep nos} {record interval} {par file} {output file} < Rscripts/simulations_par.R

Where 

- `sim` = `N` (neutral) | `S` (sanctions) | `R` (recognition) | `SR` (both)
- `rep nos` = range of replicate i.d. numbers, as `1-100`; `101-200`, &c
- `record interval` = time between write-out, in simulation generations
- `par file` = path to file with parameters for the sims, such as `sim_par_example.R`
- `output file` = path and (base) filename for files recording simulation results

The script will read simulation parameters from the parameter file and write results to the path and file specified as the output file. 

