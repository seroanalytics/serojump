# Code for "Modelling serological data using reversible jump Markov chain Monte Carlo"

## Overview

This repository contains code and resources for a simulation recovery projected aimed at demostrating the ability of a reversible jump MCMC algorithm at recovering correlates of protection from simualted data. The project is a work in progress.

## Contents

-**src/**: This directory contains the source code for the rj-mcmc algorithm.
-**vingettes/**: This directory contains R Markdown files needed to recreate the analysis. Namely file `rjmc_serosim_ces_notd.Rmd`
-**R/**: This directory contains R files used to run the rjmcmc and post-process the results. The `serosim/` folder contains code to simuate the serological data used for the simulation recovery. The folders `cesNoCOP_notd` and `cesCOP_notd` are the folders containing the data used in the paper.

-**outputs/**: This directory contains all the outputs of the project.  `sim_data/`folder contains the data and figures for the simulated data.`fits/` folder contains all the results and figures from teh RJ-MCMC algorithm.The folders `cesNoCOP_notd` and `cesCOP_notd` are the folders containing the data used in the paper.
-**LICENSE**: This file contains the project's license information.

## Contact

For questions or inquiries about the project, please contact david.hodgson@lshtm.ac.uk
