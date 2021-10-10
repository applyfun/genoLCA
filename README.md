# Latent class analysis of simulated genotypes


## Overview

**This simulation analysis uses the simdata function of the poLCA package to simulate "SNPs" as categorical variables.
The simulation assumes three underlying traits or classes that are associated with the simulated SNP genotypes. 
The allele frequency between groups of cases comprising the three simulated classes and that of controls differs based on a constructed probability matrix.**

The parameters that are set to vary in the simulations are hard-coded in the scripts. These are; sample size (numbers of individual cases), and difference in allele frequency between classes and therefore strength of association with the underlying trait.

The proportions of cases to controls is 50:50. The proportion of cases in each of the 3 classes in the overall sample is  0.25 / 0.15 / 0.1

The number of SNPs simulated is 300, of which half (50%) are simulated to be pleiotropic.
  
## Setup

The core package used for these simulations is poLCA

```r
install.packages("poLCA")
library(poLCA)

```

### Simulate the data

To simulate the genotype matrix, run the  simulate_poLCA_data.r file as a job using the run_simulate_data.sh submission file.

### Run the latent class analysis

The models run as separate jobs on the cluster. Create the job submission files for each one using the make_jobs.bash script, which uses the simulation-template.sh file to create separate job submission files. You can submit those jobs to the scheduler using the run_jobs.sh script.

### Describe and plot the results

The plot_describe_simulation_results.r summarises and plots the results of the latent class models. It is designed to be run interactively.
