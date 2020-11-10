# Latent class analysis of simulated genotypes


## Analysis overview

** This simulation analysis uses the simdata function of the poLCA package to simulate "SNPs" as categorical variables, where the allele frequency between groups of cases of the 4 simulated classes and that of controls differs based on a pre-defined probability matrix. 
**

The parameters that are set to vary in the simulations are hard-coded in the scripts. These are; sample size (numbers of individual cases), and difference in allele frequency between classes and therefore strength of association with the underlying trait.

The proportions of cases to controls is 50:50. The proportion of cases in each of the 3 classes in the overall sample is  0.25 / 0.15 / 0.1

The number of SNPs simulated is 300, of which half (50%) are simulated to be pleiotropic.
  
The core package used for these simulations is poLCA

```r
install.packages("poLCA")

```
