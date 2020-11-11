######### Simulate genotypes, describe and run LCA
######### Ryan Arathimos
######### 06/06/2020

#library(depmixS4)
library(poLCA)
library(parallel)
library(foreach)
library(doParallel)
library(janitor)
library(broom)
library(MASS)

set.seed(123456)

data_dir <- "~/brc_scratch/data/lca_simulation/"
output_dir <- "~/brc_scratch/output/lca_simulation/"
scripts_dir <- "~/brc_scratch/scripts/lca_simulation/"

### sample sizes to simulate

n_series <- c("10000","20000","40000","80000","160000")

### simulateion information - N SNPs, pleotropy percentage, N classes etc

nsnp <- 300 #vector indicating how many snps

percent_snps_pleiotropic <- 50 # % in whole numbers

nid <- 320000 #create x samples/individuals

nclass <- 4

alpha = c(0.5, 0.25, 0.15, 0.1) ## proportions of each latent class - prevalence of each latent class in the sample

### bias levels names - bias to allele frequencies between "cases" and "controls", hence determine streght of SNP-trait assoc

bias_levels <- c("bias_level_1","bias_level_2","bias_level_3")

######################

nresp <- rep(3, nsnp)

all_probs <- list()

nsnps_not_pleiotropic <- (nsnp*percent_snps_pleiotropic)/100
nsnps_pleiotropic <- (nsnp*(100-percent_snps_pleiotropic))/100

#########################################################################
################## CREATE MATRIX OF PROBABILITIES #######################

### initialise empty lists for storage

bias_level1 <- list()
bias_level2 <- list()
bias_level3 <- list()

probs <- list()

### create probability matrix for each bias level which will form a simulated dataset

for(i in 1:nsnps_not_pleiotropic) {
  
  ### the pleiotropy mutiplier here is 1

  pleio_multiplier <- c(1,1,1) 

  ### pleiotropic snps - create matrix of probabilities
  
  class1_probs <- matrix(c(0.5, 0.3, 0.2), ncol=3) #the control class in all bias levels
  
  bias_class2 <- c(-0.05, 0.025, 0.025) * pleio_multiplier[1] #probability bias for each response in class 2
  bias_class3 <- c(0.05, -0.025, -0.025) * pleio_multiplier[2] #these are adjustments to the probability of a SNP having a '1' genotype
  bias_class4 <- c(0.05, -0.025, -0.025) * pleio_multiplier[3] #these are adjustments to the probability of a SNP having a '1' genotype
  
  class2_probs <- class1_probs + bias_class2
  class3_probs <- class1_probs + bias_class3
  class4_probs <- class1_probs + bias_class4
  
  bias_level1[[i]] <- rbind(class1_probs, class2_probs, class3_probs, class4_probs)

  bias_class2 <- c(-0.02, 0.01, 0.01)  * pleio_multiplier[1] #probability bias for each response in class 2
  bias_class3 <- c(0.02, -0.01, -0.01)  * pleio_multiplier[2]
  bias_class4 <- c(0.02, -0.01, -0.01)  * pleio_multiplier[3] #these are adjustments to the probability of a SNP having a '1' genotype
  
  class2_probs <- class1_probs + bias_class2
  class3_probs <- class1_probs + bias_class3
  class4_probs <- class1_probs + bias_class4
  
  bias_level2[[i]] <- rbind(class1_probs, class2_probs, class3_probs, class4_probs)  
  
  bias_class2 <- c(-0.01, 0.005, 0.005) * pleio_multiplier[1] #probability bias for each response in class 2
  bias_class3 <- c(0.01, -0.005, -0.005) * pleio_multiplier[2]
  bias_class4 <- c(0.01, -0.005, -0.005) * pleio_multiplier[3] #these are adjustments to the probability of a SNP having a '1' genotype
  
  class2_probs <- class1_probs + bias_class2
  class3_probs <- class1_probs + bias_class3
  class4_probs <- class1_probs + bias_class4
  
  bias_level3[[i]] <- rbind(class1_probs, class2_probs, class3_probs, class4_probs)
  
}

all_probs_not_pleiotropic <- c(list(bias_level1), list(bias_level2), list(bias_level3))

### now create matrix of probabilities for non pleiotropic SNPs

for(i in 1:nsnps_pleiotropic) {
    
    ### the pleiotropy multiplier here is randomised across the 3 classes

    pleio_multiplier <- list(c(0,0,1), c(0,1,0), c(1,0,0))[[sample(c(1:3), 1)]]
  
    ### pleiotropic snps - create matrix of probabilities
    class1_probs <- matrix(c(0.5,0.3,0.2),ncol=3)
    
    bias_class2 <- c(-0.05, 0.025, 0.025) * pleio_multiplier[1] #probability bias for each response in class 2
    bias_class3 <- c(0.05, -0.025, -0.025) * pleio_multiplier[2] #these are adjustments to the probability of a SNP having a '1' genotype
    bias_class4 <- c(0.05, -0.025, -0.025) * pleio_multiplier[3] #these are adjustments to the probability of a SNP having a '1' genotype
    
    class2_probs <- class1_probs + bias_class2
    class3_probs <- class1_probs + bias_class3
    class4_probs <- class1_probs + bias_class4
    
    bias_level1[[i]] <- rbind(class1_probs, class2_probs, class3_probs, class4_probs)
    
    bias_class2 <- c(-0.02, 0.01, 0.01)  * pleio_multiplier[1] #probability bias for each response in class 2
    bias_class3 <- c(0.02, -0.01, -0.01)  * pleio_multiplier[2]
    bias_class4 <- c(0.02, -0.01, -0.01)  * pleio_multiplier[3] #these are adjustments to the probability of a SNP having a '1' genotype
    
    class2_probs <- class1_probs + bias_class2
    class3_probs <- class1_probs + bias_class3
    class4_probs <- class1_probs + bias_class4
    
    bias_level2[[i]] <- rbind(class1_probs, class2_probs, class3_probs, class4_probs)  
    
    bias_class2 <- c(-0.01, 0.005, 0.005) * pleio_multiplier[1] #probability bias for each response in class 2
    bias_class3 <- c(0.01, -0.005, -0.005) * pleio_multiplier[2]
    bias_class4 <- c(0.01, -0.005, -0.005) * pleio_multiplier[3] #these are adjustments to the probability of a SNP having a '1' genotype
    
    class2_probs <- class1_probs + bias_class2
    class3_probs <- class1_probs + bias_class3
    class4_probs <- class1_probs + bias_class4
    
    bias_level3[[i]] <- rbind(class1_probs, class2_probs, class3_probs, class4_probs)
    
}  
  
all_probs_pleiotropic <- c(list(bias_level1), list(bias_level2), list(bias_level3))
  
all_probs <- mapply(c, all_probs_not_pleiotropic, all_probs_pleiotropic, SIMPLIFY = FALSE) 
  
print("Created probability matrix...!")


#########################################################################
###################### SIMULATE GENOTYPE DATA ###########################

### simulate data using poLCA simdata function

for (k in 1:length(bias_levels)) {
  
  probs <- all_probs[[k]]
  
  lca1 <- poLCA.simdata(N = nid, probs = probs, nclass = nclass, ndv = nsnp, nresp = 3,  P=alpha)

  geno_df <- lca1$dat #SNP data

  #rename columns to SNPs

  names(geno_df) <- paste0("SNP",seq(1:length(geno_df)))

  #convert simulated SNPs to numerics        
  indx <- sapply(geno_df, is.factor)
  geno_df[indx] <- lapply(geno_df[indx], function(x) as.numeric(as.character(x)) )
  
  geno_df_trueclass <- geno_df
  geno_df_trueclass$trueclass <- lca1$trueclass
  
  #save SNP matrix
  saveRDS(geno_df, file = paste0(data_dir,"simulated_poLCA_SNP_data_bias_level_",k,".rds"))

  #save truth
  saveRDS(lca1, file = paste0(data_dir,"simulated_poLCA_data_bias_level_",k,".rds"))

  print("Saved data!")
  
  ######################################################

  ### describe simulated dataset and check with regression models and correlations matrices the simulated vars
      
    freq_list <- list()
      
    for(i in names(geno_df)) {
      
        tb_class <- data.frame()
        
          for(class_k in c(1:4)) {

          tbx <- tabyl(geno_df_trueclass[[i]][geno_df_trueclass$trueclass==class_k])

          tb_class <- rbind(tb_class, tbx[,3])
          
          }
        names(tb_class) <- c("V1","V2","V3")
        freq_list[[i]] <- tb_class
        
      }
      
    freq_data <- do.call(rbind, freq_list)
     
    write.csv(freq_data, file=paste0(output_dir,"freq_data_snps_bias_level_",k,".csv"))
    
    print("Done with freqs!")
    
    ### regress and check P, beta and SE + Rsquared

    dfx = data.frame(geno_df,lca1$trueclass)

    dfx$trueclass <- as.factor(dfx$`lca1.trueclass`)

    snps_loop <- names(geno_df)[1:150]
      
    ### enter parallel session
      
    cl <- makeCluster(detectCores()-3)
    registerDoParallel(cl)

    print("Entering parallel loop!")
    
    all_res <- foreach(i=1:length(snps_loop), .packages="MASS" ) %dopar% {
      
        j <- snps_loop[i]
        
        all_latent_res <- data.frame()

          w <- "trueclass"

  	      formulax <- paste0(j," ~ ",w)
  	      
  	      mod <- glm(formulax, data=dfx)

  	      res_i <- data.frame(SNP=as.character(j),
  	      					  Latent_trait=as.character(w),
  	                  BETA_class2=coef(summary(mod))[2,1],
  	                  SE_class2=coef(summary(mod))[2,2],
  	                  Pvalue_class2=coef(summary(mod))[2,4],
                      BETA_class3=coef(summary(mod))[3,1],
                      SE_class3=coef(summary(mod))[3,2],
                      Pvalue_class3=coef(summary(mod))[3,4],
  	      					  BETA_class3=coef(summary(mod))[4,1],
  	      					  SE_class3=coef(summary(mod))[4,2],
  	      					  Pvalue_class3=coef(summary(mod))[4,4],
  	      					  N=nobs(mod) )

    }

  stopImplicitCluster()

  print("Simulated data!")

  gc()

  ### convert to dataframe and save  

  all_res <- as.data.frame(do.call(rbind, all_res))

  write.csv(all_res, file = paste0(output_dir, "regressions_snps_on_trueclass_bias_level_", k, ".csv"))

}

warnings()

#end

##########################################################################
################# POST SIMULATION DESCRIPTIVE STATS ######################

### calculate average class SNP frequencies across three datasets

average_class_freq <- data.frame()

### for non pleiotropic SNPs, i.e. 1:600

for (k in 1:length(bias_levels)) {

	frqs1 <- read.csv(paste0(output_dir,"freq_data_snps_bias_level_",k,".csv"))[1:600,]

  temp2 <- strsplit(as.character(frqs1$X), "[.]")
  
  frqs1$SNP <- lapply(temp2, `[[`, 1)
  
  frqs1$class <- lapply(temp2, `[[`, 2)

  for(class_k in c(1:4)) {

  	c1 <- round(mean(frqs1$V1[frqs1$class==class_k]),3)
  	c2 <- round(mean(frqs1$V2[frqs1$class==class_k]),3)
  	c3 <- round(mean(frqs1$V3[frqs1$class==class_k]),3)

  	average_class_freq <- rbind(average_class_freq, c(c1,c2,c3))

     }

}

print(average_class_freq)

### for pleiotropic SNPs

average_class_freq <- data.frame()

for (k in 1:length(bias_levels)) {

	frqs1 <- read.csv(paste0(output_dir,"freq_data_snps_bias_level_",k,".csv"))[601:1200,]
  
  temp2 <- strsplit(as.character(frqs1$X), "[.]")
  
  frqs1$SNP <- lapply(temp2, `[[`, 1)
  
  frqs1$class <- lapply(temp2, `[[`, 2)

  for(class_k in c(1:4)) {

  	c1 <- round(mean(frqs1$V1[frqs1$class==class_k]),3)
  	c2 <- round(mean(frqs1$V2[frqs1$class==class_k]),3)
  	c3 <- round(mean(frqs1$V3[frqs1$class==class_k]),3)

  	average_class_freq <- rbind(average_class_freq, c(c1,c2,c3))

     }

}

print(average_class_freq)
