### LCA of simulated genotypes
### R ARATHIMOS
### 100720
 
#library(depmixS4)
library(poLCA)
library(parallel)
library(foreach)
library(doParallel)
library(janitor)
library(broom)

set.seed(123456)

data_dir <- "~/brc_scratch/data/genotype_lca/"
output_dir <- "~/brc_scratch/output/genotype_lca/"
scripts_dir <- "~/brc_scratch/scripts/genotype_lca/"

n_series <- c("10000","20000","40000","80000","160000")

### aprox effect sizes for a continuous measure of depression from this paper:
### https://www.nature.com/articles/ng.3552/tables/1

#snp_latent_trait_associations <- c("bias_level_1","bias_level_2","bias_level_3") #strength of the snp association with the latent trait (effect size)

### pull in args

args	<- 	commandArgs(trailingOnly = TRUE)

snp_latent_trait_associations	<- 	toString(args[1]) #Phenotypic exposure or outcome of interest

latent_vars <- "trueclass"

for(k in 1:length(snp_latent_trait_associations)) {

  assoc_effect <- snp_latent_trait_associations[k]

  dfx <- readRDS(file = paste0(data_dir,"simulated_poLCA_data_",assoc_effect,".rds")   )  #load simulated data
  
  ### split "control" class - i.e. largest class, from other classes

  indx1 <- tail(names(sort(table(dfx$trueclass))),1)
  print(paste0("The largest class and therefore the control class to be excluded is ", indx1))
  
  dfx <- dfx$dat[which(dfx$trueclass!=indx1),]
  dim(dfx)
  
  ### initiate cluster for parallel process

  cl <- makeCluster(detectCores()-4)
  registerDoParallel(cl)

  all_lca <- foreach(j=1:length(n_series), .packages='poLCA') %dopar% {
    
    set.seed(123456)

    ### split latent traits from SNPs

    geno_df <-dfx[,! colnames(dfx) %in% latent_vars]

    ### convert df of genotypes to factors

    indx <- sapply(geno_df, is.numeric)

    geno_df[indx] <- lapply(geno_df[indx], function(x) as.factor(x))
    
    ### take a sample of size of nid or the max sample size available if 160000 specified

    nid <- as.numeric(as.character(n_series[j]))

    if (nid==160000) {

      geno_df2 <- geno_df[sample(row.names(geno_df), size=NROW(geno_df)),] } else {
      geno_df2 <- geno_df[sample(row.names(geno_df), size=nid),] 	

    }
    
    ### create formula for LCA

    measurevar <- "1"
    groupvars  <- names(geno_df)
    
    f <- as.formula(paste0("cbind(", paste(groupvars, collapse=","),")~", measurevar))
    
    ### run LCA

    res <- poLCA(formula=f, data=geno_df2, nclass = 3, maxiter = 4000, nrep=20)
    
    assign(paste0("res_",j), res)
    
    return(get(paste0("res_",j)))

  }

  stopImplicitCluster()

  for(i in 1:length(n_series)) {

    capture.output(print(all_lca[[i]]), file = paste0(output_dir,"lca_output_N_",n_series[i],"_assoc_",assoc_effect,".txt"))

  }

  print(all_lca)

  save.image(paste0(output_dir,"lca_workspace_effect_",assoc_effect,".Rdata"))

}

#end
