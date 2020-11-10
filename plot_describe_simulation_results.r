### Describe and plot LCA simulation results
### R ARATHIMOS
### 220620

library(scales)
library(parallel)
library(foreach)
library(doParallel)
library(janitor)
library(broom)
library(reshape2)
library(ggplot2)

#install.packages("viridis")
library(viridis)

set.seed(123456)

data_dir <- "~/brc_scratch/data/genotype_lca/"
output_dir <- "~/brc_scratch/output/genotype_lca/"
scripts_dir <- "~/brc_scratch/scripts/genotype_lca/"

### sample sizes of different models

n_series <- c("10000", "20000", "40000", "80000", "160000")

bias_levels <- as.character(c("bias_level_1", "bias_level_2", "bias_level_3"))

### alpaha; proportions of each latent class in the cases subset - only 3 classes permited here

alpha = c(0.2, 0.3, 0.5) 

### LCA output 

#P      sizes of each latent class; equal to the mixing proportions in the basic latent class
#model, or the mean of the priors in the latent class regression model.

#probs.se     standard errors of estimated class-conditional response probabilities, in the same
#format as probs

#probs     is estimated class-conditional response probabilities.


#######################################

all_mod_errors <- list()
all_runtimes_list <- list()

for (num_i in 1:length(bias_levels)) {

	i <- as.character(bias_levels[num_i])

	res <- load(paste0(data_dir,"lca_workspace_effect_",i,".Rdata"))

	i <- as.character(bias_levels[num_i]) #there is an "i" in the imported rdata - overwrite

	print(paste0("Loaded results for ", i))

	proportion_errors <- list()
  runtime_list <- list()

	for (k in 1:length(n_series)) {

		### calculate deviation from true class proportion at population level

		maxdeviation <- max(all_lca[[k]]$P) - max(alpha)

		mindeviation <- min(all_lca[[k]]$P) - min(alpha)

		mediandeviation <- median(all_lca[[k]]$P) - median(alpha)

		errordeviation <- abs(mindeviation) + abs(maxdeviation) + abs(mediandeviation)

		percenterrordeviation <- errordeviation*100

    if(all_lca[[k]]$maxiter == all_lca[[k]]$numiter) { ml_fail <- "Yes" } else {ml_fail <- "No"}

		percenterrordeviationsample <- as.character(c(percenterrordeviation, n_series[k], ml_fail))

		proportion_errors[[k]] <- percenterrordeviationsample

    runtime_list[[k]] <- as.numeric((all_lca[[k]]$time))

	}
	
	proportion_errors <- as.data.frame(do.call(rbind, proportion_errors))
	
  all_runtimes_list[[i]] <- as.data.frame(do.call(rbind, runtime_list))

	names(proportion_errors) <- c("percentage_error","sample_size","ML_fail")

	#assign(paste0("proportion_errors_", i), proportion_errors)

	all_mod_errors[[i]] <- proportion_errors

}

all_mod_errors <- do.call(rbind, all_mod_errors)

all_runtime_df <- do.call(rbind, all_runtimes_list)

#tidy
all_mod_errors$bias  <- row.names(all_mod_errors)
all_mod_errors$`SNP strength` <- as.character(substr(all_mod_errors$bias, 1, nchar(all_mod_errors$bias) - 2 ))
all_mod_errors$`SNP strength`[all_mod_errors$`SNP strength`=="bias_level_1" ] <- "Simulation 1 - Strong"
all_mod_errors$`SNP strength`[all_mod_errors$`SNP strength`=="bias_level_2" ] <- "Simulation 2 - Moderate"
all_mod_errors$`SNP strength`[all_mod_errors$`SNP strength`=="bias_level_3" ] <- "Simulation 3 - Weak"

all_runtime_df$bias  <- row.names(all_runtime_df)
all_runtime_df$`SNP strength` <- as.character(substr(all_runtime_df$bias,1,nchar(all_runtime_df$bias)-2))
all_runtime_df$`SNP strength`[all_runtime_df$`SNP strength`=="bias_level_1" ] <- "Simulation 1 - Strong"
all_runtime_df$`SNP strength`[all_runtime_df$`SNP strength`=="bias_level_2" ] <- "Simulation 2 - Moderate"
all_runtime_df$`SNP strength`[all_runtime_df$`SNP strength`=="bias_level_3" ] <- "Simulation 3 - Weak"
all_runtime_df$sample_size <- n_series

### manual rescale times in days to hours - may need editing depending on output

all_runtime_df$V1[all_runtime_df$bias=="bias_level_3.4" | 
all_runtime_df$bias=="bias_level_3.5" | all_runtime_df$bias=="bias_level_2.5"] <- all_runtime_df$V1[all_runtime_df$bias=="bias_level_3.4" | 
all_runtime_df$bias=="bias_level_3.5" | all_runtime_df$bias=="bias_level_2.5"] * 24

###plot

pal <- viridisLite::viridis(3)
print(pal)

setwd(output_dir)

png(paste0('runtime_per_model_hours_plot.png'), unit='px', res = 300, width = 2000, height = 1500)

  ggplot(all_runtime_df) + 
    aes(as.numeric(sample_size), V1, color = `SNP strength`) + 
    geom_point() +
    geom_line() + 
    theme_bw() +
    ylab("Runtime (hours)") +
    xlab("Sample size model") +
    scale_color_manual(values = pal) +
    theme(axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) 

dev.off()


### plot of percentage error in class proportions

all_mod_errors$sample_size <- as.numeric(as.character(all_mod_errors$sample_size))

all_mod_errors$percentage_error <- as.numeric(as.character(all_mod_errors$percentage_error))

base_breaks <- function(n = 5){

  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }

}

png(paste0('error_class_proportions_plot.png'), unit = 'px', res = 300, width = 2000, height = 1500)

  ggplot(all_mod_errors) + 
    aes(sample_size, percentage_error, color = `SNP strength`) + 
    geom_point() +
    geom_line() + 
    theme_bw() +
    ylab("% error in proportion of classes") +
    xlab("Sample size") +
    scale_color_manual(values = pal) +
    theme(axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) +
    geom_text(data=subset(all_mod_errors, ML_fail == "Yes"), 
                aes(sample_size, percentage_error), label = "*", size=  12, 
                nudge_y = 0.1, show.legend = F, colour = "black" ) 
  

dev.off()

  
### plots of misclassified in to trueclass

all_trueclass_errors <- list()

for (num_i in 1:length(bias_levels)) {
  
  i <- as.character(bias_levels[num_i])
  
  res <- load(paste0(data_dir,"lca_workspace_effect_",i,".Rdata"))
  
  i <- as.character(bias_levels[num_i]) #there is an "i" in the imported rdata - overwrite
  
  print(paste0("Loaded results for ", i))
  
  trueclass_errors <- list()
  
  for (k in 1:length(n_series)) {
    
    ### calculate deviation from true class proportion at population level

    trueclass_deviation <- max(all_lca[[k]]$predclass) - max(alpha)

    if(all_lca[[k]]$maxiter == all_lca[[k]]$numiter) { ml_fail <- "Yes" } else {ml_fail <- "No"}

    percenterrordeviationsample <- as.character(c(percenterrordeviation, n_series[k], ml_fail))

    proportion_errors[[k]] <- percenterrordeviationsample
    
  }
  
  proportion_errors <- as.data.frame(do.call(rbind, proportion_errors))
  
  names(proportion_errors) <- c("percentage_error","sample_size","ML_fail")

  #assign(paste0("proportion_errors_", i), proportion_errors)

  all_mod_errors[[i]] <- proportion_errors
  
}

#end
