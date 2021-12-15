####### Violin Plots ###


#### Author: Erica Baller

### 3/16/2021

###pre: requires that you have run the spin tests for your models, including the positive and negative spin tests
###post: 2 violin plots per model. Violin plot x axis - yeo networks, y axis - proportion of vertices within networks. dotted lines - mean from permutations
   #### plot 1 will contain a thick black line for the *actual values from your analysis so you can compare how far above or below your value is from permuted
   #### plot 2 will contain red lines detailing where the actual values from your analysis from the POSITIVE domain, blue for NEGATIVE
   #### You will also get a table that tells you the p value (uncorrected) for each network within each model.
   #### for my final plot, I added a * to these manually
### uses: Creates violin plots and analyses to help make sense of the results from permutation spin tests. The question we are asking is:
    ### what is the likelihood that the number of vertices within a network is significant rather than due to chance
### dependencies: Using R 3.6.3 but any R will do. Libraries to include listed below

library(tidyr)
library(ggplot2)
library(reshape)

#################
### set home directory
#homedir <- "/Users/eballer/BBL/imco/pmacs/PMACS_remote"
homedir <- "/project/imco"
source(paste0(homedir, '/baller/scripts//imco_functions.R'))

#outdir_name <- "spin_stats_gam_sex_lm_exec_acc_all_pos_neg"
outdir_name <- "spin_stats_gam_sex_gam_exec_acc_all_pos_neg"

models_for_stats = c("gam_age", "pos_gam_age", "neg_gam_age", "gam_sex", "pos_gam_sex", "neg_gam_sex","lm_exec_accuracy", "pos_lm_exec_accuracy", "neg_lm_exec_accuracy", "gam_exec_accuracy", "pos_gam_exec_accuracy", "neg_gam_exec_accuracy", "mean_coupling")
models_for_plots = c("gam_age", "gam_sex", "lm_exec_accuracy", "gam_exec_accuracy")

network_names <- c("VIS", "MOT", "DA", "VA", "LIM", "FP", "DM")
num_spins = 2000
#get Yeo colors from function, these values were set manually - this can be obtained through https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011
   #I typed the RGB values into a rgb -> hex converter and stored the values here. Works!

yeo_colors <- get_yeo7_colors()

#for storing statistics at the end
stats <- data.frame(matrix(nrow = 7, ncol = length(models_for_stats)))
names(stats) <- models_for_stats
row.names(stats) <- network_names

for (model in models_for_stats) {
  
  if (model == "mean_coupling") {
    lh_spin_df <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/lh_spin_test_mean_coupling.csv"), sep = ","))
    rh_spin_df <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/rh_spin_test_mean_coupling.csv"), sep = ","))
  } else {
    lh_spin_df <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/lh_spin_test_", model, "_proportions.csv"), sep = ","))
    rh_spin_df <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/rh_spin_test_", model, "_proportions.csv"), sep = ","))
  }
  
  
  #take means of left and right
  actual_results <- (lh_spin_df[,1] + rh_spin_df[,1])/2
  
  #dataframes for all spins
  spin_without_target_col <- cbind(lh_spin_df[,2:1001],rh_spin_df[,2:1001])
  
  ####### Stats ########
  for (i in 1:7){
    #store p values
    #equivalen to: stats$model[i] <- (length(which(spin_without_target_col[i,] > actual_results[i]))/2000)
    eval(parse(text=as.character(paste0("stats$", model, "[", i ,"] <- (length(which(spin_without_target_col[",i, ",] > actual_results[", i, "]))/", num_spins, ")"))))
  }
}

print(stats)
write.csv(stats, paste0(homedir, "/baller/results/coupling_accuracy/spin_test_results/", outdir_name, ".csv"))

violin_plot_means_mean_coupling(homedir = homedir, network_names = network_names, num_spins = num_spins)
violin_plot_pos_and_neg_lines(homedir = homedir, models = models_for_plots, network_names = network_names, num_spins = num_spins)
violin_plot_means(homedir = homedir, models = models_for_plots, network_names = network_names, num_spins = num_spins)

#violin_plot_pos_and_neg_lines_with_color_gradation(homedir = homedir, models = models_for_plots, network_names = network_names, num_spins = num_spins)
