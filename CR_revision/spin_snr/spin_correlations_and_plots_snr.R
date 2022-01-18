##############################################
########### Spin Test Correlations for #######
#########Intermodal Coupling Paper SNR########

####### Author: Erica Baller
#### Date: 3/9/2021

#######
##pre: right and left 10242 x 1000 matrices from matlab SpinPermuFS snr, mean coupling R& L
##post: right and left correlations
## uses: Takes output of spin test, and calcualted the correlation of snr to mean coupling
    #### 1) Read in snr spin results and mean maps
    #### 2) Convert -1s to NAs for medial wall
    #### 3) Foreach permutation (r and l separately), correlate with mean map
    #### 4) Store Rs and calc significance
    #### 5) Plot

### dependencies: ggplot2, bigmemory, vroom


#library(bigmemory.sri)
library(ggplot2)
library(tidyr)



#################
### set home directory
#homedir <- "/Users/eballer/BBL/imco/pmacs/PMACS_remote"
homedir <- "/project/imco"

source(paste0(homedir, "/baller/scripts/imco_functions.R"))

#initialize
hemis <- c("lh", "rh")
permNum <- 1000
models = c("snr_med_wall_-1")

################
### Read in matrices 
#mean coupling
  lh_mean_coupling <- read.table(paste0(homedir, "/baller/results/CR_revision/lh_mean_coupling_med_wall_-1.csv"))
  rh_mean_coupling <- read.table(paste0(homedir, "/baller/results/CR_revision/rh_mean_coupling_med_wall_-1.csv"))
  
 # all_mean_coupling <- c(lh_mean_coupling$V1, rh_mean_coupling$V1)

  
#make -1s into NA
  lh_mean_coupling_NA <- na_if(lh_mean_coupling, -1)
  rh_mean_coupling_NA <- na_if(rh_mean_coupling, -1)
  
  #actual maps
  
  lh_model <- read.table(paste0(homedir, "/baller/results/CR_revision/lh_snr_map_med_wall_-1.csv"))
  rh_model <- read.table(paste0(homedir, "/baller/results/CR_revision/rh_snr_map_med_wall_-1.csv"))
  
  #assign na
  lh_model_NA <- na_if(lh_model, -1)
  rh_model_NA <- na_if(rh_model, -1)
  
  
  #spins
  lh_spin <- t(read.table(paste0(homedir, "/baller/results/CR_revision/coupling_accuracy/spin_test_results/lh_spin_test_snr_med_wall_-1_output.csv"), sep = ","))
  rh_spin <- t(read.table(paste0(homedir, "/baller/results/CR_revision/coupling_accuracy/spin_test_results/rh_spin_test_snr_med_wall_-1_output.csv"), sep = ","))

  #assign NA where medial wall, -1 is medial wall
  lh_spin_NA <- na_if(lh_spin, -1)
  rh_spin_NA <- na_if(rh_spin, -1)
  
   #bring together, with original values as first column
  lh_act_results_and_spin <- cbind(lh_model_NA, lh_spin_NA)
  rh_act_results_and_spin <- cbind(rh_model_NA, rh_spin_NA)
  
  #correlations
  #go through each hemisphere, go through each perm, and go through each network
  
  lh_hemi_spin_calculations_r <- data.frame(vector(length = (permNum + 1)))
  rh_hemi_spin_calculations_r <- data.frame(vector(length = (permNum + 1)))
  
  for (hemi in hemis){

    for (perm in 1:(permNum + 1)){
     
        #correlation, use complete obs
        correlation <- eval(parse(text = as.character(paste0("cor(", hemi, "_act_results_and_spin[,", perm, "], ", 
                                                             hemi, "_mean_coupling_NA$V1, method = \"spearman\", use = \"complete.obs\")"))))
        
        #store r 
        eval(parse(text = as.character(paste0(hemi, "_hemi_spin_calculations_r[", perm, ",1] = ", correlation))))
     }
  }
  
  lh_p_permuted <- length(which(lh_hemi_spin_calculations_r[2:1001,1] > lh_hemi_spin_calculations_r[1,1]))/1000 
  rh_p_permuted <- length(which(rh_hemi_spin_calculations_r[2:1001,1] > rh_hemi_spin_calculations_r[1,1]))/1000  
  
  avg_p_permuted <- (lh_p_permuted + rh_p_permuted)/2
  
  save_df <- c(lh_p_permuted, rh_p_permuted, avg_p_permuted)
  names(save_df) <- c("lh_p_permuted", "rh_p_permuted", "avg_p_permuted")
  
  write.table(save_df, paste0(homedir, "/baller/results/CR_revision/meanxsnr_spin_results.csv"),quote = F, col.name = F)

  
  ##### SNR 25/35
  
  lh_model_25 <- read.table(paste0(homedir, "/baller/results/CR_revision/lh_snr_25_map_med_wall_-1.csv"))
  rh_model_25 <- read.table(paste0(homedir, "/baller/results/CR_revision/rh_snr_25_map_med_wall_-1.csv"))
  
  #assign na
  lh_model_25_NA <- na_if(lh_model_25, -1)
  rh_model_25_NA <- na_if(rh_model_25, -1)
  
  
  #spins
  lh_spin_25 <- t(read.table(paste0(homedir, "/baller/results/CR_revision/coupling_accuracy/spin_test_results/lh_spin_test_snr_25_med_wall_-1_output.csv"), sep = ","))
  rh_spin_25 <- t(read.table(paste0(homedir, "/baller/results/CR_revision/coupling_accuracy/spin_test_results/rh_spin_test_snr_25_med_wall_-1_output.csv"), sep = ","))
  
  #assign NA where medial wall, -1 is medial wall
  lh_spin_25_NA <- na_if(lh_spin_25, -1)
  rh_spin_25_NA <- na_if(rh_spin_25, -1)
  
  
  lh_model_35 <- read.table(paste0(homedir, "/baller/results/CR_revision/lh_snr_35_map_med_wall_-1.csv"))
  rh_model_35 <- read.table(paste0(homedir, "/baller/results/CR_revision/rh_snr_35_map_med_wall_-1.csv"))
  
  #assign na
  lh_model_35_NA <- na_if(lh_model_35, -1)
  rh_model_35_NA <- na_if(rh_model_35, -1)
  
  #spins
  lh_spin_35 <- t(read.table(paste0(homedir, "/baller/results/CR_revision/coupling_accuracy/spin_test_results/lh_spin_test_snr_35_med_wall_-1_output.csv"), sep = ","))
  rh_spin_35 <- t(read.table(paste0(homedir, "/baller/results/CR_revision/coupling_accuracy/spin_test_results/rh_spin_test_snr_35_med_wall_-1_output.csv"), sep = ","))
  
  #assign NA where medial wall, -1 is medial wall
  lh_spin_35_NA <- na_if(lh_spin_35, -1)
  rh_spin_35_NA <- na_if(rh_spin_35, -1)
  
  
###### do all spins together
  
  #df
  
  all_spin_calc_r <- data.frame(vector(length = 1001))
  
  #take concat of left and right
  all_mean_coupling <- c(lh_mean_coupling_NA$V1,rh_mean_coupling_NA$V1)
  all_actual_snr <- c(lh_model_NA$V1,rh_model_NA$V1)
  all_spin_snr <- rbind(lh_spin_NA, rh_spin_NA)

  big_df <- cbind(all_actual_snr, all_spin_snr)
  
  #run correlations
  for (perm in 1:1001) {
    all_spin_calc_r[perm,1] <- cor(big_df[,perm], all_mean_coupling, method = "spearman", use = "complete.obs")
  }

  permuted_p <- length(which(all_spin_calc_r[2:1001,1] > all_spin_calc_r[1,1]))/1000 #p=0.34, it is the top value

  hist(all_spin_calc_r[,1], main = paste0("Hist spin meanxsnr, R = ", round(all_spin_calc_r[1,1],3), " p= ", round(permuted_p,4)))
  
  
  ###### SNR 25/35
  
  ### 25
  all_spin_calc_r_25 <- data.frame(vector(length = 1001))
  
  #take concat of left and right
  all_actual_snr_25 <- c(lh_model_25_NA$V1,rh_model_25_NA$V1)
  all_spin_snr_25 <- rbind(lh_spin_25_NA, rh_spin_25_NA)
  
  big_df_25 <- cbind(all_actual_snr_25, all_spin_snr_25)
  
  #run correlations
  for (perm in 1:1001) {
    all_spin_calc_r_25[perm,1] <- cor(big_df_25[,perm], all_mean_coupling, method = "spearman", use = "complete.obs")
  }
  
  permuted_p_25 <- length(which(all_spin_calc_r_25[2:1001,1] > all_spin_calc_r_25[1,1]))/1000 #p=0.34, it is the top value
  hist(all_spin_calc_r_25[,1], main = paste0("Hist spin meanxsnr, R = ", round(all_spin_calc_r_25[1,1],3), " p= ", round(permuted_p_25,4)))
  
  ### 35
  
  all_spin_calc_r_35 <- data.frame(vector(length = 1001))
  
  #take concat of left and right
  all_actual_snr_35 <- c(lh_model_35_NA$V1,rh_model_35_NA$V1)
  all_spin_snr_35 <- rbind(lh_spin_35_NA, rh_spin_35_NA)
  
  big_df_35 <- cbind(all_actual_snr_35, all_spin_snr_35)
  
  #run correlations
  for (perm in 1:1001) {
    all_spin_calc_r_35[perm,1] <- cor(big_df_35[,perm], all_mean_coupling, method = "spearman", use = "complete.obs")
  }
  
  permuted_p_35 <- length(which(all_spin_calc_r_35[2:1001,1] > all_spin_calc_r_35[1,1]))/1000 #p=0.34, it is the top value
  hist(all_spin_calc_r_35[,1], main = paste0("Hist spin meanxsnr, R = ", round(all_spin_calc_r_35[1,1],3), " p= ", round(permuted_p_35,4)))
  
  
  
  
  ######## 40 and 50
  #actual maps
  
  lh_model_40 <- read.table(paste0(homedir, "/baller/results/CR_revision/lh_snr_40_map_med_wall_-1.csv"))
  rh_model_40 <- read.table(paste0(homedir, "/baller/results/CR_revision/rh_snr_40_map_med_wall_-1.csv"))
  
  #assign na
  lh_model_40_NA <- na_if(lh_model_40, -1)
  rh_model_40_NA <- na_if(rh_model_40, -1)
  
  
  #spins
  lh_spin_40 <- t(read.table(paste0(homedir, "/baller/results/CR_revision/coupling_accuracy/spin_test_results/lh_spin_test_snr_40_med_wall_-1_output.csv"), sep = ","))
  rh_spin_40 <- t(read.table(paste0(homedir, "/baller/results/CR_revision/coupling_accuracy/spin_test_results/rh_spin_test_snr_40_med_wall_-1_output.csv"), sep = ","))
  
  #assign NA where medial wall, -1 is medial wall
  lh_spin_40_NA <- na_if(lh_spin_40, -1)
  rh_spin_40_NA <- na_if(rh_spin_40, -1)
  
  
  lh_model_50 <- read.table(paste0(homedir, "/baller/results/CR_revision/lh_snr_50_map_med_wall_-1.csv"))
  rh_model_50 <- read.table(paste0(homedir, "/baller/results/CR_revision/rh_snr_50_map_med_wall_-1.csv"))
  
  #assign na
  lh_model_50_NA <- na_if(lh_model_50, -1)
  rh_model_50_NA <- na_if(rh_model_50, -1)
  
  
  #spins
  lh_spin_50 <- t(read.table(paste0(homedir, "/baller/results/CR_revision/coupling_accuracy/spin_test_results/lh_spin_test_snr_50_med_wall_-1_output.csv"), sep = ","))
  rh_spin_50 <- t(read.table(paste0(homedir, "/baller/results/CR_revision/coupling_accuracy/spin_test_results/rh_spin_test_snr_50_med_wall_-1_output.csv"), sep = ","))
  
  
  #assign NA where medial wall, -1 is medial wall
  lh_spin_50_NA <- na_if(lh_spin_50, -1)
  rh_spin_50_NA <- na_if(rh_spin_50, -1)
  
  #df
  
  
  ### 40
  all_spin_calc_r_40 <- data.frame(vector(length = 1001))
  
  #take concat of left and right
  all_actual_snr_40 <- c(lh_model_40_NA$V1,rh_model_40_NA$V1)
  all_spin_snr_40 <- rbind(lh_spin_40_NA, rh_spin_40_NA)
  
  big_df_40 <- cbind(all_actual_snr_40, all_spin_snr_40)
  
  #run correlations
  for (perm in 1:1001) {
    all_spin_calc_r_40[perm,1] <- cor(big_df_40[,perm], all_mean_coupling, method = "spearman", use = "complete.obs")
  }
  
  permuted_p_40 <- length(which(all_spin_calc_r_40[2:1001,1] > all_spin_calc_r_40[1,1]))/1000 #p=0.34, it is the top value
  hist(all_spin_calc_r_40[,1], main = paste0("Hist spin meanxsnr, R = ", round(all_spin_calc_r_40[1,1],3), " p= ", round(permuted_p_40,4)))
  
  ### 50
  
  all_spin_calc_r_50 <- data.frame(vector(length = 1001))
  
  #take concat of left and right
  all_actual_snr_50 <- c(lh_model_50_NA$V1,rh_model_50_NA$V1)
  all_spin_snr_50 <- rbind(lh_spin_50_NA, rh_spin_50_NA)
  
  big_df_50 <- cbind(all_actual_snr_50, all_spin_snr_50)
  
  #run correlations
  for (perm in 1:1001) {
    all_spin_calc_r_50[perm,1] <- cor(big_df_50[,perm], all_mean_coupling, method = "spearman", use = "complete.obs")
  }
  
  permuted_p_50 <- length(which(all_spin_calc_r_50[2:1001,1] > all_spin_calc_r_50[1,1]))/1000 #p=0.175, P = 0.107, it is the top value
  hist(all_spin_calc_r_50[,1], main = paste0("Hist spin mean x snr, R = ", round(all_spin_calc_r_50[1,1],3), " p= ", round(permuted_p_50,4)))
       #####)
  
  
  ### plots
#  plot(all_actual_snr, all_mean_coupling, pch=20, main = "SNR vs Mean Couplin, R = 0.34, P < 0.001", xlab = "SNR", ylab = "Mean Coupling (T)", abline(cor(all_mean_coupling,all_actual_snr, use="complete.obs", method = "spearman"), col = "red"))
 # visreg(lm(all_mean_coupling ~ all_actual_snr),main = "SNR vs Mean Coupling, R=0.34, P < 0.001", xlab = "SNR", ylab = "Mean Coupling (T)", points=list(col="black"))
#  visreg(lm(all_mean_coupling ~ all_actual_snr_50),main = "SNR>=50 vs Mean Coupling, R=0.17, P = 0.11", xlab = "SNR", ylab = "Mean Coupling (T)", points=list(col="black"))