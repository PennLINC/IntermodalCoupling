##################################################
### Accuracy Scripts 02/09/2021 ###
##################################################

####################
##### Summary ######
####################

#input: asc files, and pnc demographics, cnb, clinical files
#output: 10242 length vector csvs with T and/or p values for vertex-wide regression
#uses: goes vertex by vertex and does regression (coupling by age, sex, cognition, etc), pulls out T and p from these values and sticks it in a vector. The vector can then be used for visualization in matlab
#dependencies: R (3.6.3 is my current default in pmacs)

####################
###  Libraries   ###
####################

library(mgcv)
library(dplyr)
library(ggplot2)
library(visreg)

#####################################################################################
####              Makes the 831x10242 matrices, both left and right              ####
#####################################################################################

#set home directory, switch this depending on whether running from PMACS or from home directory
#homedir <- "/Users/eballer/BBL/imco/pmacs/PMACS_remote/"
homedir <- "/project/imco/"

# read in demos
subjDemos <- read.csv(paste0(homedir, "/baller/subjectLists/n831_alff_cbf_finalSample_imageOrder.csv"))

#some verification preprocessing
subjDemos$sex <- as.factor(subjDemos$sex)
subjDemos$race <- as.factor(subjDemos$race)
subjDemos$race2 <- as.factor(subjDemos$race2)

#add osex category for use in gam later
subjDemos$osex <- ordered(subjDemos$sex)

#add psych bifactor scores
psych <- read.csv(paste0(homedir, "/pnc/clinical/n1601_goassess_itemwise_bifactor_scores_20161219.csv"), header = TRUE)

#remove 4factorv2 from title
names(psych) <-gsub("_4factorv2", "", names(psych))

#merge
subjDemos <- merge(subjDemos, psych, by = "bblid")

#cognitive data
cog <- read.csv(paste0(homedir, "/pnc/cnb/n1601_cnb_factor_scores_tymoore_20151006.csv"))
accuracy <- subset(cog, select = c("bblid", "F1_Exec_Comp_Res_Accuracy"))

#merge
subjDemos <- merge(subjDemos, accuracy, by = "bblid")

#drop the scanid.y
subjDemos <- subset(subjDemos, select = -scanid.y)

#rename scanid.x to scanid
names(subjDemos) <- gsub("scanid.x", "scanid", names(subjDemos)) 

#make list of bblid/scanid
bblid_scanid <- paste0(subjDemos$bblid, "_", subjDemos$datexscanid)

#####################
##### Left Side #####
#####################

#initiate matrix for storage
numrows <- dim(subjDemos)[1]
lh_matrix <- matrix(nrow = numrows, ncol = 10242)

#go through each subject, grab 5th column in asc, transpose and stick in matrix
# output is 831 x 10242 matrix
for (subj in 1:831) {
  
  bblid <- subjDemos$bblid[subj]
  datexscanid <- subjDemos$datexscanid[subj]
  file_path <- paste0(homedir, "/couplingSurfaceMaps/alffCbf/lh/stat/", bblid, "_", datexscanid, "_lh.coupling_coef_alff_cbf.fwhm15.fsaverage5.asc")
  alff_data <- read.table(file_path, stringsAsFactors = FALSE)
  lh_matrix[subj,] <- t(alff_data$V5)
  
}

#append with demographics
subjDemos_with_lh_matrix <- cbind(subjDemos, lh_matrix)

#write output
write.table(subjDemos_with_lh_matrix, file = paste0(homedir, "/baller/results/coupling_accuracy/subjDemos_with_lh_", numrows, "x10242.csv"), sep = ",")

#####################
#### Right Side #####
#####################
#initiate matrix for storage
numrows <- dim(subjDemos)[1]
rh_matrix <- matrix(nrow = numrows, ncol = 10242)

#go through each subject, grab 5th column in asc, transpose and stick in matrix
# output is 831 x 10242 matrix
for (subj in 1:831) {
  
  bblid <- subjDemos$bblid[subj]
  datexscanid <- subjDemos$datexscanid[subj]
  file_path <- paste0(homedir, "/couplingSurfaceMaps/alffCbf/rh/stat/", bblid, "_", datexscanid, "_rh.coupling_coef_alff_cbf.fwhm15.fsaverage5.asc")
  alff_data <- read.table(file_path, stringsAsFactors = FALSE)
  rh_matrix[subj,] <- t(alff_data$V5)
  
}

#append with demographics
subjDemos_with_rh_matrix <- cbind(subjDemos, rh_matrix)

#write output
write.table(subjDemos_with_rh_matrix, file = paste0(homedir, "/baller/results/coupling_accuracy/subjDemos_with_rh_", numrows, "x10242.csv"), sep = ",")


####-----------------------------------------------------------------------------####
####---------------------------End of Part 1- Making matrices---_----------------####
####-----------------------------------------------------------------------------####

#####################################################################################
####                   Run regression, both left and right                       ####
#####################################################################################

#make easier to reference names
lh_imco <- subjDemos_with_lh_matrix #can also read directly from files if you'd like
rh_imco <- subjDemos_with_rh_matrix


#####################################################
#                       lm/gams                     #
#####################################################

#initialize vectors for models

hemis <- c("lh", "rh") #hemispheres
models <- c("age", "sex", "exec_accuracy")

coeffs <- c("p", "t") #p or t value
corrs <- c("uncor", "fdr") #correction

for (hemi in hemis){
  for (model in models) {
    for (coeff in coeffs) {
      for (corr in corrs) {
        vector_init_cmd <- paste0(hemi, "_gam_", model, "_", coeff, "_", corrs, " <- vector(length = 10242)")
        print(vector_init_cmd)
        eval(parse(text=as.name(vector_init_cmd)))
      }
    }
  }
}


#make linear models as well
#lm_models <- c("age", "accuracy", "speed", "efficiency")
lm_models <- c("age", "sex","exec_accuracy")
for (hemi in hemis) {
  for (model in lm_models) {
    for (coeff in coeffs) {
      for (corr in corrs) {
        vector_init_cmd <- paste0(hemi, "_lm_", model, "_", coeff, "_", corrs, "<- vector(length= 10242)")
        eval(parse(text=as.name(vector_init_cmd)))
      }
    }
  }
}


#######################
######## Left #########
#######################

#get # of items in df for calculation of column)
numcolumns <- dim(lh_imco)[2]
#run gams models and store info in respective vectors
for (i in 1:10242) {
  curcol = (numcolumns - 10242 + i) # will start you counting at the right part of the df
  age_sex_model <- gam(lh_imco[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                         osex + s(ageAtScan1, k = 4, fx = T), data=lh_imco)
  
  ## accuracy
 
  exec_accuracy_model <- gam(lh_imco[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                               osex + s(ageAtScan1, k = 4, fx = T) + F1_Exec_Comp_Res_Accuracy, data=lh_imco)
  
  #lm
  age_lm_model <- lm(lh_imco[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1, 
                     data=lh_imco)
  sex_lm_model <- lm(lh_imco[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + osex, 
                     data=lh_imco)
  
  ## acc
  
  exec_accuracy_lm_model <- lm(lh_imco[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1 + F1_Exec_Comp_Res_Accuracy, 
                               data=lh_imco)
  

  
  #put pvalue in it's appropriate lm
  lh_gam_age_p_uncor[i] <- summary(age_sex_model)$s.table[1,4] #smooth term for ageAtScan1
  lh_gam_sex_p_uncor[i] <- summary(age_sex_model)$p.table[4,4] #linear term
  
  lh_gam_exec_accuracy_p_uncor[i] <- summary(exec_accuracy_model)$p.table[5,4] #accuracy term
 
  #lm to assess directionality
  lh_lm_age_p_uncor[i] <- summary(age_lm_model)$coeff[4,4]
  lh_lm_sex_p_uncor[i] <- summary(sex_lm_model)$coeff[4,4]
  
  lh_lm_exec_accuracy_p_uncor[i] <- summary(exec_accuracy_lm_model)$coeff[5,4]
  
  #put tvalue in it's appropriate lm
  lh_gam_age_t_uncor[i] <- summary(age_sex_model)$s.table[1,3] #smooth term for ageAtScan1
  lh_gam_sex_t_uncor[i] <- summary(age_sex_model)$p.table[4,3] #linear term
  
  lh_gam_exec_accuracy_t_uncor[i] <- summary(exec_accuracy_model)$p.table[5,3] #accuracy term

  #lm to assess directionality
  lh_lm_age_t_uncor[i] <- summary(age_lm_model)$coeff[4,3]
  lh_lm_sex_t_uncor[i] <- summary(sex_lm_model)$coeff[4,3] #linear term
  
  lh_lm_exec_accuracy_t_uncor[i] <- summary(exec_accuracy_lm_model)$coeff[5,3]

}

#####################
###### RIGHT ########
#####################
#get # of items in df for calculation of column)
numcolumns <- dim(rh_imco)[2]
#run gams models and store info in respective vectors
for (i in 1:10242) {
  curcol = (numcolumns - 10242 + i) # will start you counting at the right part of the df
  age_sex_model <- gam(rh_imco[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                         osex + s(ageAtScan1, k = 4, fx = T), data=rh_imco)
  
  ## accuracy
  exec_accuracy_model <- gam(rh_imco[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                               osex + s(ageAtScan1, k = 4, fx = T) + F1_Exec_Comp_Res_Accuracy, data=rh_imco)
  
 
  
  
  #lm
  age_lm_model <- lm(rh_imco[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1, 
                     data=rh_imco)
  sex_lm_model <- lm(rh_imco[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + osex, 
                     data=rh_imco)
  
  ## acc

  exec_accuracy_lm_model <- lm(rh_imco[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1 + F1_Exec_Comp_Res_Accuracy, 
                               data=rh_imco)
  
  
  
  #put pvalue in it's appropriate lm
  rh_gam_age_p_uncor[i] <- summary(age_sex_model)$s.table[1,4] #smooth term for ageAtScan1
  rh_gam_sex_p_uncor[i] <- summary(age_sex_model)$p.table[4,4] #linear term
  
  rh_gam_exec_accuracy_p_uncor[i] <- summary(exec_accuracy_model)$p.table[5,4] #accuracy term
  
  #lm to assess directionality
  rh_lm_age_p_uncor[i] <- summary(age_lm_model)$coeff[4,4]
  rh_lm_sex_p_uncor[i] <- summary(sex_lm_model)$coeff[4,4]
  
  rh_lm_exec_accuracy_p_uncor[i] <- summary(exec_accuracy_lm_model)$coeff[5,4]
  
  #put tvalue in it's appropriate lm
  rh_gam_age_t_uncor[i] <- summary(age_sex_model)$s.table[1,3] #smooth term for ageAtScan1
  rh_gam_sex_t_uncor[i] <- summary(age_sex_model)$p.table[4,3] #linear term
  
  rh_gam_exec_accuracy_t_uncor[i] <- summary(exec_accuracy_model)$p.table[5,3] #accuracy term

  #lm to assess directionality
  rh_lm_age_t_uncor[i] <- summary(age_lm_model)$coeff[4,3]
  rh_lm_sex_t_uncor[i] <- summary(sex_lm_model)$coeff[4,3] #linear term
  
  rh_lm_exec_accuracy_t_uncor[i] <- summary(exec_accuracy_lm_model)$coeff[5,3]
  
}








#################################################################################
#################################################################################

#####################################################
#                     results                       #
#####################################################

#### FDR correction ####
for (hemi in hemis) {
  for (model in models) {
    hemi_model_p_unc <- paste0(hemi, "_gam_", model, "_p_uncor") 
    hemi_model_p_fdr <- paste0(hemi, "_gam_", model, "_p_fdr")    
    
    print(hemi_model_p_unc)
    
    #correct p values
    pfdr <- eval(substitute(p.adjust(i, method="fdr"), list(i = as.name(hemi_model_p_unc))))
    
    #figure out which values are < 0.05 and add to pfdr matrix
    pfdr <- as.data.frame(pfdr)
    pfdr$sig <- ifelse(pfdr<0.05, 1, 0)
    pfdr$sig_noNA <- ifelse(is.na(pfdr$sig), 0, pfdr$sig)
    names(pfdr) <- c("pfdr", "sig05", "sig05_noNA")
    hemi_model_p_fdr <- as.data.frame(pfdr[,1]) #sig05
    
    
    
    #multiply T values by fdr vector to get the list of Ts that are fdr corrected
    hemi_model_t_unc <- paste0(hemi, "_gam_", model, "_t_uncor")
    hemi_model_t_fdr <- paste0(hemi, "_gam_", model, "_t_fdr")
    t_df <- eval(substitute(as.data.frame(i), list(i = as.name(hemi_model_t_unc))))
    names(t_df) <- c("tval")
    t_df$tfdr <- pfdr[,3] * t_df$tval
    hemi_model_t_fdr <- as.data.frame(t_df[,2])
    
    
    #######################
    #### write tables #####
    #######################
    
    ## uncorrected ##
    
    ### p
    filename <- paste0(homedir, "/baller/results/coupling_accuracy/", hemi_model_p_unc, ".csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_p_unc, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ### t
    
    filename <- paste0(homedir, "/baller/results/coupling_accuracy/", hemi_model_t_unc, ".csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_t_unc, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ## corrected ##
    
    ### p
    filename <- paste0(homedir, "/baller/results/coupling_accuracy/", hemi, "_gam_", model, "_p_fdr05.csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_p_fdr, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ### t
    filename <- paste0(homedir, "/baller/results/coupling_accuracy/", hemi, "_gam_", model, "_t_fdr05.csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_t_fdr, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
  }
}

######  linear model alone #####

#### FDR correction ####
for (hemi in hemis) {
  for (model in lm_models) {
    hemi_model_p_unc <- paste0(hemi, "_lm_", model, "_p_uncor") 
    hemi_model_p_fdr <- paste0(hemi, "_lm_", model, "_p_fdr")    
    
    print(hemi_model_p_unc)
    
    #correct p values
    pfdr <- eval(substitute(p.adjust(i, method="fdr"), list(i = as.name(hemi_model_p_unc))))
    
    #figure out which values are < 0.05 and add to pfdr matrix
    pfdr <- as.data.frame(pfdr)
    pfdr$sig <- ifelse(pfdr<0.05, 1, 0)
    pfdr$sig_noNA <- ifelse(is.na(pfdr$sig), 0, pfdr$sig)
    names(pfdr) <- c("pfdr", "sig05", "sig05_noNA")
    hemi_model_p_fdr <- as.data.frame(pfdr[,1]) #sig05
    
    
    
    #multiply T values by fdr vector to get the list of Ts that are fdr corrected
    hemi_model_t_unc <- paste0(hemi, "_lm_", model, "_t_uncor")
    hemi_model_t_fdr <- paste0(hemi, "_lm_", model, "_t_fdr")
    t_df <- eval(substitute(as.data.frame(i), list(i = as.name(hemi_model_t_unc))))
    names(t_df) <- c("tval")
    t_df$tfdr <- pfdr[,3] * t_df$tval
    hemi_model_t_fdr <- as.data.frame(t_df[,2])
    
    
    #######################
    #### write tables #####
    #######################
    
    ## uncorrected ##
    
    ### p
    filename <- paste0(homedir, "/baller/results/coupling_accuracy/", hemi_model_p_unc, ".csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_p_unc, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ### t
    
    filename <- paste0(homedir, "/baller/results/coupling_accuracy/", hemi_model_t_unc, ".csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_t_unc, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ## corrected ##
    
    ### p
    filename <- paste0(homedir, "/baller/results/coupling_accuracy/", hemi, "_lm_", model, "_p_fdr05.csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_p_fdr, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ### t
    filename <- paste0(homedir, "/baller/results/coupling_accuracy/", hemi, "_lm_", model, "_t_fdr05.csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_t_fdr, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
  }
}

###########
## Plots ##
###########


### Find Gam age FDR corrected
### just for age
lh_gam_age_t_fdr05 <- read.csv(paste0(homedir, '/baller/results/coupling_accuracy/lh_gam_age_t_fdr05.csv'), header = F)
rh_gam_age_t_fdr05 <- read.csv(paste0(homedir, '/baller/results/coupling_accuracy/rh_gam_age_t_fdr05.csv'), header = F)

#turn 0s into na
lh_gam_to_keep <- lh_gam_age_t_fdr05$V1
lh_gam_to_keep[lh_gam_to_keep==0] <- NA

rh_gam_to_keep <- rh_gam_age_t_fdr05$V1
rh_gam_to_keep[rh_gam_to_keep==0] <- NA

#concatenate rows to keep
lh_and_rh_to_keep <- c(lh_gam_to_keep, rh_gam_to_keep)
columns_to_drop <- which(is.na(lh_and_rh_to_keep))

# concatenate 831x10242 right and left matrices, and drop the columns with NA
lh_and_rh_matrix <- data.frame(cbind(lh_matrix, rh_matrix))
lh_and_rh_matrix_fdr_corrected = subset(lh_and_rh_matrix, select = -(columns_to_drop)) 

lh_and_rh_rowmeans <- rowMeans(lh_and_rh_matrix_fdr_corrected) #make vector of row means
age = subjDemos$ageAtScan1
#age_x2_for_plotting <- c(subjDemos$ageAtScan1, subjDemos$ageAtScan1) #make vector of ages

#plot and save
jpeg(paste0(homedir,'baller/results/images/Mean_coupling_by_age_rplot_fdr.jpg'))
plot(age, lh_and_rh_rowmeans, ylab = "Coupling (Z)", xlab = "Age (In Years)", main = "Mean Coupling by Age")
abline(lm(lh_and_rh_rowmeans~age), col = 'red')
dev.off()

### Find Gam age uncorrected
### just for age
lh_gam_age_t_uncor <- read.csv(paste0(homedir, '/baller/results/coupling_accuracy/lh_gam_age_t_uncor.csv'), header = F)
rh_gam_age_t_uncor <- read.csv(paste0(homedir, '/baller/results/coupling_accuracy/rh_gam_age_t_uncor.csv'), header = F)

#turn 0s into na
lh_gam_to_keep <- lh_gam_age_t_uncor$V1
lh_gam_to_keep[lh_gam_to_keep==0] <- NA

rh_gam_to_keep <- rh_gam_age_t_uncor$V1
rh_gam_to_keep[rh_gam_to_keep==0] <- NA

#concatenate rows to keep
lh_and_rh_to_keep <- c(lh_gam_to_keep, rh_gam_to_keep)
columns_to_drop <- which(is.na(lh_and_rh_to_keep))

# concatenate 831x10242 right and left matrices, and drop the columns with NA
lh_and_rh_matrix <- data.frame(cbind(lh_matrix, rh_matrix))
lh_and_rh_matrix_uncorrected = subset(lh_and_rh_matrix, select = -(columns_to_drop)) 

lh_and_rh_rowmeans <- rowMeans(lh_and_rh_matrix_uncorrected) #make vector of row means
age = subjDemos$ageAtScan1
#age_x2_for_plotting <- c(subjDemos$ageAtScan1, subjDemos$ageAtScan1) #make vector of ages

#plot and save
jpeg(paste0(homedir,'baller/results/images/Mean_coupling_by_age_rplot_uncor.jpg'))
plot(age, lh_and_rh_rowmeans, ylab = "Coupling (Z)", xlab = "Age (In Years)", main = "Mean Coupling by Age")
abline(lm(lh_and_rh_rowmeans~age), col = 'red')
dev.off()

#visreg
jpeg(paste0(homedir,'baller/results/images/Mean_coupling_by_age_rplot_uncor_visreg.jpg'))
df <- data.frame(cbind(lh_and_rh_rowmeans, age))
names(df) <- c("Means", "Age")
fit <- lm(Means~Age, data=df)
visreg(fit, "Age")
dev.off()
