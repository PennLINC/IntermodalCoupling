###########
## Plots ##
###########


###pre: gam age results, both fdr and uncorrected. rh and lh matrices
###post: jpegs of linear and gam models, and derivatives 
###uses: makes the scatter plot of age by couping means for Figure 2
### dependencies: R 3.6.3

library(mgcv)
library(dplyr)
library(ggplot2)
library(visreg)
library(gratia)

#set home directory, switch this depending on whether running from PMACS or from home directory
#homedir <- "/Users/eballer/BBL/imco/pmacs/PMACS_remote/"
homedir <- "/project/imco/"
numrows = 831


# read in lh and rh matrices, which have demos appended
lh_df <- read.table(paste0(homedir, "/baller/results/coupling_accuracy/subjDemos_with_lh_", numrows, "x10242.csv"), sep = ",")
rh_df <- read.table(paste0(homedir, "/baller/results/coupling_accuracy/subjDemos_with_rh_", numrows, "x10242.csv"), sep = ",")

#just get the coupling values
lh_matrix <- lh_df[,(dim(lh_df)[2]-10241):dim(lh_df)[2]]
rh_matrix <- rh_df[,(dim(rh_df)[2]-10241):dim(rh_df)[2]]

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
age = lh_df$ageAtScan1

#plot and save
jpeg(paste0(homedir,'baller/results/images/Mean_coupling_by_age_rplot_fdr05.jpg'))
plot(age, lh_and_rh_rowmeans, ylab = "Coupling (Z)", xlab = "Age (In Years)", main = "Mean Coupling by Age")
abline(lm(lh_and_rh_rowmeans~age), col = 'red')
dev.off()

#visreg

jpeg(paste0(homedir,'baller/results/images/Mean_coupling_by_age_rplot_fdr05_visreg.jpg'))
fit <- lm(Means~Age, data=df)
visreg(fit, "Age", ylab = "Mean Coupling (Z)", xlab = "Age (in Years)",
       line=list(col="red"),
       points=list(col="black"))
dev.off()

jpeg(paste0(homedir, 'baller/results/images/Mean_coupling_by_age_rplot_fdr05_visreg_gam.jpg'))
fit_gam <- gam(Means~s(Age, k = 4, fx=T), data=df)
visreg(fit_gam, "Age", ylab = "Mean Coupling (Z)", xlab = "Age (in Years)",
       line=list(col="red"),
       points=list(col="black"))
dev.off()


### SNR 50
### Find Gam age FDR corrected
### just for age
lh_gam_age_t_fdr05 <- read.csv(paste0(homedir, '/baller/results/coupling_accuracy/lh_gam_age_t_fdr05.csv'), header = F)
rh_gam_age_t_fdr05 <- read.csv(paste0(homedir, '/baller/results/coupling_accuracy/rh_gam_age_t_fdr05.csv'), header = F)


#snr
lh_snr <- read.csv(paste0(homedir, '/baller/results/CR_revision/lh_snr_map_med_wall_-1.csv'), header = F)
rh_snr <- read.csv(paste0(homedir, '/baller/results/CR_revision/rh_snr_map_med_wall_-1.csv'), header = F)

#turn 0s into na
lh_gam_to_keep <- lh_gam_age_t_fdr05$V1
lh_gam_to_keep[lh_gam_to_keep==0] <- NA

rh_gam_to_keep <- rh_gam_age_t_fdr05$V1
rh_gam_to_keep[rh_gam_to_keep==0] <- NA

#put NAs in also if snr50 NA
lh_gam_to_keep[lh_snr$V1<50] <- NA
rh_gam_to_keep[rh_snr$V1<50] <- NA

#concatenate rows to keep
lh_and_rh_to_keep <- c(lh_gam_to_keep, rh_gam_to_keep)
columns_to_drop <- which(is.na(lh_and_rh_to_keep))

# concatenate 831x10242 right and left matrices, and drop the columns with NA
lh_and_rh_matrix <- data.frame(cbind(lh_matrix, rh_matrix))
lh_and_rh_matrix_fdr_corrected = subset(lh_and_rh_matrix, select = -(columns_to_drop)) 

lh_and_rh_rowmeans <- rowMeans(lh_and_rh_matrix_fdr_corrected) #make vector of row means
age = lh_df$ageAtScan1

#plot and save
jpeg(paste0(homedir,'baller/results/CR_revision/images/Mean_coupling_by_age_rplot_fdr05_snr50.jpg'))
plot(age, lh_and_rh_rowmeans, ylab = "Coupling (Z)", xlab = "Age (In Years)", main = "Mean Coupling by Age")
abline(lm(lh_and_rh_rowmeans~age), col = 'red')
dev.off()

#visreg

jpeg(paste0(homedir,'baller/results/CR_revision/images/Mean_coupling_by_age_rplot_fdr05_visreg_snr50.jpg'))
fit <- lm(Means~Age, data=df)
visreg(fit, "Age", ylab = "Mean Coupling (Z)", xlab = "Age (in Years)",
       line=list(col="red"),
       points=list(col="black"))
dev.off()

jpeg(paste0(homedir, 'baller/results/CR_revision/images/Mean_coupling_by_age_rplot_fdr05_visreg_gam_snr50.jpg'))
fit_gam <- gam(Means~s(Age, k = 4, fx=T), data=df)
visreg(fit_gam, "Age", ylab = "Mean Coupling (Z)", xlab = "Age (in Years)",
       line=list(col="red"),
       points=list(col="black"))
dev.off()

#################
## Derivatives ##
#################
d<-derivatives(fit_gam,n=1000)
d_plot <- draw(d)
print(d_plot)
ggsave(plot = d_plot,filename = paste0(homedir, "baller/results/CR_revision/images/derivative_plot_age_gam_fdr_corrected_snr50.png"),device = "png",width = 180,height = 120,units = "mm")

d<- d %>%
  mutate(sig = !(0 >lower & 0 < upper)) #Ages where the CI does not include zero
cat(sprintf("\nSignificant change: %1.2f - %1.2f\n",min(d$data[d$sig==T]),max(d$data[d$sig==T]))) #this only work




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
age = lh_df$ageAtScan1
#age_x2_for_plotting <- c(subjDemos$ageAtScan1, subjDemos$ageAtScan1) #make vector of ages

#plot and save
jpeg(paste0(homedir,'baller/results/images/Mean_coupling_by_age_rplot_uncor.jpg'))
plot(age, lh_and_rh_rowmeans, ylab = "Coupling (Z)", xlab = "Age (In Years)", main = "Mean Coupling by Age")
abline(lm(lh_and_rh_rowmeans~age), col = 'red')
dev.off()

#visreg
df <- data.frame(cbind(lh_and_rh_rowmeans, age))
names(df) <- c("Means", "Age")

jpeg(paste0(homedir,'baller/results/images/Mean_coupling_by_age_rplot_uncor_visreg.jpg'))
fit <- lm(Means~Age, data=df)
visreg(fit, "Age", ylab = "Mean Coupling (Z)", xlab = "Age (in Years)",
       line=list(col="red"),
       points=list(col="black"))
dev.off()

jpeg(paste0(homedir, 'baller/results/images/Mean_coupling_by_age_rplot_uncor_visreg_gam.jpg'))
fit_gam <- gam(Means~s(Age, k = 4, fx=T), data=df)
visreg(fit_gam, "Age", ylab = "Mean Coupling (Z)", xlab = "Age (in Years)",
       line=list(col="red"),
       points=list(col="black"))
dev.off()


######### SNR 50
### Find Gam age uncorrected
### just for age
lh_gam_age_t_uncor <- read.csv(paste0(homedir, '/baller/results/coupling_accuracy/lh_gam_age_t_uncor.csv'), header = F)
rh_gam_age_t_uncor <- read.csv(paste0(homedir, '/baller/results/coupling_accuracy/rh_gam_age_t_uncor.csv'), header = F)

lh_snr <- read.csv(paste0(homedir, '/baller/results/CR_revision/lh_snr_map_med_wall_-1.csv'), header = F)
rh_snr <- read.csv(paste0(homedir, '/baller/results/CR_revision/rh_snr_map_med_wall_-1.csv'), header = F)

#turn 0s into na
lh_gam_to_keep <- lh_gam_age_t_uncor$V1
lh_gam_to_keep[lh_gam_to_keep==0] <- NA

rh_gam_to_keep <- rh_gam_age_t_uncor$V1
rh_gam_to_keep[rh_gam_to_keep==0] <- NA

#put NAs in also if snr50 NA
lh_gam_to_keep[lh_snr$V1<50] <- NA
rh_gam_to_keep[rh_snr$V1<50] <- NA

#concatenate rows to keep
lh_and_rh_to_keep <- c(lh_gam_to_keep, rh_gam_to_keep)
columns_to_drop <- which(is.na(lh_and_rh_to_keep))

# concatenate 831x10242 right and left matrices, and drop the columns with NA
lh_and_rh_matrix <- data.frame(cbind(lh_matrix, rh_matrix))
lh_and_rh_matrix_uncorrected = subset(lh_and_rh_matrix, select = -(columns_to_drop)) 

lh_and_rh_rowmeans <- rowMeans(lh_and_rh_matrix_uncorrected) #make vector of row means
age = lh_df$ageAtScan1
#age_x2_for_plotting <- c(subjDemos$ageAtScan1, subjDemos$ageAtScan1) #make vector of ages

#plot and save
jpeg(paste0(homedir,'baller/results/CR_revision/images/Mean_coupling_by_age_rplot_uncor_snr50.jpg'))
plot(age, lh_and_rh_rowmeans, ylab = "Coupling (Z)", xlab = "Age (In Years)", main = "Mean Coupling by Age")
abline(lm(lh_and_rh_rowmeans~age), col = 'red')
dev.off()

#visreg
df <- data.frame(cbind(lh_and_rh_rowmeans, age))
names(df) <- c("Means", "Age")

jpeg(paste0(homedir,'baller/results/CR_revision/images/Mean_coupling_by_age_rplot_uncor_visreg_snr50.jpg'))
fit <- lm(Means~Age, data=df)
visreg(fit, "Age", ylab = "Mean Coupling (Z)", xlab = "Age (in Years)",
       line=list(col="red"),
       points=list(col="black"))
dev.off()

write.csv(df, paste0(homedir, '/baller/results/CR_revision/images/data_frame_to_do_deriv_on_locally_bc_no_gratia_pmacs.csv'), row.names = F)
jpeg(paste0(homedir, 'baller/results/CR_revision/images/Mean_coupling_by_age_rplot_uncor_visreg_gam_snr50.jpg'))
fit_gam <- gam(Means~s(Age, k = 4, fx=T), data=df)
visreg(fit_gam, "Age", ylab = "Mean Coupling (Z)", xlab = "Age (in Years)",
       line=list(col="red"),
       points=list(col="black"))
dev.off()