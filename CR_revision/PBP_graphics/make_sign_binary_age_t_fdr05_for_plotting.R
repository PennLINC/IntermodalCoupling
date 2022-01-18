#need to use the sign of the lm for this, swap for visualization purposes
lh_lm <- read.csv(file = "/project/imco/baller/results/CR_revision/coupling_accuracy/lh_lm_age_t_fdr05.csv", header = F)
lh_lm$binary=-1
lh_lm$binary[lh_lm$V1>0] <- 1
rh_lm <- read.csv(file = "/project/imco/baller/results/CR_revision/coupling_accuracy/rh_lm_age_t_fdr05.csv", header = F)
rh_lm$binary = -1
rh_lm$binary[rh_lm$V1>0] <- 1
write.table(lh_lm$binary, file = "/project/imco/baller/results/CR_revision/coupling_accuracy/lh_lm_sign_binary_age_t_fdr05.csv", row.names = F, col.names = F)
write.table(rh_lm$binary, file = "/project/imco/baller/results/CR_revision/coupling_accuracy/rh_lm_sign_binary_age_t_fdr05.csv", row.names = F, col.names = F)

###do the actual multiplication

lh_agexmask_10 <- read.csv("/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_10_lh.csv", header = F)
rh_agexmask_10 <- read.csv("/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_10_rh.csv", header = F)

lh_agexmask_25 <- read.csv("/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_25_lh.csv", header = F)
rh_agexmask_25 <- read.csv("/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_25_rh.csv", header = F)

lh_agexmask_35 <- read.csv("/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_35_lh.csv", header = F)
rh_agexmask_35 <- read.csv("/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_35_rh.csv", header = F)

lh_agexmask_50 <- read.csv("/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_50_lh.csv", header = F)
rh_agexmask_50 <- read.csv("/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_50_rh.csv", header = F)


#multiply to get sign

lh_agexmask_10$signed <- lh_agexmask_10$V1 * lh_lm$binary
rh_agexmask_10$signed <- rh_agexmask_10$V1 * rh_lm$binary

lh_agexmask_25$signed <- lh_agexmask_25$V1 * lh_lm$binary
rh_agexmask_25$signed <- rh_agexmask_25$V1 * rh_lm$binary

lh_agexmask_35$signed <- lh_agexmask_35$V1 * lh_lm$binary
rh_agexmask_35$signed <- rh_agexmask_35$V1 * rh_lm$binary

lh_agexmask_50$signed <- lh_agexmask_50$V1 * lh_lm$binary
rh_agexmask_50$signed <- rh_agexmask_50$V1 * rh_lm$binary

### write output

write.table(lh_agexmask_10$signed, file = "/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_10_lh_lm_signed.csv", row.names = F, col.names = F)
write.table(rh_agexmask_10$signed, file = "/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_10_rh_lm_signed.csv", row.names = F, col.names = F)

write.table(lh_agexmask_25$signed, file = "/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_25_lh_lm_signed.csv", row.names = F, col.names = F)
write.table(rh_agexmask_25$signed, file = "/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_25_rh_lm_signed.csv", row.names = F, col.names = F)

write.table(lh_agexmask_35$signed, file = "/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_35_lh_lm_signed.csv", row.names = F, col.names = F)
write.table(rh_agexmask_35$signed, file = "/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_35_rh_lm_signed.csv", row.names = F, col.names = F)

write.table(lh_agexmask_50$signed, file = "/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_50_lh_lm_signed.csv", row.names = F, col.names = F)
write.table(rh_agexmask_50$signed, file = "/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_50_rh_lm_signed.csv", row.names = F, col.names = F)




