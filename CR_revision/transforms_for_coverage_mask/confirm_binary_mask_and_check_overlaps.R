lh_cbfxalff_mask <- read.table("/project/imco/baller/transforms_for_surface_projection/lh_cbfxalff_mask_binary", header = F)
rh_cbfxalff_mask <- read.table("/project/imco/baller/transforms_for_surface_projection/rh_cbfxalff_mask_binary", header = F)

lh_cbfxalff_mask$binary <- ifelse(lh_cbfxalff_mask$V1 == 0, 0, 1) #n = 9495 in mask
rh_cbfxalff_mask$binary <- ifelse(rh_cbfxalff_mask$V1 == 0, 0, 1) #n = 9505 in mask

write.table(lh_cbfxalff_mask$binary, file = "/project/imco/baller/transforms_for_surface_projection/lh_cbfxalff_mask_binary_R.csv", row.names = F, col.names = F, quote = F)
write.table(rh_cbfxalff_mask$binary, file = "/project/imco/baller/transforms_for_surface_projection/rh_cbfxalff_mask_binary_R.csv", row.names = F, col.names = F, quote = F)

snr_lh <- read.table("/project/imco/baller/results/CR_revision/couplingxsnr_maps/snrxmask_50_lh.csv", header = F) #n=8324 that are not 0 or NA
snr_rh <- read.table("/project/imco/baller/results/CR_revision/couplingxsnr_maps/snrxmask_50_rh.csv", header = F) #n = 8297 that are not 0 or NA

#check for overlap between areas ruled out for crappy SNR and the poorly sampled mask
snr_lh$cbfxalff_masked <- snr_lh$V1 * lh_cbfxalff_mask$binary #n = 8245 
snr_rh$cbfxalff_masked <- snr_rh$V1 * rh_cbfxalff_mask$binary #n = 8201

mean_coupling_lh <- read.table("/project/imco/baller/results/CR_revision/couplingxsnr_maps/mean_couplingxmask_50_lh.csv", header = F) #n = 8240
mean_coupling_rh <- read.table("/project/imco/baller/results/CR_revision/couplingxsnr_maps/mean_couplingxmask_50_rh.csv", header = F) #n = 8209

mean_coupling_lh$cbfxalff_masked <- mean_coupling_lh$V1*lh_cbfxalff_mask$binary #n = 8161 (n=79 vertices less)
mean_coupling_rh$cbfxalff_masked <- mean_coupling_rh$V1*rh_cbfxalff_mask$binary #n = 8113 (96 vertices less)

write.table(mean_coupling_lh$cbfxalff_masked, file="/project/imco/baller/transforms_for_surface_projection/lh_mean_couplingxcbfxalff_mask_binary_snr50.csv", row.names = F, col.names = F, quote = F)
write.table(mean_coupling_rh$cbfxalff_masked, file="/project/imco/baller/transforms_for_surface_projection/rh_mean_couplingxcbfxalff_mask_binary_snr50.csv", row.names = F, col.names = F, quote = F)

age_lh <- read.table("/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_50_lh_lm_signed.csv", header = F) #n = 6327 not NA or 0
age_rh <- read.table("/project/imco/baller/results/CR_revision/couplingxsnr_maps/agexmask_50_rh_lm_signed.csv", header = F) #n = 6203

age_lh$cbfalff_masked <- age_lh$V1*lh_cbfxalff_mask$binary #n = 6320 (7 vertices less)
age_rh$cbfalff_masked <- age_rh$V1*rh_cbfxalff_mask$binary #n=6199 (4 vertices less)

sex_lh <- read.table("/project/imco/baller/results/CR_revision/couplingxsnr_maps/sexxmask_50_lh.csv", header = F) #n= 1101
sex_rh <- read.table("/project/imco/baller/results/CR_revision/couplingxsnr_maps/sexxmask_50_rh.csv", header = F) #n = 1764

sex_lh$cbfalff_masked <- sex_lh$V1*lh_cbfxalff_mask$binary #n=1067 (34 vertices less)
sex_rh$cbfalff_masked <- sex_rh$V1*rh_cbfxalff_mask$binary #n=1733 (31 vertices less)

ea_lh <- read.table("/project/imco/baller/results/CR_revision/couplingxsnr_maps/eaxmask_50_lh.csv", header = F) #n = 340
ea_rh <- read.table("/project/imco/baller/results/CR_revision/couplingxsnr_maps/eaxmask_50_rh.csv", header = F) #n = 348

ea_lh$cbfalff_masked <- ea_lh$V1*lh_cbfxalff_mask$binary #n = 340 (no change)
ea_rh$cbfalff_masked <- ea_rh$V1*rh_cbfxalff_mask$binary #n = 348 (no change)
