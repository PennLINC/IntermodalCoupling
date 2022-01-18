###############################################
### PUT -1 In SNR MAP  ###
###############################################

### Author: Erica Baller
### Date 12/22/21

### pre: statistical snr map (l and r), each 10242 vertices
### post: l and r maps, with -1s in medial wall
### uses: Need snr maps to have -1 in medial wall for spin, which I'll then use to calculate correlations rather than network enrichment 
### dependencis: and R will do. I used 3.2.5

## set abs and relative paths. Must toggle before running locally/on cluster
#homedir = "/Users/eballer/BBL/imco/pmacs/PMACS_remote"
homedir = "/project/imco"

source(paste0(homedir, '/baller/scripts/imco_functions.R'))

#### read in l and r snr maps
lh_snr_map <- read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/tSNR_Avg_lh.csv"), header = F) 
rh_snr_map <- read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/tSNR_Avg_rh.csv"), header = F)


#### get medial wall indices and change snr maps to have -1 in those spots
parcel_mapping <- get_parcel_mapping_yeo(7)
lh_numerical_map <- parcel_mapping[[1]]
rh_numerical_map <- parcel_mapping[[2]]

lh_medial_wall_nums <- which(lh_numerical_map == 8)
rh_medial_wall_nums <- which(rh_numerical_map == 8)
lh_snr_map$V1[lh_medial_wall_nums] = -1
rh_snr_map$V1[rh_medial_wall_nums] = -1

#set the mask to treat areas not included as NA, rather than 0, because 0 is a value for SNR
lh_snr_map$mask_25 <- ifelse(lh_snr_map$V1 <25, -1, lh_snr_map$V1)
rh_snr_map$mask_25 <- ifelse(rh_snr_map$V1 <25, -1, rh_snr_map$V1)
lh_snr_map$mask_35 <- ifelse(lh_snr_map$V1 <35, -1, lh_snr_map$V1)
rh_snr_map$mask_35 <- ifelse(rh_snr_map$V1 <35, -1, rh_snr_map$V1)
lh_snr_map$mask_40 <- ifelse(lh_snr_map$V1 <40, -1, lh_snr_map$V1)
rh_snr_map$mask_40 <- ifelse(rh_snr_map$V1 <40, -1, rh_snr_map$V1)
lh_snr_map$mask_45 <- ifelse(lh_snr_map$V1 <45, -1, lh_snr_map$V1)
rh_snr_map$mask_45 <- ifelse(rh_snr_map$V1 <45, -1, rh_snr_map$V1)
lh_snr_map$mask_50 <- ifelse(lh_snr_map$V1 <50, -1, lh_snr_map$V1)
rh_snr_map$mask_50 <- ifelse(rh_snr_map$V1 <50, -1, rh_snr_map$V1)

#### save new maps
#save maps
write.table(x = lh_snr_map$V1, file = (paste0(homedir, '/baller/results/CR_revision/lh_snr_map_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
write.table(x = rh_snr_map$V1, file = (paste0(homedir, '/baller/results/CR_revision/rh_snr_map_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
write.table(x = lh_snr_map$mask_25, file = (paste0(homedir, '/baller/results/CR_revision/lh_snr_25_map_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
write.table(x = rh_snr_map$mask_25, file = (paste0(homedir, '/baller/results/CR_revision/rh_snr_25_map_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
write.table(x = lh_snr_map$mask_35, file = (paste0(homedir, '/baller/results/CR_revision/lh_snr_35_map_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
write.table(x = rh_snr_map$mask_35, file = (paste0(homedir, '/baller/results/CR_revision/rh_snr_35_map_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
write.table(x = lh_snr_map$mask_40, file = (paste0(homedir, '/baller/results/CR_revision/lh_snr_40_map_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
write.table(x = rh_snr_map$mask_40, file = (paste0(homedir, '/baller/results/CR_revision/rh_snr_40_map_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
write.table(x = lh_snr_map$mask_45, file = (paste0(homedir, '/baller/results/CR_revision/lh_snr_45_map_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
write.table(x = rh_snr_map$mask_45, file = (paste0(homedir, '/baller/results/CR_revision/rh_snr_45_map_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
write.table(x = lh_snr_map$mask_50, file = (paste0(homedir, '/baller/results/CR_revision/lh_snr_50_map_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
write.table(x = rh_snr_map$mask_50, file = (paste0(homedir, '/baller/results/CR_revision/rh_snr_50_map_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
