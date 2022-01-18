###############################################
### PUT NA In SNR MAP  ###
###############################################

### Author: Erica Baller
### Date 12/22/21

### pre: statistical snr map (l and r), each 10242 vertices
### post: l and r maps, with NAs in medial wall
### uses: Need snr maps to have NA in medial wall for spin, which I'll then use to calculate correlations rather than network enrichment 
### dependencis: and R will do. I used 3.2.5

## set abs and relative paths. Must toggle before running locally/on cluster
#homedir = "/Users/eballer/BBL/imco/pmacs/PMACS_remote"
homedir = "/project/imco"

source(paste0(homedir, '/baller/scripts/imco_functions.R'))

#### read in l and r snr maps
lh_snr_map <- read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/tSNR_Avg_lh.csv"), header = F) 
rh_snr_map <- read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/tSNR_Avg_rh.csv"), header = F)


#### get medial wall indices and change snr maps to have NA in those spots
parcel_mapping <- get_parcel_mapping_yeo(7)
lh_numerical_map <- parcel_mapping[[1]]
rh_numerical_map <- parcel_mapping[[2]]

lh_medial_wall_nums <- which(lh_numerical_map == 8)
rh_medial_wall_nums <- which(rh_numerical_map == 8)
lh_snr_map$V1[lh_medial_wall_nums] = NA
rh_snr_map$V1[rh_medial_wall_nums] = NA

#### save new maps
#save maps
write.table(x = lh_snr_map$V1, file = (paste0(homedir, '/baller/results/CR_revision/lh_snr_map_med_wall_NA.csv')), quote = F, row.names = F, col.names = F)
write.table(x = lh_snr_map$V1, file = (paste0(homedir, '/baller/results/CR_revision/rh_snr_map_med_wall_NA.csv')), quote = F, row.names = F, col.names = F)
