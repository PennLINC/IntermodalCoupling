###########################################
### Make Trinarized Maps for Spin Test  ###
###########################################

### Author: Erica Baller
### Date 3/10/21

### pre: statistical maps (l and r), each 10242 vertices
### post: 6 matrices with 1s if value fdr corrected, 0 if not, and -1 if vertex in medial wall of yeo map
# lh and rh pos and neg together, positive alone, and negative alone
### uses: In order to spin my vertices in a way that allows me to later assess their yeo membership, I need all of my fdr values indicated, as well as which values are in the medial wall.  This script takes lh and rh fdr corrected vectors generated in previous scripts, booleanizes them into in/out of fdr map. Then, -1s are added for vertices in the medial wall. 
### dependencis: and R will do. I used 3.2.5

## set abs and relative paths. Must toggle before running locally/on cluster
#homedir = "/Users/eballer/BBL/imco/pmacs/PMACS_remote"
homedir = "/project/imco"

source(paste0(homedir, '/baller/scripts/imco_functions.R'))

#### set # parcels in case I want to do 7 or 17 or something else in the future
parcel_type = "Yeo" 
parcel_num = 7 
input_parcel_array_length = 10242

#if I'd like to make mean maps as well
mean_maps = 1 
#will loop through each of these
analyses <- c("coupling_accuracy") #, "cbf_accuracy", "alff_accuracy") 
models <- c("gam_age", "gam_sex", "gam_exec_accuracy") #gam_age", "lm_age", "gam_sex", "lm_sex", "lm_accuracy", "lm_exec_accuracy")
corrs <- c("fdr05") #correction

parcel_mapping <- get_parcel_mapping_yeo(7)
lh_numerical_map <- parcel_mapping[[1]]
rh_numerical_map <- parcel_mapping[[2]]

for (analysis in analyses) {
  for (model in models) {
    for (corr in corrs) {  
      
      ### set results path
      stat_path <- paste0("/", analysis, "/")
      print(stat_path)
      result_path <- paste0(model, "_t_", corr)
      print(result_path)
      
      ### set paths
      ## input
      lh_stat_map <- read.csv(paste0(homedir, "/baller/results/", stat_path, "lh_", result_path, ".csv"), header = F)
      rh_stat_map <- read.csv(paste0(homedir, "/baller/results/", stat_path, "rh_", result_path, ".csv"), header = F)
      
      ## output
      lh_outdir <- paste0(homedir, "/baller/results/", stat_path, "lh_", result_path, "_", parcel_type, parcel_num, "_1_0_-1.csv")
      rh_outdir <- paste0(homedir, "/baller/results/", stat_path, "rh_", result_path, "_", parcel_type, parcel_num, "_1_0_-1.csv")
      
      lh_outdir_pos <- paste0(homedir, "/baller/results/", stat_path, "lh_pos_", result_path, "_", parcel_type, parcel_num, "_1_0_-1.csv")
      rh_outdir_pos <- paste0(homedir, "/baller/results/", stat_path, "rh_pos_", result_path, "_", parcel_type, parcel_num, "_1_0_-1.csv")
      
      lh_outdir_neg <- paste0(homedir, "/baller/results/", stat_path, "lh_neg_", result_path, "_", parcel_type, parcel_num, "_1_0_-1.csv")
      rh_outdir_neg <- paste0(homedir, "/baller/results/", stat_path, "rh_neg_", result_path, "_", parcel_type, parcel_num, "_1_0_-1.csv")
      
      #convert stat map to boolean
      lh_stat_boolean <- ifelse(lh_stat_map == 0, 0, 1)
      rh_stat_boolean <- ifelse(rh_stat_map == 0, 0, 1)
      
      lh_stat_boolean_pos <- ifelse(lh_stat_map > 0, 1, 0)
      rh_stat_boolean_pos <- ifelse(rh_stat_map > 0, 1, 0)
      
      lh_stat_boolean_neg <- ifelse(lh_stat_map < 0, 1, 0)
      rh_stat_boolean_neg <- ifelse(rh_stat_map < 0, 1, 0)
      
      #convert NAs to 0
      lh_stat_boolean[is.na(lh_stat_boolean)] <- 0
      rh_stat_boolean[is.na(rh_stat_boolean)] <- 0
      
      lh_stat_boolean_pos[is.na(lh_stat_boolean_pos)] <- 0
      rh_stat_boolean_pos[is.na(rh_stat_boolean_pos)] <- 0
      
      lh_stat_boolean_neg[is.na(lh_stat_boolean_neg)] <- 0
      rh_stat_boolean_neg[is.na(rh_stat_boolean_neg)] <- 0
      
      
      #substitute -1 for those locations with medial wall stuff
      lh_medial_wall_nums <- which(lh_numerical_map == 8)
      rh_medial_wall_nums <- which(rh_numerical_map == 8)
      
      lh_stat_boolean[lh_medial_wall_nums] <- -1
      rh_stat_boolean[rh_medial_wall_nums] <- -1
      
      lh_stat_boolean_pos[lh_medial_wall_nums] <- -1
      rh_stat_boolean_pos[rh_medial_wall_nums] <- -1
      
      lh_stat_boolean_neg[lh_medial_wall_nums] <- -1
      rh_stat_boolean_neg[rh_medial_wall_nums] <- -1
      
      
      #write output
      write.table(x = lh_stat_boolean, file = lh_outdir, quote = F, row.names = F, col.names = F)
      write.table(x = rh_stat_boolean, file = rh_outdir, quote = F, row.names = F, col.names = F)
      
      write.table(x = lh_stat_boolean_pos, file = lh_outdir_pos, quote = F, row.names = F, col.names = F)
      write.table(x = rh_stat_boolean_pos, file = rh_outdir_pos, quote = F, row.names = F, col.names = F)
      
      write.table(x = lh_stat_boolean_neg, file = lh_outdir_neg, quote = F, row.names = F, col.names = F)
      write.table(x = rh_stat_boolean_neg, file = rh_outdir_neg, quote = F, row.names = F, col.names = F)
      
    }
  }
}
###########################
##### for mean maps #######
###########################

if (mean_maps == 1){
  #read in maps
  
  lh_map <- read.csv(paste0(homedir, '/baller/results/mean_maps/n831_lh.coupling_coef_alff_cbf.fwhm15.fsaverage5.csv'), header = F)
  rh_map <- read.csv(paste0(homedir, '/baller/results/mean_maps/n831_rh.coupling_coef_alff_cbf.fwhm15.fsaverage5.csv'), header = F)
  #set medial wall to -1
  lh_medial_wall_nums <- which(lh_numerical_map == 8)
  rh_medial_wall_nums <- which(rh_numerical_map == 8)
  lh_map$V1[lh_medial_wall_nums] = -1
  rh_map$V1[rh_medial_wall_nums] = -1
 
  #save maps
  write.table(x = lh_map$V1, file = (paste0(homedir, '/baller/results/lh_mean_coupling_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
  write.table(x = rh_map$V1, file = (paste0(homedir, '/baller/results/rh_mean_coupling_med_wall_-1.csv')), quote = F, row.names = F, col.names = F)
}
