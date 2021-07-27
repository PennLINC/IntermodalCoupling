require('visreg')
require('mgcv')
require('tableone')
require('dplyr')

make_demographics_table<- function(data_frame) {
  
  #subset demographics
  listVars <- c("Age", "Sex", "Race", "Maternal Ed") #Race 1 = caucasian, Maternal Ed = years, age = years
  demo <- data.frame(data_frame$ageAtScan1, data_frame$sex, data_frame$race2, data_frame$medu1)
  names(demo) <- c(listVars)
  
  #Change categorical values to have names
  demo$Race <- ifelse(demo$Race == 1, "Caucasian", "Non-caucasian")
  demo$Sex <- ifelse(demo$Sex == 1, "Male", "Female")
 
  #Define Categorical Variables
  cat_variables <- c("Sex", "Race")
  title <- paste0("IMCO Demographics, n = ", dim(demo)[1])
  
  #create demographics table
  demo_table <- CreateTableOne(vars = listVars, data = demo, factorVars = cat_variables)
  print(demo_table, showAllLevels = TRUE)
}

get_parcel_mapping_yeo <- function(parcel_num){
  
  #pre: input parcel #, either 7 or 17
  #post: list lh and rh yeo networks that map onto code #s
  #uses: easy way to translate the weird numerical maps in fsaverage 5 space into something we are more familiar with
  #dependencies: Any R will do, I used 3.2.5
  
  ## Set Yeo info
  #### set # parcels in case I want to do 7 or 17 or something else in the future
  parcel_type = "Yeo" 
  parcel_num = parcel_num 
  input_parcel_array_length = 10242
  
  # read in yeo fsaverage5 vectors
  parcelID <- read.csv(paste0(homedir, "/baller/processed_data/yeo_network_data/NetworkIDnumbers", parcel_type, parcel_num, ".csv"), header = F)
  parcelName <- t(read.csv(paste0(homedir, "/baller/processed_data/yeo_network_data/NetworkNames", parcel_type, parcel_num, ".csv"), header = F))
  
  lh_parcel_nums <- read.csv(paste0(homedir, "/baller/processed_data/yeo_network_data/lh_", input_parcel_array_length, "_vertex_nums_", parcel_type, parcel_num, ".csv"), header = F)
  rh_parcel_nums <- read.csv(paste0(homedir, "/baller/processed_data/yeo_network_data/rh_", input_parcel_array_length, "_vertex_nums_", parcel_type, parcel_num, ".csv"), header = F)
  
  # map Yeo numbers to parcels
  #make a column of numbers for mapping
  parcelID$network_num <- c(1:dim(parcelID)[1])
  
  #add extra row to parcelID, not clear why this didn't come from Yeo labels, maybe cerebellum?... 8 will equal 65793
  # comment this out if not using yeo 
  parcelID<- rbind(parcelID, c(65793, 8))
  
  #make vector for lh and rh with mapping
  lh_numerical_map <- lh_parcel_nums
  rh_numerical_map <- rh_parcel_nums
  
  #foreach vertex, which contains a bunch of numbers, match it to the appropriate column, and take the network num (i.e. yeo 2, which would correspond to Motor), associated with it
  lh_numerical_map[] <- lapply(lh_parcel_nums, function(x) parcelID$network_num[match(x, parcelID$V1)])
  rh_numerical_map[] <- lapply(rh_parcel_nums, function(x) parcelID$network_num[match(x, parcelID$V1)])

  lh_and_rh_numerical_map_list <- list(lh_numerical_map$V1, rh_numerical_map$V1)
  return(lh_and_rh_numerical_map_list)
  
}

melt_df_for_violin_plot_yeo7 <- function(df, network_names, num_spins){
  #melt dataframe so it is in a good format for violin plotting.
  #melt df so it is in a good position to be plotted
  melted_df_network_name <- rep(x = network_names, each = num_spins)
  melted_df_network_num <- rep(x = seq(1:7), each = num_spins)
  melted_df_spin_results <- rbind(t(df[1,]), 
                                  t(df[2,]),
                                  t(df[3,]),
                                  t(df[4,]),
                                  t(df[5,]),
                                  t(df[6,]),
                                  t(df[7,]))
  
  melted_df <- as.data.frame(cbind(melted_df_network_name, melted_df_network_num,melted_df_spin_results))
  names(melted_df) <- c("network_name", "network_num", "spin")
  melted_df$spin <- as.numeric(as.character(melted_df$spin))
  return(melted_df)
}

violin_plot_pos_and_neg_lines <- function (homedir, models, network_names, num_spins){
  for (model in models) {

    print(model)
    #for storing statistics at the end
    lh_spin_df <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/lh_spin_test_", model, "_proportions.csv"), sep = ",")
    )
    lh_spin_pos <-data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/lh_spin_test_pos_", model, "_proportions.csv"), sep = ",")
    )
    lh_spin_neg <-data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/lh_spin_test_neg_", model, "_proportions.csv"), sep = ",")
    )
    rh_spin_df <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/rh_spin_test_", model, "_proportions.csv"), sep = ",")
    )
    rh_spin_pos <-data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/rh_spin_test_pos_", model, "_proportions.csv"), sep = ",")
    )
    rh_spin_neg <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/rh_spin_test_neg_", model, "_proportions.csv"), sep = ",")
    )
    
    #take means of left and right
    actual_results <- (lh_spin_df[,1] + rh_spin_df[,1])/2
    actual_results_pos <- ((lh_spin_pos[,1] + rh_spin_pos[,1])/2)
    actual_results_neg <- ((lh_spin_neg[,1] + rh_spin_neg[,1])/2)
    
   # print(actual_results)
  #  print(actual_results_pos)
  #  print(actual_results_neg)
    #dataframes for all spins, as well as positive and negative
    spin_without_target_col <- cbind(lh_spin_df[,2:1001],rh_spin_df[,2:1001])
   # print(min(spin_without_target_col))
  #  print(spin_without_target_col[4,5:10])
    
    melted_df <- melt_df_for_violin_plot_yeo7(spin_without_target_col, network_names, num_spins)
  
  #  print(melted_df[1:4,2:3])
    #with mean lines, different fonts
    #save images
    plot_violin <- ggplot(melted_df, aes(x = factor(network_name, level = network_names), y = spin, fill = network_name)) +  
      scale_fill_manual(values=yeo_colors) + 
      geom_violin(trim = TRUE) + 
      xlab("Yeo 7 Network") + ylab(paste0("Proportion")) +
      ylim(0,NA) + 
      geom_violin(trim=FALSE) + 
      theme_classic() + 
      theme(legend.position = "none",
            legend.title = element_blank(),
            axis.text.x = element_text(size = 10, colour = "black"),
            axis.text.y = element_text(size = 10, colour = "black"),
            axis.title.y = element_text(size = 10),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 10)) +
      stat_summary(fun.y = mean, geom = "errorbar", 
                  aes(ymax = ..y.., ymin = ..y.., group = factor(network_name)),
                  width = 0.5, linetype = "dashed", position = position_dodge(0.9)) + 
      geom_segment(aes(x = 0.5, y = actual_results_pos[1], xend = 1.5, yend = actual_results_pos[1]), color="red") + 
      geom_segment(aes(x = 1.5, y = actual_results_pos[2], xend = 2.5, yend = actual_results_pos[2]), color="red") +
      geom_segment(aes(x = 2.5, y = actual_results_pos[3], xend = 3.5, yend = actual_results_pos[3]), color="red") +
      geom_segment(aes(x = 3.5, y = actual_results_pos[4], xend = 4.5, yend = actual_results_pos[4]), color="red") +
      geom_segment(aes(x = 4.5, y = actual_results_pos[5], xend = 5.5, yend = actual_results_pos[5]), color="red") +
      geom_segment(aes(x = 5.5, y = actual_results_pos[6], xend = 6.5, yend = actual_results_pos[6]), color="red") +
      geom_segment(aes(x = 6.5, y = actual_results_pos[7], xend = 7.5, yend = actual_results_pos[7]), color="red") +
      geom_segment(aes(x = 0.5, y = actual_results_neg[1], xend = 1.5, yend = actual_results_neg[1]), color="blue") + 
      geom_segment(aes(x = 1.5, y = actual_results_neg[2], xend = 2.5, yend = actual_results_neg[2]), color="blue") +
      geom_segment(aes(x = 2.5, y = actual_results_neg[3], xend = 3.5, yend = actual_results_neg[3]), color="blue") +
      geom_segment(aes(x = 3.5, y = actual_results_neg[4], xend = 4.5, yend = actual_results_neg[4]), color="blue") +
      geom_segment(aes(x = 4.5, y = actual_results_neg[5], xend = 5.5, yend = actual_results_neg[5]), color="blue") +
      geom_segment(aes(x = 5.5, y = actual_results_neg[6], xend = 6.5, yend = actual_results_neg[6]), color="blue") +
      geom_segment(aes(x = 6.5, y = actual_results_neg[7], xend = 7.5, yend = actual_results_neg[7]), color="blue")
    #  ggtitle(paste0("Spin Test Perm: ", model))
    ggsave(plot=plot_violin, filename = paste0(homedir, "/baller/results/images/spin_", model, "_pos_and_neg_lines_t_fdr05.png"), width = 4.81, height = 4.81)
    ggsave(plot=plot_violin, filename = paste0(homedir, "/baller/results/images/spin_", model, "_pos_and_neg_lines_t_fdr05.pdf"), width = 4.81, height = 4.81)
   }
}

violin_plot_pos_and_neg_lines_with_color_gradation <- function (homedir, models, network_names, num_spins){
  for (model in models) {
    
    print(model)
    #for storing statistics at the end
    lh_spin_df <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/lh_spin_test_", model, "_proportions.csv"), sep = ",")
    )
    lh_spin_pos <-data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/lh_spin_test_pos_", model, "_proportions.csv"), sep = ",")
    )
    lh_spin_neg <-data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/lh_spin_test_neg_", model, "_proportions.csv"), sep = ",")
    )
    rh_spin_df <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/rh_spin_test_", model, "_proportions.csv"), sep = ",")
    )
    rh_spin_pos <-data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/rh_spin_test_pos_", model, "_proportions.csv"), sep = ",")
    )
    rh_spin_neg <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/rh_spin_test_neg_", model, "_proportions.csv"), sep = ",")
    )
    
    #take means of left and right
    actual_results <- (lh_spin_df[,1] + rh_spin_df[,1])/2
    actual_results_pos <- ((lh_spin_pos[,1] + rh_spin_pos[,1])/2)
    actual_results_neg <- ((lh_spin_neg[,1] + rh_spin_neg[,1])/2)
    
    # print(actual_results)
    #  print(actual_results_pos)
    #  print(actual_results_neg)
    #dataframes for all spins, as well as positive and negative
    spin_without_target_col <- cbind(lh_spin_df[,2:1001],rh_spin_df[,2:1001])
    # print(min(spin_without_target_col))
    #  print(spin_without_target_col[4,5:10])
    
    melted_df <- melt_df_for_violin_plot_yeo7(spin_without_target_col, network_names, num_spins)
    
    #  print(melted_df[1:4,2:3])
    #with mean lines, different fonts
    #save images
  plot_violin <- ggplot(melted_df, aes(x = factor(network_name, level = network_names), y = spin, fill = network_name)) +  
      scale_fill_manual(values=yeo_colors) + 
      geom_violin(trim = TRUE) + 
      xlab("Yeo 7 Network") + ylab(paste0("Proportion")) +
      ylim(0,NA) + 
      theme_classic() + 
      theme(legend.position = "none",
            legend.title = element_blank(),
            axis.text.x = element_text(size = 10, colour = "black"),
            axis.text.y = element_text(size = 10, colour = "black"),
            axis.title.y = element_text(size = 10),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 10)) +
      stat_summary(fun.y = mean, geom = "errorbar", 
                   aes(ymax = ..y.., ymin = ..y.., group = factor(network_name)),
                   width = 0.5, linetype = "dashed", position = position_dodge(0.9)) + 
      geom_segment(aes(x = 0.6, y = actual_results_pos[1], xend = 1.4, yend = actual_results_pos[1]), color="red") + 
      geom_segment(aes(x = 1.6, y = actual_results_pos[2], xend = 2.4, yend = actual_results_pos[2]), color="red") +
      geom_segment(aes(x = 2.6, y = actual_results_pos[3], xend = 3.4, yend = actual_results_pos[3]), color="red") +
      geom_segment(aes(x = 3.6, y = actual_results_pos[4], xend = 4.4, yend = actual_results_pos[4]), color="red") +
      geom_segment(aes(x = 4.6, y = actual_results_pos[5], xend = 5.4, yend = actual_results_pos[5]), color="red") +
      geom_segment(aes(x = 5.6, y = actual_results_pos[6], xend = 6.4, yend = actual_results_pos[6]), color="red") +
      geom_segment(aes(x = 6.6, y = actual_results_pos[7], xend = 7.4, yend = actual_results_pos[7]), color="red") +
      geom_segment(aes(x = 0.6, y = actual_results_neg[1], xend = 1.4, yend = actual_results_neg[1]), color="blue") + 
      geom_segment(aes(x = 1.6, y = actual_results_neg[2], xend = 2.4, yend = actual_results_neg[2]), color="blue") +
      geom_segment(aes(x = 2.6, y = actual_results_neg[3], xend = 3.4, yend = actual_results_neg[3]), color="blue") +
      geom_segment(aes(x = 3.6, y = actual_results_neg[4], xend = 4.4, yend = actual_results_neg[4]), color="blue") +
      geom_segment(aes(x = 4.6, y = actual_results_neg[5], xend = 5.4, yend = actual_results_neg[5]), color="blue") +
      geom_segment(aes(x = 5.6, y = actual_results_neg[6], xend = 6.4, yend = actual_results_neg[6]), color="blue") +
      geom_segment(aes(x = 6.6, y = actual_results_neg[7], xend = 7.4, yend = actual_results_neg[7]), color="blue")
    #  ggtitle(paste0("Spin Test Perm: ", model))
    ggsave(plot=plot_violin, filename = paste0(homedir, "/baller/results/images/spin_", model, "_pos_and_neg_lines_t_fdr05_test_gradation.png"), width = 4.81, height = 4.81)
    ggsave(plot=plot_violin, filename = paste0(homedir, "/baller/results/images/spin_", model, "_pos_and_neg_lines_t_fdr05_test_gradation.pdf"), width = 4.81, height = 4.81)
    
  }
}

violin_plot_means <- function (homedir, models, network_names, num_spins){
  for (model in models) {
    #for storing statistics at the end
    lh_spin_df <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/lh_spin_test_", model, "_proportions.csv"), sep = ",")
    )
    rh_spin_df <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/rh_spin_test_", model, "_proportions.csv"), sep = ",")
    )
  
    
    #take means of left and right
    actual_results <- (lh_spin_df[,1] + rh_spin_df[,1])/2
    
    #dataframes for all spins
    spin_without_target_col <- cbind(lh_spin_df[,2:1001],rh_spin_df[,2:1001])
    melted_df <- melt_df_for_violin_plot_yeo7(spin_without_target_col, network_names, num_spins)
  
    #with mean lines, different fonts
    #save images
    plot_violin <- ggplot(melted_df, aes(x = factor(network_name, level = network_names), y = spin, fill = network_name)) +  
      scale_fill_manual(values=yeo_colors) + 
      geom_violin(trim = TRUE) + 
      xlab("Yeo 7 Network") + ylab(paste0("Proportion")) +
      ylim(0,NA) + 
      theme_classic() + 
      theme(legend.position = "none",
            legend.title = element_blank(),
            axis.text.x = element_text(size = 10, colour = "black"),
            axis.text.y = element_text(size = 10, colour = "black"),
            axis.title.y = element_text(size = 10),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 10)) +
     # stat_summary(fun.y = mean, geom = "errorbar", 
      #             aes(ymax = ..y.., ymin = ..y.., group = factor(network_name)),
       #            width = 0.5, linetype = "dashed", position = position_dodge(0.9)) + 
      #geom_boxplot(width = 0.15, position = position_dodge(0.9)) + 
        geom_segment(aes(x = 0.6, y = actual_results[1], xend = 1.4, yend = actual_results[1])) + 
        geom_segment(aes(x = 1.6, y = actual_results[2], xend = 2.4, yend = actual_results[2])) +
        geom_segment(aes(x = 2.6, y = actual_results[3], xend = 3.4, yend = actual_results[3])) +
        geom_segment(aes(x = 3.6, y = actual_results[4], xend = 4.4, yend = actual_results[4])) +
        geom_segment(aes(x = 4.6, y = actual_results[5], xend = 5.4, yend = actual_results[5])) +
        geom_segment(aes(x = 5.6, y = actual_results[6], xend = 6.4, yend = actual_results[6])) +
        geom_segment(aes(x = 6.6, y = actual_results[7], xend = 7.4, yend = actual_results[7]))
  #    ggtitle(paste0("Spin Test Perm: ", model))
    ggsave(plot=plot_violin, filename = paste0(homedir, "/baller/results/images/spin_", model, "_t_fdr05.png"), width = 4.81, height = 4.81)
    ggsave(plot=plot_violin, filename = paste0(homedir, "/baller/results/images/spin_", model, "_t_fdr05.pdf"), width = 4.81, height = 4.81)
    
  }
}

violin_plot_means_mean_coupling <- function (homedir, network_names, num_spins){
 
    #for storing statistics at the end
    lh_spin_df <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/lh_spin_test_mean_coupling.csv"), sep = ","))
    rh_spin_df <- data.frame(read.table(paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/rh_spin_test_mean_coupling.csv"), sep = ","))
    
    
    #take means of left and right
    actual_results <- (lh_spin_df[,1] + rh_spin_df[,1])/2
    
    #dataframes for all spins
    spin_without_target_col <- cbind(lh_spin_df[,2:1001],rh_spin_df[,2:1001])
    melted_df <- melt_df_for_violin_plot_yeo7(spin_without_target_col, network_names, num_spins)
    
    #with mean lines, different fonts
    #save images
    plot_violin <- ggplot(melted_df, aes(x = factor(network_name, level = network_names), y = spin, fill = network_name)) +  
      scale_fill_manual(values=yeo_colors) + 
      geom_violin(trim = TRUE) + 
      xlab("Yeo 7 Network") + ylab(paste0("Mean")) +
      ylim(0,NA) + 
      theme_classic() + 
      theme(legend.position = "none",
            legend.title = element_blank(),
            axis.text.x = element_text(size = 10, colour = "black"),
            axis.text.y = element_text(size = 10, colour = "black"),
            axis.title.y = element_text(size = 10),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 10)) +
      stat_summary(fun.y = mean, geom = "errorbar", 
                   aes(ymax = ..y.., ymin = ..y.., group = factor(network_name)),
                   width = 0.5, linetype = "dashed", position = position_dodge(0.9)) + 
      #geom_boxplot(width = 0.15, position = position_dodge(0.9)) + 
      geom_segment(aes(x = 0.6, y = actual_results[1], xend = 1.4, yend = actual_results[1])) + 
      geom_segment(aes(x = 1.6, y = actual_results[2], xend = 2.4, yend = actual_results[2])) +
      geom_segment(aes(x = 2.6, y = actual_results[3], xend = 3.4, yend = actual_results[3])) +
      geom_segment(aes(x = 3.6, y = actual_results[4], xend = 4.4, yend = actual_results[4])) +
      geom_segment(aes(x = 4.6, y = actual_results[5], xend = 5.4, yend = actual_results[5])) +
      geom_segment(aes(x = 5.6, y = actual_results[6], xend = 6.4, yend = actual_results[6])) +
      geom_segment(aes(x = 6.6, y = actual_results[7], xend = 7.4, yend = actual_results[7]))
      #    ggtitle(paste0("Spin Test Perm: ", model))
      ggsave(plot=plot_violin, filename = paste0(homedir, "/baller/results/images/spin_mean_coupling_unc.png"), width = 4.81, height = 4.81)
      ggsave(plot=plot_violin, filename = paste0(homedir, "/baller/results/images/spin_mean_coupling_unc.pdf"), width = 4.81, height = 4.81)
      
}


get_yeo7_colors <- function() {
  yeo_colors <- c(
    `VIS` = "#781286",
    `MOT` = "#4682b4",
    `DA` = "#00760e",
    `VA` = "#c43afa",
    `LIM` = "#dcf8a4",
    `FP` = "#e69422",
    `DM` = "#cd3e56")
  return(yeo_colors)
}

get_yeo7_colors_faded <- function() {
  yeo_colors <- c(
    `VIS` = "#b870c2",
    `MOT` = "#709dc2",
    `DA` = "#70c27a",
    `VA` = "#ab70c2",
    `LIM` = "#eff8de",
    `FP` = "#e6be85",
    `DM` = "#c2707e")
  return(yeo_colors)
}

get_yeo_illustrator_colors <- function(){
  yeo_colors <- c(
    `VIS` = "#781286",
    `MOT` = "#709dc2",
    `DA` = "#70c27a",
    `VA` = "#ab70c2",
    `LIM` = "#eff8de",
    `FP` = "#e6be85",
    `DM` = "#c2707e")
  return(yeo_colors)
   781286
}
get_network_colors_with_stats_yeo7 <- function(stats){
  #takes the statistics for a yeo7 matrix and returns a color palatte for plotting
  #white for non-significant
  #dark color for significant
  #faded color for trend
  
}