#!/bin/sh

##################################
### PBP visualizations Wrapper ###
##################################

###### Author: Erica Baller ######
######   Date: 2/26/2021    ######

## pre: commands_for_matlab (in same directory as PBP scripts
   #- contains all the commands we want to pass to matlab to run visualizations automaticall
## post: images in /project/imco/baller/results/images/pbp
## uses: I was getting tired of having to open matlab and run each of these commands manually. This script takes a bunch of commands and just feeds them to matlab, no muss, no fuss
## dependencies: Matlab 2020b, please run "source /project/imco/baller/scripts/load_matlab before starting this


#set directories
homedir='/project/imco'
#set homedir = '/Users/eballer/BBL/imco/pmacs/PMACS_remote/'
wking_dir='/baller/scripts/CR_revision/PBP_graphics'
command_file="$homedir/$wking_dir/commands_for_matlab"
echo $command_file

#initialize matlab


#######
### loop ###
while IFS= read -r line; do
	echo "$line"
	matlab -nosplash -nodesktop -nodisplay -r "$line; exit"

done < $command_file


