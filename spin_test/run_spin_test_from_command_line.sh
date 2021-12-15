#!/bin/sh

##################################
###    Spin Test    Wrapper    ###
##################################

###### Author: Erica Baller ######
######   Date: 2/26/2021    ######

## pre: commands_for_matlab 
   #- trinarized maps
## post: spin results in /project/imco/baller/results/coupling_accuracy/spin_test_results
## uses: I was getting tired of having to open matlab and run each of these commands manually. This script takes a bunch of commands and just feeds them to matlab, no muss, no fuss
## dependencies: Matlab 2020b, please run "source /project/imco/baller/scripts/load_matlab before starting this


#set directories
homedir='/project/imco'
#set homedir = '/Users/eballer/BBL/imco/pmacs/PMACS_remote/'
wking_dir='/baller/scripts/spin_test'
command_file="$homedir/$wking_dir/commands_for_matlab"
echo $command_file

#initialize matlab


#######
### loop ###
while IFS= read -r line; do
	echo "$line"
	echo "here"
	matlab -nosplash -nodesktop -nodisplay -r "$line; exit"

done < $command_file


