#!/bin/bash
# set freesurfer version
FREESURFER_HOME=/share/apps/freesurfer/6.0.0

# useage
function Usage {
    cat <<USAGE

Usage:
  -i:  input list of bblid, scanid, image to project, and path to coreg file. required.
  -s:  freesurfer processed directory. required. 
  -o:  output directory. required.

Example: ./vol2surf_wrapper.sh -i /data/joy/BBL/tutorials/exampleData/vol2surf/subjList.csv -s /data/joy/BBL/studies/pnc/processedData/structural/freesurfer53/ -o /home/lbeard/vol2surf

USAGE
    exit 1
}

# check for arguments
[ $# -lt 1 ] && Usage

# read command line arguments
while getopts "i:s:o:" OPT
	do
	case $OPT in
		h) #help
		Usage >&2
		exit 0
		;;
		i) #input file
		subj_list=${OPTARG}
		;;
		s) #freesurfer directory
		SUBJECTS_DIR=${OPTARG}
		;;	
		o) #output directory
		targDir=${OPTARG}
		;;
		\?) # getopts issues an error message
		Usage >&2
		exit 1
		;;
	esac
done

# check for required arguments
[ -z $subj_list ] || [ -z $SUBJECTS_DIR ] || [ -z $targDir ] && echo "subject list, freesurfer directory and output directory are required." && Usage

#submit projection script
qsub -q himem.q -V -N vol2surf -v subj_list=${subj_list} -v SUBJECTS_DIR=${SUBJECTS_DIR} -v targDir=${targDir} ./vol2surf.sh
