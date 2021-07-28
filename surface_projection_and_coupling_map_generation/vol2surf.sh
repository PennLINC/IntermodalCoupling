#!/bin/bash

# set parameters
hemi=(lh rh)
bold=$(tput bold)
normal=$(tput sgr0)

# conversion, projection, resampling
for i in `cat $subj_list`; do
	bblid=`echo $i | cut -d "," -f 1`
	datexscanid=`echo $i | cut -d "," -f 2`
	funcImage=`echo $i | cut -d "," -f 3`
	boldDir=`echo $i | cut -d "," -f 4`

	echo "${bold}start $bblid/$datexscanid"
	
	cd $targDir
	
	if [ ! -e "${bblid}/${datexscanid}" ]; then
		mkdir -p ${bblid}/${datexscanid} && cd ${bblid}/${datexscanid}
	else
		cd ${bblid}/${datexscanid}
	fi

	# convert XCP coreg file to freesurfer format	
	lta_convert --infsl ${boldDir}/${bblid}_${datexscanid}_seq2struct.mat \
	--outreg ${targDir}/${bblid}/${datexscanid}/${bblid}_${datexscanid}_seq2struct.dat \
	--subject ${bblid}/${datexscanid} \
	--src ${boldDir}/${bblid}_${datexscanid}_referenceVolume.nii.gz \
	--trg ${boldDir}/${bblid}_${datexscanid}_targetVolume.nii.gz \
	>> ${targDir}/${bblid}/${datexscanid}/${bblid}_${datexscanid}_projection.log 2>&1

	# write script to transform the matrix: (switch 2nd/3rd column; multiply column 2 by -1)
	Rscript /data/joy/BBL/tutorials/code/vol2surf/transformMatrix.R
	mv  ${targDir}/${bblid}/${datexscanid}/matrixFinalTransform.dat ${targDir}/${bblid}/${datexscanid}/${bblid}_${datexscanid}_matrixFinalTransform.dat
	rm ${targDir}/${bblid}/${datexscanid}/${bblid}_${datexscanid}_seq2struct.dat
	
	# surface projection #	
	for j in ${hemi[@]}; do	

		# projection from volume to freesurfer native surfaces				
		echo "   ${normal}start vol2surf ${j}"	
		mri_vol2surf --mov $funcImage \
		--hemi ${j} --projfrac 0.5 --interp trilinear --noreshape \
		--reg ${targDir}/${bblid}/${datexscanid}/${bblid}_${datexscanid}_matrixFinalTransform.dat \
		--o ${targDir}/${bblid}/${datexscanid}/${bblid}_${datexscanid}_${j}_surf.mgh \
		>> ${targDir}/${bblid}/${datexscanid}/${bblid}_${datexscanid}_projection.log 2>&1		
		echo "   done vol2surf ${j}"

		# resampling from freesurfer native surface to freesurfer fsaverage5 atlas surface 	
		echo "   start surf2surf  ${j}"
		mri_surf2surf --srcsubject ${bblid}/${datexscanid} \
		--trgsubject ico --trgicoorder 5 --hemi ${j} \
		--srcsurfval ${targDir}/${bblid}/${datexscanid}/${bblid}_${datexscanid}_${j}_surf.mgh \
		--trgsurfval ${targDir}/${bblid}/${datexscanid}/${bblid}_${datexscanid}_${j}_fs5_surf.mgh \
		>> ${targDir}/${bblid}/${datexscanid}/${bblid}_${datexscanid}_projection.log 2>&1
		echo "   done surf2surf ${j}"
	
		done

	echo "${bold}done ${bblid}/${datexscanid}"
	echo ""

done
