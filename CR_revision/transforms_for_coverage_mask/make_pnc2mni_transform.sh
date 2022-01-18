#!/bin/bash

#moving cbfxalff into mni space, then projecting to fsaverage5 space

#Parcellation and Template Paths
transforms_path='/project/imco/baller/transforms_for_surface_projection/'
pnc_path=${transforms_path}/space/PNC/PNC_transforms
cbfxalff_mask=$transforms_path/cbfxalff_mask.nii
cbfxalff_mask_mni=$transforms_path/cbfxalff_mask_mni.nii
#mni_2mm_path=$transforms_path/space/MNI/MNI-2x2x2.nii.gz
mni_2mm_path=${FSLDIR}//data/standard/MNI152_T1_2mm_brain.nii.gz
affine_mat=$pnc_path/PNC-MNI_1Affine.mat
pnc2mni_warp=$pnc_path/PNC-MNI_0Warp.nii.gz


# set parameters
hemis=(lh rh)

#################################
## Create 2mm MNI_path ##
#################################
echo "applying antsApplyTransforms"
#antsApplyTransforms -e 3 -d 3 -i ${cbfxalff_mask} -o ${cbfxalff_mask_mni} -r ${mni_2mm_path} -t ${pnc2mni_warp} -t ${affine_mat} -n MultiLabel
################################
## Project vol to surface ##
################################

for hemi in ${hemis[@]}; do		
#	mri_vol2surf --mov  $cbfxalff_mask_mni  --regheader fsaverage5 --hemi ${hemi} --o ${hemi}.surf.mgh  --projfrac-avg 0 1 0.1 --surf white
	mris_convert -c ${transforms_path}/${hemi}.surf.mgh ${SUBJECTS_DIR}/fsaverage5/surf/${hemi}.white ${transforms_path}/${hemi}.surf.asc	
#	more ${transforms_path}/${hemi}.surf.asc | cut -d ' ' -f 5 > ${transforms_path}/${hemi}_cbfxalff_mask_binary
#mri_surf2surf --srcsubject ${transforms_path} --trgsubject ico --trgicoorder 5 ---srcsurfval ${transforms_path}/cbfxalff_mask_mni_${hemi}_surf.mgh -trgsurfval ${transforms_path}/cbfxalff_mask_mni_${hemi}_fs5_surf.mgh --hemi ${hemi} 
done
################################
## Surf to fsaverage5 ##
################################
#mri_surf2surf 
