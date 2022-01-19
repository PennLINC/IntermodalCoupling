#!/bin/bash

# Things to change based on project: SUBJECTS_DIR, fwhm, scripts, sublist, nameOfWrapper (surfCoup_wrapper.sh)
SUBJECTS_DIR=/data/joy/BBL/studies/pnc/processedData/structural/freesurfer53/
fwhm=15
scripts=/data/jux/BBL/projects/coupling/coupling_test_eb_20200918/code
sublist=/data/jux/BBL/projects/coupling/coupling_test_eb_20200918/subjectLists/n5_surfCoupling_list
logdir=$scripts/logdir
outdir=/data/jux/BBL/projects/coupling/coupling_test_eb_20200918/surfCoupling

mkdir $logdir 2>/dev/null
rm -f $logdir/*
queue="himem.q"
time=$(date | cut -d " " -f4 | sed "s+:++g")


nsub=$(cat $sublist | wc -l )
qsub -V -N coupling_$time -t 1 -cwd -o $logdir -e $logdir -q $queue $scripts/surfCoup_wrapper.sh -s $sublist -f $fwhm -o $outdir
qsub -V -t 2-$nsub -cwd -hold_jid coupling_$time -o $logdir -e $logdir -q $queue $scripts/surfCoup_wrapper.sh -s $sublist -f $fwhm -o $outdir
