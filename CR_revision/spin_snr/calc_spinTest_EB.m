% addpaths
addpath(genpath('/project/imco/baller/scripts/spin_test/'));
ProjectFolder = '/project/imco/baller/results/coupling_accuracy/';
%model = 'gam_sex'
model = 'gam_exec_accuracy'

% set outdir
outdir= '/project/imco/baller/results/coupling_accuracy/spin_test_results/';
outFn=strcat([outdir, '/', model, '_spin_results_yeo_1_0_-1_output']);

SpinPermuFS([ProjectFolder, '/lh_', model, '_t_fdr05_Yeo7_1_0_-1.csv'],[ProjectFolder, '/rh_', model, '_t_fdr05_Yeo7_1_0_-1.csv'], 1000, outFn, model);

