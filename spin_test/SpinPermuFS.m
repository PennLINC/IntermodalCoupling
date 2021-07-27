function SpinPermuFS(readleft,readright,permno,wsname,model)
% Compute designated # of permutations/spins of the input surface data
% in FreeSurfer fsaverage5.
% FORMAT SpinPermuFS(readleft,readright,permno)
% readleft     - the filename of left surface data to spin 
% readright    - the filename of right surface data to spin 
% permno       - the number of permutations
% wsname       - the name of a workspace file including all spun data to be saved
% model        - gam_accuracy, lm_accuracy, gam_sex, etc
% Example   SpinPermuFS('../data/depressionFSdataL.csv','../data/depressionFSdataR.csv',100,'../data/rotationFS.mat')
% will spin prebuilt data, neurosynth map associated with 'depression', 100
% times, and save the workspace file of all spun data in ../data/rotationFS.mat
% Aaron Alexander-Bloch & Siyuan Liu 
% SpinPermuFS.m, 2018-04-22
% The implementation of generating random rotations originally described in our paper — 
% rotating the coordinates of vertices at angles uniformly chosen between zero and 360 degrees
% about each of the x (left-right), y (anterior-posterior) and z (superior-inferior) axes —
% introduces a preference towards oversampling certain rotations. 
% Thus, we modified the code to incorporate an approach, Lefèvre et al. (2018), 
% that samples uniformly from the space of possible rotations. The updated
% uniform sampling prodcedure does not require AxelRot.m anymore.
% Updated on 2018-07-18


outdir= '/project/imco/baller/results/coupling_accuracy/spin_test_results/';

%Set up paths
fshome = getenv('FREESURFER_HOME');
fsmatlab = sprintf('%s/matlab',fshome);
path(path,fsmatlab);
%read the data saved in csv
datal=importdata(readleft);
datar=importdata(readright);


%%extract the correspoding sphere surface coordinates for rotation
% AP changed the filepath 10/26/20, i don't know how cubic managed to set up a fshome path not amenable to this
%[verticesl, ~] = freesurfer_read_surf(fullfile('/project/imco/surfaces/fsaverage5/surf/lh.sphere'));
%[verticesr, ~] = freesurfer_read_surf(fullfile('/project/imco/surfaces/fsaverage5/surf/rh.sphere'));
[verticesl, ~] = freesurfer_read_surf(fullfile('/appl/freesurfer-7.1.1/subjects/fsaverage5/surf/lh.sphere'));
[verticesr, ~] = freesurfer_read_surf(fullfile('/appl/freesurfer-7.1.1/subjects/fsaverage5/surf/rh.sphere'));


rng(0);
%Use rng to initialize the random generator for reproducible results.
%initialize variables to save rotation
bigrotl=[];
bigrotr=[];
distfun = @(a,b) sqrt(bsxfun(@minus,bsxfun(@plus,sum(a.^2,2),sum(b.^2,1)),2*(a*b)));
%function to calculate Euclidian distance
I1 = eye(3,3);
I1(1,1)=-1;
bl=verticesl;
br=verticesr;
%permutation starts
for j=1:permno
    j
    %the updated uniform sampling procedure
    A = normrnd(0,1,3,3);
    [TL, temp] = qr(A);
    TL = TL * diag(sign(diag(temp)));
    if(det(TL)<0)
        TL(:,1) = -TL(:,1);
    end
    %reflect across the Y-Z plane for right hemisphere
    TR = I1 * TL * I1;
    bl =bl*TL;
    br = br*TR;    
    
    %Find the pair of matched vertices with the min distance and reassign
    %values to the rotated surface.
    distl=distfun(verticesl,bl');
    distr=distfun(verticesr,br');
    [~, Il]=min(distl,[],2);
    [~, Ir]=min(distr,[],2);
    %save rotated data
    bigrotl=[bigrotl; datal(Il)'];
    bigrotr=[bigrotr; datar(Ir)'];
    % it is also feasible to save Il Ir and apply them to different datasets
    % for repeated use
    %If annotation file is used, annotation file for each rotation could be
    %saved by write_annotation.m of FreeSurfer
end
save(wsname,'bigrotl','bigrotr')
writematrix(bigrotl, [outdir, '/lh_spin_test_', model, '_output.csv']);
writematrix(bigrotr, [outdir, '/rh_spin_test_', model, '_output.csv']);
%save bigrotl and bigrotr in a workspace file for the null distribution
%use it in pvalvsNull.m to caclulate pvalue
