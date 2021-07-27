function PBP_vertWiseEffect_Erica(LH,RH,name)% pretty picture code, AAB 4/2018 - AP 5/1/20 - Updated to threshold according to input - 1/12/21
% data should be vectors, 10242 in length if fsaverage5 is used
% if using higher resolution, then change accordingly
% depencies include: matlab freesrufer functions, subaxis.m (matlab central), inferno color scale (matlab central - for Sam ;)

%%% SET THRESHOLDS AS DESIRED HERE: only fill in each threshold as needed (no need to set both if you only want to threshold one end)
% Values at or above this set to gray
 %Uthresh=-2;
% Values at or below this set to gray
 LThresh=2;
%%%


addpath(genpath('/appl/freesurfer-6.0.0/matlab/'));
addpath(genpath('/project/imco/baller/scripts/subaxis/'));
addpath(genpath('/project/imco/baller/scripts/Colormaps/Colormaps (5)/Colormaps/'));
%{
addpath(genpath('/cbica/projects/alpraz_EI/scripts/tools/'));
ProjectFolder = '/cbica/projects/pinesParcels/data/SingleParcellation';
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage5';
%}

plot_text='';
[vertices, faces] = freesurfer_read_surf('/project/imco/surfaces/fsaverage5/surf/lh.inflated');
%using lh.gray will make more anatomical looking plot but harder to see into sulci
right = readtable(RH,'TreatAsEmpty','NA','ReadVariableNames',false);
datar = table2array(right);
left = readtable(LH,'TreatAsEmpty','NA','ReadVariableNames',false);
datal = table2array(left);

%left=load(LHvec);
%right=load(RHvec);
%datal=left;
%datar=right;
minval = min(min(datal),min(datar)) %useful for colorbar later	

%set NaN to 0
%I generally have the midcut region set to NaN
%in the csv files that I read in
indexNaNrh = find(isnan(datar));
indexNaNlh = find(isnan(datal));
datar(indexNaNrh)=0;
datal(indexNaNlh)=0;
datalr=[datal; datar];
%invoke thresholding 1/12/21
if exist('Uthresh','Var') == 1;
	AboveThresh= datalr > Uthresh;
	datalr(AboveThresh)=0;
end
if exist('LThresh','Var') ==1;
	BelowThresh= datalr < LThresh;
	datalr(BelowThresh)=0;
end
%%% set color scale
% 1/12/21 - for p values, visualizing 1/p might be more effective. comment out line below and  uncomment subsequent line to nix this approach.
%datalr=1./datalr;


% 12/1/21 tiny bit of code to deal with 1/0 in matlab
InfIndex=find(datalr==Inf);
% Infinity values to 0
datalr(InfIndex)=0;

%AP% set to make white zero on all maps
maxabs=prctile(abs(datalr),88);
%mincol= minval-.00001 
%maxcol=maxabs
%mincol=-maxabs
maxcol=max(datalr)
mincol=min(datalr)
%change above to set max/min manually or by other means
%custommap=colormap('plasma'); %or whatever
% for white at 0
%custommap=colormap(b2r(-1,1));
%custommap=colormap('jet');
custommap=colormap('plasma')
custommap(1,:)=[0.75 0.75 0.75];


data=datalr(1:10242);
asub = subaxis(4,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);
%asub = subplot(4,2,1)
% note use of subaxis is to ged rid of white space around brains 
% if you don't care about that, it's faster and less likely to cause
% issues if you use subplot instead
% if so, bet rid of all of the posnew stuff below

aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting gouraud; %phong; 
material metal %shiny %metal; 
shading flat;
camlight;
alpha(1)

asub = subaxis(4,2,3, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
rotate(aplot, [0 0 1], 180)
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting gouraud; %phong; 
material metal %shiny %metal; 
shading flat;
camlight;
alpha(1)
set(gcf,'Color','w')

asub = subaxis(4,2,5, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
rotate(aplot, [0 0 1], 225)
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting gouraud; %phong; 
material metal %shiny %metal; 
shading flat;
camlight;
alpha(1)
set(gcf,'Color','w')

asub = subaxis(4,2,7, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0, 'MT', 0.0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
axis vis3d off;
rotate(aplot, [0 1 0], 270)
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
lighting gouraud; %phong; 
material metal %shiny %metal; 
shading flat;
camlight;
alpha(1)
set(gcf,'Color','w')

 pos = get(asub, 'Position');
 posnew = pos; posnew(2) = posnew(2) + 0.04; set(asub, 'Position', posnew);
 %white space again

%plot title 
title(plot_text)
set(get(gca,'title'),'Position',[332 119 3])

%%% right hemisphere
data=datalr(10243:20484);

[vertices, faces] = freesurfer_read_surf('/project/imco/surfaces/fsaverage5/surf/rh.inflated');

asub = subaxis(4,2,2, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
rotate(aplot, [0 0 1], 180)
colormap(custommap)
caxis([mincol; maxcol]);
%caxis([NAval; max_data])
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting phong; %gouraud
material metal %shiny %metal; 
shading flat;
camlight;
alpha(1)
%colormap(mycol)


 pos = get(asub, 'Position');
 posnew = pos; posnew(1) = posnew(1) - 0.22; set(asub, 'Position', posnew);

asub = subaxis(4,2,4, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting gouraud; %phong; 
material metal %shiny %metal; 
shading flat;
camlight;
alpha(1)
set(gcf,'Color','w')
 pos = get(asub, 'Position');
 posnew = pos; posnew(1) = posnew(1) - 0.22; set(asub, 'Position', posnew);

asub = subaxis(4,2,6, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);

aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
rotate(aplot, [0 0 1], -45)
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting gouraud; %phong; 
material metal %shiny %metal; 
shading flat;
camlight;
alpha(1)
set(gcf,'Color','w')
 pos = get(asub, 'Position');
 posnew = pos; posnew(1) = posnew(1) - 0.22; set(asub, 'Position', posnew);


%%%
asub = subaxis(4,2,8, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);

aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
axis vis3d off;
rotate(aplot, [0 1 0], 270)
rotate(aplot, [1 0 0], 180)
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
lighting gouraud; %phong; 
material metal %shiny %metal; 
shading flat;
camlight;
alpha(1)
%set(gcf,'Color',[.2 .2 .2])
set(gcf,'Color',[1,1,1])
 pos = get(asub, 'Position');
 posnew = pos; posnew(2) = posnew(2) + 0.04; set(asub, 'Position', posnew);
 pos = get(asub, 'Position');
 posnew = pos; posnew(1) = posnew(1) - 0.22; set(asub, 'Position', posnew);
%%%


acbar = colorbar('EastOutside')
set(acbar, 'position', [0.40 0.270 0.02 0.20])


% going lower rez for now, but giant vector rendering was beaut
print('-dpng','-r600',['~/' char(name)])
