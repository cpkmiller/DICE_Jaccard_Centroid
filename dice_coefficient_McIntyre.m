clc
clear all
%% dice_coefficient_McIntyre.m
% Jennifer Muller, written Jul 20 2021

% ROI
% (STN_L.nii, GPi_L.nii, redNucleus_L.nii, etc.)
mROI = {'GPi_L.nii' 'GPi_R.nii' 'STN_L.nii' 'STN_R.nii' 'RN_L.nii' 'RN_R.nii'};
cmROI = {'CM_GPI L.nii' 'CM_GPI R.nii' 'CM_STN L.nii' 'CM_STN R.nii' 'CM_Red Nucleus L.nii' 'CM_Red Nucleus R.nii'};
wROI = {'Wu_L_GPi.nii' 'Wu_R_GPi.nii' 'Wu_L_STN.nii' 'Wu_R_STN.nii' 'Wu_L_RN.nii' 'Wu_R_RN.nii'};
%wROI =  {'GPi_L_SynDeformed.nii' 'GPi_R_SynDeformed.nii' 'STN_L_SynDeformed.nii' 'STN_R_SynDeformed.nii' 'RN_L_SynDeformed.nii' 'RN_R_SynDeformed.nii'};
csROI = {'CS_GPI L.nii' 'CS_GPI R.nii' 'CS_STN L.nii' 'CS_STN R.nii' 'CS_Red Nucleus L.nii' 'CS_Red Nucleus R.nii'}; 
atlasROI = {'GPi_L_SynDeformed.nii' 'GPi_R_SynDeformed.nii' 'STN_L_SynDeformed.nii' 'STN_R_SynDeformed.nii' 'RN_L_SynDeformed.nii' 'RN_R_SynDeformed.nii'}; 

% subject directory
dir = '/Users/jxm491/Dropbox/Mac (3)/Desktop/New_PD04';

cd(dir)

% McIntyre directory
McIntyre = 'McIntyre/';

% Wu directory
Wu = 'Wu/';

% Miller directory
Miller = 'CM Segmentation/';

% Cranial Suite
CS = 'Cranial Suite Segmentation/';

% Atlas directory
Atlas_dir = 'Atlas/';

%% read and store files


dice_coef = zeros(3,numel(mROI));

for count=1:numel(mROI);
    
    %McIntyre=McInt
    
m = MRIread(strcat(McIntyre,mROI{count}));
m1 = m.vol;
[r c z] = size(m1);
for i = 1:z
    out1(:,:,i) = im2bw(m1(:,:,i));
end
m.vol = out1;

w = MRIread(strcat(Wu,wROI{count}));
w1 = w.vol;
[r c z] = size(w1);
for i = 1:z
    out2(:,:,i) = im2bw(w1(:,:,i));
end
w.vol = out2;

mil = MRIread(strcat(Miller,cmROI{count}));
mil1 = mil.vol;
[r c z] = size(mil1);
for i = 1:z
    out3(:,:,i) = im2bw(mil1(:,:,i));
end
mil.vol = out3;

cran = MRIread(strcat(CS,csROI{count}));
cran1 = cran.vol;
[r c z] = size(cran1);
for i = 1:z
    out4(:,:,i) = im2bw(cran1(:,:,i));
end
cran.vol = out4;

Atlas = MRIread(strcat(Atlas_dir,atlasROI{count}));
Atlas1 = Atlas.vol;
[r c z] = size(Atlas1);
for i = 1:z
    out5(:,:,i) = im2bw(Atlas1(:,:,i));
end
Atlas.vol = out5;

%% dice coefficients

dc_union = []

seg1 = m.vol(:);
seg2 = w.vol(:);
seg3 = mil.vol(:);
seg4 = cran.vol(:);
seg5 = Atlas.vol(:);
s1=regionprops(m.vol);
s2=regionprops(w.vol);
s3=regionprops(mil.vol);
s4=regionprops(cran.vol);
s5=regionprops(Atlas.vol);
tt=s5.Centroid;

%McIntyre vs. Wu
voxelsnumber1 = sum(seg1);
voxelsnumber2 = sum(seg2);
commonarea = sum(seg1&seg2);

dice_mcintyre_wu = (2*(commonarea)/(voxelsnumber1+voxelsnumber2));
jaccard_mcintyre_wu = jaccard(seg1, seg2);
dist_mcintyre_wu = norm(s1.Centroid-s2(1).Centroid);

%McIntyre vs. Miller
voxelsnumber1 = sum(seg1);
voxelsnumber3 = sum(seg3);
commonarea = sum(seg1&seg3);
dist_mcintyre_miller = norm(s1.Centroid-s3.Centroid);
dice_mcintyre_miller = (2*(commonarea)/(voxelsnumber1+voxelsnumber3));
jaccard_mcintyre_miller = jaccard(seg1, seg3);

%McIntyre vs. CS
voxelsnumber1 = sum(seg1);
voxelsnumber4 = sum(seg4);
commonarea = sum(seg1&seg4);
dist_mcintyre_CS = norm(s1.Centroid-s4.Centroid);
dice_mcintyre_CS = (2*(commonarea)/(voxelsnumber1+voxelsnumber4));
jaccard_mcintyre_CS = jaccard(seg1, seg4);

voxelsnumber1 = sum(seg1);
voxelsnumber5 = sum(seg5);
commonarea = sum(seg1&seg5);
dist_mcintyre_atlas = norm(s1.Centroid-tt);
dice_mcintyre_atlas = (2*(commonarea)/(voxelsnumber1+voxelsnumber5));
jaccard_mcintyre_atlas = jaccard(seg1, seg5);

voxelsnumber2 = sum(seg2);
voxelsnumber3 = sum(seg3);
commonarea = sum(seg2&seg3);
dist_wu_miller = norm(s2(1).Centroid-s3.Centroid);
dice_wu_miller = (2*(commonarea)/(voxelsnumber2+voxelsnumber3));
jaccard_wu_miller = jaccard(seg2,seg3);

voxelsnumber2 = sum(seg2);
voxelsnumber4 = sum(seg4);
commonarea = sum(seg2&seg4);
dist_wu_CS = norm(s2(1).Centroid-s4.Centroid);
dice_wu_CS = (2*(commonarea)/(voxelsnumber2+voxelsnumber4));
jaccard_wu_CS = jaccard(seg2, seg4);

voxelsnumber2 = sum(seg2);
voxelsnumber5 = sum(seg5);
commonarea = sum(seg2&seg5);
dist_wu_atlas = norm(s2(1).Centroid-tt);
dice_wu_atlas = (2*(commonarea)/(voxelsnumber2+voxelsnumber5));
jaccard_wu_atlas = jaccard(seg2, seg5);

voxelsnumber3 = sum(seg3);
voxelsnumber4 = sum(seg4);
commonarea = sum(seg3&seg4);
dist_miller_CS = norm(s3.Centroid-s4.Centroid);
dice_miller_CS = (2*(commonarea)/(voxelsnumber3+voxelsnumber4));
jaccard_miller_CS = jaccard(seg3,seg4);

voxelsnumber3 = sum(seg3);
voxelsnumber5 = sum(seg5);
commonarea = sum(seg3&seg5);
dist_miller_atlas = norm(s3.Centroid-tt);
dice_miller_atlas = (2*(commonarea)/(voxelsnumber3+voxelsnumber5));
jaccard_miller_atlas = jaccard(seg3,seg5);

voxelsnumber4 = sum(seg4);
voxelsnumber5 = sum(seg5);
commonarea = sum(seg4&seg5);
dist_CS_atlas = norm(s4.Centroid-tt);
dice_CS_atlas = (2*(commonarea)/(voxelsnumber4+voxelsnumber5));
jaccard_CS_atlas = jaccard(seg4,seg5);

%A=[1 dice_mw dice_mmil dice_mCS;
%dice_mw 1 dice_wmil dice_wCS;
%dice_mmil dice_wmil 1 dice_milCS;
%dice_mCS dice_wCS dice_milCS 1]

dice(:,count)=[dice_mcintyre_miller; dice_mcintyre_CS; dice_miller_CS; dice_mcintyre_wu; dice_wu_miller; dice_wu_CS; dice_mcintyre_atlas; dice_wu_atlas; dice_miller_atlas; dice_CS_atlas];
JAC(:,count)=[jaccard_mcintyre_miller; jaccard_mcintyre_CS; jaccard_miller_CS; jaccard_mcintyre_wu; jaccard_wu_miller; jaccard_wu_CS; jaccard_mcintyre_atlas; jaccard_wu_atlas; jaccard_miller_atlas; jaccard_CS_atlas];
dist(:,count)=[dist_mcintyre_miller; dist_mcintyre_CS; dist_miller_CS; dist_mcintyre_wu; dist_wu_miller; dist_wu_CS; dist_mcintyre_atlas; dist_wu_atlas; dist_miller_atlas; dist_CS_atlas];
dice_lab(:,count)=[1;2;3;4;5;6;7;8;9;10];
seg_lab(:,count)={'mcintyre_miller'; 'mcintyre_CS'; 'miller_CS'; 'mcintyre_wu'; 'wu_miller'; 'wu_CS'; 'mcintyre_atlas'; 'wu_atlas'; 'miller_atlas'; 'CS_atlas'};
region_lab(:,count)={mROI{count}; mROI{count}; mROI{count}; mROI{count}; mROI{count}; mROI{count}; mROI{count}; mROI{count}; mROI{count}; mROI{count}};

clear m.vol(:);
clear mil.vol(:);
clear CS.vol(:);
clear w.vol
clear Atlas.vol(:);

clear s1.Centroid
clear s2.Centroid
clear s3.Centroid
clear s4.Centroid
clear s5.Centroid

end

D=reshape(dice,[60 1]);
Region_Labels = reshape(region_lab, [60 1]);
Segmentation_Labels = reshape(seg_lab, [60 1]);
J=reshape(JAC,[60 1]);
Dist=reshape(dist,[60 1]);
seg=reshape(seg_lab,[60 1]);
reg=reshape(region_lab,[60 1]);

for k=1:10
    Davg(k,1)=((D(k)+D(k+10))/2);
end
for j=21:30
    Davg2(j,1)=((D(j)+D(j+10))/2);
end
for j=41:50
    Davg3(j,1)=((D(j)+D(j+10))/2);
end

Final_Dice = [Davg; nonzeros(Davg2); nonzeros(Davg3)]

for k=1:10
    Javg(k,1)=((J(k)+J(k+10))/2);
end
for j=21:30
    Javg2(j,1)=((J(j)+J(j+10))/2);
end
for j=41:50
    Javg3(j,1)=((J(j)+J(j+10))/2);
end

Final_Jaccard = [Javg; nonzeros(Javg2); nonzeros(Javg3)]

for k=1:10
    Dist_avg(k,1)=((Dist(k)+Dist(k+10))/2);
end
for j=21:30
    Dist_avg2(j,1)=((Dist(j)+Dist(j+10))/2);
end
for j=41:50
    Dist_avg3(j,1)=((Dist(j)+Dist(j+10))/2);
end

Final_Distance = [Dist_avg; nonzeros(Dist_avg2); nonzeros(Dist_avg3)]

A=[num2cell(D) num2cell(J) num2cell(Dist) Region_Labels Segmentation_Labels]
A=sortrows(A,[5,4])