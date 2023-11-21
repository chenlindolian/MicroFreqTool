% Fiber ODF (oritation distibution function) peak maps were reconstructed
% using multi-shell multi-tissue CSD
% preprocessed_mean_bzero_b02g0GenericAffine.mat is the transform matrix
% derived from rigid coregistration of b0 image to GRE space
% Calculate angle for H0 direction using peak maps and cross angle between
% two fiber populations
% H0 = [sind(theta)*cosd(phi) sind(theta)*sind(phi) cosd(theta)]';

clc,clear all
addpath('~/Desktop/WM/code/SupportFunction/nifti/')
addpath('~/Desktop/WM/code/Exampledata/ExampleImage')

%% B0 used in the CSD algorithm
H0 = [0,0,1]'; % H0 in lab frame

Params.TAng = [1,0,0;0,1,0;0,0,1];

load('preprocessed_mean_bzero_b02g0GenericAffine.mat');

Params.TAng = reshape(AffineTransform_double_3_3(1:9),[3 3])*Params.TAng;

H0sub = Params.TAng'*H0;     % H0 in image/subject space

nii = load_untouch_nii('peak1.nii.gz');
Peak1 = nii.img;

nii = load_untouch_nii('peak1.nii.gz');
Peak2 = nii.img;

nii = load_untouch_nii('CSF.nii.gz');
EveCSF = nii.img;

nii = load_untouch_nii('GM.nii.gz');
EveGM = nii.img;

nii = load_untouch_nii('WM.nii.gz');
EveWM = nii.img;

nii = load_untouch_nii('Atlas.nii');
dGMatlas = nii.img;

BrainMask = (EveGM>=0.5) | (EveWM>=0.5) | (EveCSF>=0.5);
for sliceii = 1:126
    BrainMask(:,:,sliceii) = imfill(BrainMask(:,:,sliceii), 'holes');
end

maskdGM = zeros(size(dGMatlas));

roiIndex = [77,78;... %CN
            79,80;... %PT
            81,82;... %GP
            83,84;... %Tha
            91,92;...%RN
            93,94];  %SN

for roiIndexii = 1:length(roiIndex)
    
    maskROI = (dGMatlas  == roiIndex(roiIndexii, 1) | dGMatlas  == roiIndex(roiIndexii, 2) );
   
    maskdGM = maskdGM | maskROI;
end

maskWM = BrainMask & (EveWM >= 0.5) & ~maskdGM;

[Index1,Index2,Index3] = ind2sub(size(maskWM),find(maskWM > 0)); 

angleCrossMap = zeros(size(maskWM));
theta = zeros(size(maskWM));
fai = zeros(size(maskWM));

for i = 1:length(Index3)
    vec1 = squeeze(Peak1(Index1(i),Index2(i),Index3(i),:)); 
    vec2 = squeeze(Peak2(Index1(i),Index2(i),Index3(i),:)); 
    anglePeak1 = acosd(dot(vec1,H0sub)/(norm(vec1)*norm(H0sub)));
    anglePeak2 = acosd(dot(vec2,H0sub)/(norm(vec2)*norm(H0sub)));
    angleCross = acosd(dot(vec1,vec2)/(norm(vec1)*norm(vec2)));
    length1 = norm(vec1);
    length2 = norm(vec2);
    if length2 > 0.3 * length1
       angleCrossMap(Index1(i),Index2(i),Index3(i)) = angleCross;
       if angleCross > 90
          angleCrossMap(Index1(i),Index2(i),Index3(i)) = 180-angleCross;
       end
    end

    theta(Index1(i),Index2(i),Index3(i)) = anglePeak1;
    fai(Index1(i),Index2(i),Index3(i)) = real(asind((cosd(anglePeak2)-cosd(angleCross)*cosd(anglePeak1))/(sind(angleCross)*sind(anglePeak1))));

end

save('angleB0','theta','fai','angleCrossMap');