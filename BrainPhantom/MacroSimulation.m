%% Making bulk susceptibility induced frequency shift for numerical head phantom 
% using T1 data and diffusion data(mutli shell diffusion data/HCP) acquired on 7T/3T
% CSF/GM/WM data using FSL FAST segmentation result
% deep GM using ANTs co-registration segmentation result
% bulk anisotropic susceptibility of WM voxel were modulate by OD and
% transformed to common subject frame using DTI eigenvector
% assign isotropic susceptibility, S0, R1 based on Table S1

clc,clear all
addpath('~/Desktop/WM/code/SupportFunction/nifti/')
addpath('~/Desktop/WM/code/Exampledata/ExampleImage/')

%% Load image
nii = load_untouch_nii('CSF.nii.gz');
EveCSF = nii.img;

nii = load_untouch_nii('GM.nii.gz');
EveGM = nii.img;

nii = load_untouch_nii('WM.nii.gz');
EveWM = nii.img;

nii = load_untouch_nii('R2star.nii.gz');
EveR2star = nii.img;

nii = load_untouch_nii('dt_vector1.nii.gz');
Eigenvec1 = nii.img;

nii = load_untouch_nii('dt_vector2.nii.gz');
Eigenvec2 = nii.img;

nii = load_untouch_nii('dt_vector3.nii.gz');
Eigenvec3 = nii.img;

nii = load_untouch_nii('T1atlas.img');
T1atlas = nii.img;
T1atlas = flip(T1atlas,2);
T1atlas(:,:,1) = T1atlas(:,:,2);

nii = load_untouch_nii('QSMatlas');
dGMatlas = nii.img;

nii = load_untouch_nii('vesselMask.img');
vesselatlas = nii.img;
vesselatlas = flip(vesselatlas,2);

odi_struct = load_untouch_nii('PhantomNODDI_odi.nii'); % orientation dispersion
odi = flip(odi_struct.img,2);

BrainMask = (EveGM>=0.5) | (EveWM>=0.5) | (EveCSF>=0.5);
for sliceii = 1:126
    BrainMask(:,:,sliceii) = imfill(BrainMask(:,:,sliceii), 'holes');
end

dims = size(BrainMask);
%% brain phantom isotropic susceptibility
chi_i = zeros([dims,3,3]);
chi_i0 = zeros([dims,3,3]);
Eve_chi0 = zeros(dims);
Eve_chi1 = zeros(dims);
Eve_GREMag0 = zeros(dims);
Eve_R1 = zeros(dims);
%% part 1
% For background tissue: Air/Bone/Fat, Note that Fat has chemical shift and bulk susceptibilty 
% 250: Air -->   0.36 ppm          (Ref, Schenck, 1996, Med Phys)
% 249: Fat -->   -8.5 ppm         -9.05 + 0.4, 0.3 ~ 0.5 ppm relative to water (Ref. Dimov 2015, Wei 2017, )
% 248, 254: Bone -->  -11.56 ppm        -9.05 - 4*pi*(0.19 ~ 0.21) (Ref. Hopkins and Wehrli, 1997, MRM)
% Water/CSF --> -9.05 ppm   
% Veins --> -9.05 + 0.3 pp,         not included in the model yet
roiName_1 = {'Air'; 'Fat'; 'Bone'};

roiIndex = [250, -1; ...    % Air
            249, -1; ...    % Fat
            248, 254];      % Bone
        
roiValue_chi = [0.36, -8.5, -11.56];  % in ppm
roiValue_Mag = [0, 0.4, 0.001];        % normalized magnitude
roiValue_R1 = [0 0.7 0.7];
 
for roiIndexii = 1:length(roiName_1)
    
    maskROI = (T1atlas == roiIndex(roiIndexii, 1) | T1atlas == roiIndex(roiIndexii, 2) );
    
    Eve_chi0(maskROI) = roiValue_chi(roiIndexii); 
    
    Eve_GREMag0(maskROI) = roiValue_Mag(roiIndexii);
    
    Eve_R1(maskROI) = roiValue_R1(roiIndexii);
end
              
Eve_GREMag0(BrainMask) = 0;
Eve_R1(BrainMask) = 0; 
Eve_chi0(BrainMask) = -9.05;  

%% part 2
% set WM/GM/CSF susceptibility, WM, -0.03; GM: 0.04; relative to CSF of 0, in ppm
% set WM/GM/CSF magnitude, WM = 0.95, GM: 1, CSF: 0.9, in a.u.

% considering partial volume too
Eve_chi1 = Eve_chi1 + EveWM.*-0.03;  % WM
Eve_chi1 = Eve_chi1 + EveGM.*0.04;  % GM
Eve_chi1 = Eve_chi1 + EveCSF.*0;  % CSF

Eve_GREMag0 = Eve_GREMag0 + EveWM.*0.7;  % WM
Eve_GREMag0 = Eve_GREMag0 + EveGM.*0.85;     % GM
Eve_GREMag0 = Eve_GREMag0 + EveCSF.*1;  % CSF

Eve_R1 = Eve_R1+ EveWM.*0.9;  % WM
Eve_R1 = Eve_R1 + EveGM.*0.6;     % GM
Eve_R1 = Eve_R1 + EveCSF.*0.25;  % CSF

% small ROIs
roiName = {'CN'; 'PT'; 'GP'; 'Thal';'RN';'SN'};   % Inside of Brain           
roiNameShort = roiName;

roiIndex = [77,78;... %CN
            79,80;... %PT
            81,82;... %GP
            83,84;... %Tha
            91,92;...%RN
            93,94];  %SN

roiValue_chi = [0.08, 0.07, 0.18, 0.05, 0.13, 0.15];  % in ppm
roiValue_Mag = [0.8, 0.75, 0.7, 0.75, 0.75, 0.75];         % normalized magnitude
roiValue_R1 = [0.6, 0.6, 0.8, 0.6, 0.8, 0.8];         % R1 in S-1
vessel_chi = 0.3; 
vessel_Mag = 0.8;
vessel_R1 = 0.5;
maskdGM = zeros(dims);
for roiIndexii = 1:length(roiName)
    
    maskROI = (dGMatlas  == roiIndex(roiIndexii, 1) | dGMatlas  == roiIndex(roiIndexii, 2) );
    
    Eve_chi1(maskROI) = roiValue_chi(roiIndexii);       % relative to CSF    
    Eve_GREMag0(maskROI) = roiValue_Mag(roiIndexii);

    Eve_R1(maskROI) = roiValue_R1(roiIndexii);

    maskdGM = maskdGM | maskROI;
end

vesselMask = (vesselatlas == 1);

Eve_chi1(vesselMask) = vessel_chi;
Eve_GREMag0(vesselMask) = vessel_Mag;
Eve_R1(vesselMask) = vessel_R1;

Eve_chi = Eve_chi0 + Eve_chi1;  % total chi_i
%% add noise
%% 
Eve_chi1 = (Eve_chi1 + 0.001*randn(size(Eve_chi1))).*BrainMask;             % add physiological noise in brain

chi_i(:,:,:,1,1) = Eve_chi1;

chi_i(:,:,:,2,2) = Eve_chi1;

chi_i(:,:,:,3,3) = Eve_chi1;

maskWM = BrainMask & (EveWM >= 0.5) & ~maskdGM;
maskCSF = BrainMask & (EveCSF == 1);
%% part 3, susceptibility transformation based on eigenvectors from diffusion image
chi_para = 0.014;
chi_perp = -0.009;
myelin_chia = [chi_para, 0, 0;0 chi_perp 0;0 0 chi_perp];

chi_a = zeros([dims,3,3]);

[Index1,Index2,Index3] = ind2sub(size(maskWM),find(maskWM > 0)); 

for i = 1:length(Index3)
    eigenvec1 = reshape(squeeze(Eigenvec1(Index1(i),Index2(i),Index3(i),:)),[1,3]); 
    eigenvec2 = reshape(squeeze(Eigenvec2(Index1(i),Index2(i),Index3(i),:)),[1,3]); 
    eigenvec3 = reshape(squeeze(Eigenvec3(Index1(i),Index2(i),Index3(i),:)),[1,3]); 
    vec = [eigenvec1;eigenvec2;eigenvec3];

    %% Modualte susceptibility anisotropy by orientation dispersion index
    OD = odi(Index1(i),Index2(i),Index3(i));
   
    chi_a(Index1(i),Index2(i),Index3(i),:,:)=vec'*myelin_chia*vec*(1-OD);
end

Simchi = chi_i + chi_a; % tissue chi

chi_i0(:,:,:,1,1) = Eve_chi1 + Eve_chi0;

chi_i0(:,:,:,2,2) = Eve_chi1 + Eve_chi0;

chi_i0(:,:,:,3,3) = Eve_chi1 + Eve_chi0;

Simchitotal = chi_i0 + chi_a; % total chi

EveR2star = EveR2star.*BrainMask; % R2star

save('BrainSimMacroSus.mat','Eve_GREMag0','EveR2star','Eve_R1','Simchi','Simchitotal','chi_i','chi_a','maskWM','Eve_chi1','BrainMask','maskCSF','vesselMask');



