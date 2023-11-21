% Simulate signal in WM using HCFM with 2 fiber populations
% H0 = [sind(theta)*cosd(phi) sind(theta)*sind(phi) cosd(theta)]';
% cross angle between two fiber populations derived from peak maps
clc,clear all
addpath('~/Desktop/WM/code/SupportFunction/nifti/')
addpath('~/Desktop/WM/code/Exampledata/ExampleImage/')
addpath('~/Desktop/WM/code/Exampledata/HCFM/HCFM2FiberSus/')

%% Create WM mask
nii = load_untouch_nii('CSF.nii.gz');
EveCSF = nii.img;

nii = load_untouch_nii('GM.nii.gz');
EveGM = nii.img;

nii = load_untouch_nii('WM.nii.gz');
EveWM = nii.img;

nii = load_untouch_nii('QSMAtlas.nii');
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
%% simulate signal voxel by voxel
Params.B0 = 7;              % Tesla
Params.gamma = 42.577e6;      % Hz/T

% myelin (required: T2, xi, xa)
Params.myelin.Mag= 0.7; 
Params.myelin.T2 = 8*1e-3;
Params.myelin.T1 = 500*1e-3;
Params.myelin.proton_density= 1; 

Params.myelin.chii = -0.10;  % myelin isotropic susceptibility (ppm)
Params.myelin.chia = -0.10;  % myelin anisotropic susceptibility (ppm)
% Params.myelin.chia = 0;

% intra axonal (required: T2)
Params.intra_axonal.Mag = 1;
Params.intra_axonal.T2 = 36*1e-3;
Params.intra_axonal.T1 = 1.5;
Params.intra_axonal.proton_density= 1; 
Params.intra_axonal.chii= 0; 

% extra axonal (required: T2)
Params.extra_axonal.Mag = 1;
Params.extra_axonal.T2 = 36*1e-3;
Params.extra_axonal.T1 = 1.5;
Params.extra_axonal.proton_density = 1; 
Params.extra_axonal.chii = 0; 

Params.E = 0.02;

Params.TEs = (1:1:30)*1e-3;

BrainMicro = zeros([size(maskWM),length(Params.TEs)]);
pad_size = [50 50 50];

numOfVoxels=length(Index3);

outputfile = 'FittedSignal_batch.mat';
load('angleB0.mat')

if exist(outputfile, 'file')
    output = load(outputfile);
    if isfield(output, 'split_end')
        % previously run and need to be restarted
        split_end = output.split_end;
        singal_normalized = output.singal_normalized;
        current_split_start = split_end + 1;
        fprintf('Resume an interrupted run from %i\n', current_split_start);
    else
        % completed
        fprintf('An output file of the same name detected.\n');
        fprintf('Choose a different output file name.\n');
        return;
    end
else
    % if this is the first run
    current_split_start = 1;
end
delete(gcp('nocreate'))
c = parcluster;
parpool(c)

logfile = fullfile(fileparts(mfilename('fullpath')), 'logd_batch1.txt');
logtxt = ['poolsize of current pool:', num2str(c.NumWorkers)];
writelog(logfile, [logtxt, '\n'], 'a');


if current_split_start == 1 % set up the fitting parameter variables if it is the first run
    singal_normalized = zeros(numOfVoxels, length(Params.TEs));
end

% set up the PARFOR Progress Monitor
progressStepSize = 2000;

tic

fprintf('%i of voxels to fit\n', numOfVoxels-current_split_start+1);

%==========================================================================
% Process all pixels
%==========================================================================
for split_start=current_split_start:progressStepSize:current_split_end
     % set up the split end
    split_end = split_start + progressStepSize - 1;
    if split_end > numOfVoxels
        split_end = numOfVoxels;
    end
    
    logfile = fullfile(fileparts(mfilename('fullpath')), 'logd_batch1.txt');
    logtxt = ['Run voxels number:', num2str(split_start),' ', datestr(datetime('now'))];
    writelog(logfile, [logtxt, '\n'], 'a'); 
  
  % fit the split
    parfor i=split_start:split_end
           angleCross = angleCrossMap(Index1(i),Index2(i),Index3(i));
           model_angle = round(angleCross/5)*5; %find round model
           mat = load(['317Fiber2BundleCrossAngle',num2str(model_angle),'.mat']);
           
           Model = [];
           Model(i).model = mat.Phantom.model;
           Model(i).mask = mat.Phantom.mask;
           Model(i).FOV = mat.Phantom.FOV;
           Model(i).chi = mat.total_chi;
          
           theta_sub = theta(Index1(i),Index2(i),Index3(i));
           phi_sub = fai(Index1(i),Index2(i),Index3(i));
           H0 = [sind(theta_sub)*cosd(phi_sub) sind(theta_sub)*sind(phi_sub) cosd(theta_sub)]';
           
           Model(i).chi = padarray(Model(i).chi,pad_size,0,'both');
           Model(i).model = padarray(Model(i).model,pad_size,1,'both');
          
           deltaB = SimulateFieldFrom3DTensor(Model(i).chi, Model(i), Params, H0);

           deltaB = deltaB(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
           Model(i).model = Model(i).model(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));

           signal = simulateSignalFrom3DField(Model(i), deltaB, Params);
    
           singal_normalized(i,:) = signal.total_normalized;
  
    end
    save(outputfile, 'split_end','singal_normalized');
    
end

%==========================================================================
delete(pool);



         
