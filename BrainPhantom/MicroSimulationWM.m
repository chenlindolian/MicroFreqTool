%% Making WM microstructure induced frequency shift for numerical head phantom 
% Combine the WM signal simulated using HCFM with 2 fiber populations to 
% compute the signal frequency f and fit the deltaf,p1,p2,Cf using 
% TE-dependent fitting method,which can be served as the
% ground truth; (f-Cf) served as the frequency shift induced by WM
% microstructure for the input of brain phantom
clc,clear all
addpath('~/Desktop/WM/code/SupportFunction/nifti/')
addpath('~/Desktop/WM/code/Exampledata/ExampleImage/')
addpath('~/Desktop/WM/code/Exampledata/HCFM/HCFM2FiberSus/')
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

%% Calculate the signal by HCFM with two fiber populations

numOfVoxels=length(Index3);

TEs = (1:1:30)*1e-3;

WMsingal_normalized = zeros(numOfVoxels, length(TEs));

load('FittedSignal_batch1.mat')

WMsingal_normalized(1:120000,:) = singal_normalized(1:120000,:);

load('FittedSignal_batch2.mat')

WMsingal_normalized(120001:240000,:) = singal_normalized(120001:240000,:) ;

load('FittedSignal_batch3.mat')

WMsingal_normalized(240001:360000,:) = singal_normalized(240001:360000,:) ;

load('FittedSignal_batch4.mat')

WMsingal_normalized(360001:480000,:) = singal_normalized(360001:480000,:) ;

load('FittedSignal_batch5.mat')

WMsingal_normalized(480001:numOfVoxels,:) = singal_normalized(480001:numOfVoxels,:);

WMsingal_normalized(isnan(WMsingal_normalized)==1) = 0;

magn_signal = abs(WMsingal_normalized);

phase_signal = zeros(numOfVoxels, length(TEs));

for ivoxel = 1: numOfVoxels

phase_signal(ivoxel,:) = phase(WMsingal_normalized(ivoxel,:));

end

freq_signal = phase_signal ./(2*pi*TEs);

BrainMicroFreq = zeros([size(maskWM),length(TEs)]);

for ivoxel= 1: numOfVoxels 

BrainMicroFreq(Index1(ivoxel),Index2(ivoxel),Index3(ivoxel),:) = freq_signal(ivoxel,:);

end

%% calculate curve for frequency difference map
x0 = [1 100 0 0];
ub = [8 300];
delta_freq = zeros(numOfVoxels,1);
p1 = zeros(numOfVoxels,1);
p2= zeros(numOfVoxels,1);
Cf = zeros(numOfVoxels,1);
freq_sim = zeros(numOfVoxels,length(TEs));

for ivoxel= 1:numOfVoxels
   x = lsqcurvefit(@Freqmodel, x0, TEs, freq_signal(ivoxel,:),[],ub);
   freq_sim(ivoxel,:) = -x(1).*(exp(-x(2)*TEs)./(x(3)*(exp(-x(2)*TEs)-1)+1))+x(4);
   delta_freq(ivoxel) = x(1);
   p1(ivoxel)=x(2);
   p2(ivoxel)=x(3);
   Cf(ivoxel) = x(4);
end

BrainDifferenceFreq = zeros(size(maskWM));
curvep1 = zeros(size(maskWM));
curvep2 = zeros(size(maskWM));
Cfmap = zeros(size(maskWM));
p1(isnan(WMsingal_normalized)==1) = 0;
p2(isnan(WMsingal_normalized)==1) = 0;

for ivoxel= 1:numOfVoxels
    BrainDifferenceFreq(Index1(ivoxel),Index2(ivoxel),Index3(ivoxel)) =  delta_freq(ivoxel);
    curvep1(Index1(ivoxel),Index2(ivoxel),Index3(ivoxel)) =  p1(ivoxel);
    curvep2(Index1(ivoxel),Index2(ivoxel),Index3(ivoxel)) =  p2(ivoxel);
    Cfmap(Index1(ivoxel),Index2(ivoxel),Index3(ivoxel)) =  Cf(ivoxel);
end
save('BrainSimMicroFreq.mat','BrainDifferenceFreq','BrainMicroFreq','delta_freq','p1','p2','Cf','curvep1','curvep2','Cfmap');

 function ydata = Freqmodel(x,iTE)
  ydata = -x(1).*(exp(-x(2)*iTE)./(x(3)*(exp(-x(2)*iTE)-1)+1))+x(4);
 end
 